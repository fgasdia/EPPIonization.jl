module EPPIonization

using Interpolations
using HDF5

export ionprofile

datadir(f) = joinpath(@__DIR__, "..", "data", f)

# For type stability
function readlut()
    lut = h5read(datadir("lookupTable.h5"), "/")
    en = lut["en"]::Vector{Float64}
    pa = lut["pa"]::Vector{Int64}
    ion = lut["ion"]::Array{Float64,3}
    masden = lut["masden"]::Vector{Float64}
    alt = lut["alt"]::Vector{Int64}
    return en, pa, ion, masden, alt
end

"""
    ionprofile(altitude, energy, energydis, pitchangle, pitchdis, massdensity)

Compute ionization profile as a function of `altitude` in km.

# Arguments

    - `altitude` [km]: altitude vector
    - `energy` [keV]: electron energy between ~3 keV and ~33 MeV
    - `energydis`: energy distribution of precipitating electrons
    - `pitchangle` [deg]: pitch angle between 0 and 90°
    - `pitchdis`: pitch angle distribution
    - `massdensity` [g/cm²]: mass density in new atmosphere defined from 0:500 km.

# Example

```julia
energy = 10.^(4:0.01:7);
energydis = exp.(-energy/1e5);  # f(E) ∝ exp(-E/100 keV)
pitchangle = 0:90;
pitchdis = ones(length(pitchangle));
ion = ionprofile(alt, energy, energydis, pitchangle, pitchdis, massdensity)
```
"""
function ionprofile(altitude, energy, energydis, pitchangle, pitchdis, massdensity)
    en, pa, ion, masden, alt = readlut()

    # Mask values up to maximum of `altitude`
    mask = alt .<= maximum(altitude)
    alt = alt[mask]
    masden = masden[mask]
    ion = ion[:,:,mask]
    
    # calculate energy bins
    edge = Vector{Float64}(undef, length(en)+1)
    edge[1] = en[1] - (en[2] - en[1])/2
    @views edge[2:length(en)] .= (en[1:end-1] .+ en[2:end])./2
    edge[end] = en[end] + (en[end] - en[end-1])/2
    den = diff(edge)

    # convert the lookup table to a new atmosphere
    ionint = Array{Float64,3}(undef, length(en), length(pa), length(massdensity))
    for j = 1:length(pa) - 1
        for i = 1:length(en)
            @inbounds ionint[i,j,:] = lookupconvert(masden, view(ion,i,j,:), massdensity, alt)
        end
    end
    ionint[:,end,:] .= 0

    # interpolate the input energy and pitch angle distributions
    enind = findall(x->minimum(energy) <= x <= maximum(energy), en)

    # interpolate in energy
    # enint=10.^interp1(log10(energy),log10(energydis),log10(en(enind)),'linear','extrap');
    en_itp = LinearInterpolation(energy, energydis, extrapolation_bc=Line())
    enint = en_itp(view(en,enind))

    # interpolate in pitch angle
    pa_itp = LinearInterpolation(pitchangle, pitchdis)
    paint = pa_itp(pa)

    # normalize in energy 
    endis = enint./sum(enint.*view(den,enind))

    # normalize in pitch angle
    padis = paint./sum(paint)

    # calculate the ionization profile, sum of ionization contribution from each energy and pitch angle component
    ionen = zeros(length(enind), length(alt))
    ionnew = zeros(length(alt))
    for i = 1:length(enind)
        for h = 1:length(pa)
            @views ionen[i,:] .+= padis[h].*ionint[enind[i],h,:]
        end
        @views ionnew .+= (endis[i]*den[enind[i]]).*ionen[i,:]
    end

    # extrapolate to altitudes above 500 km
    # Because of the interpolation method, the last element of ionnew is 0. 
    # Thus can only do interpolation in log space between 1 and end-1.
    # Extrapolation in log space to altitudes above 500 km may sometimes look werid. 
    # Extra caution!
    @views ion_itp = LinearInterpolation(alt[1:end-1], log.(ionnew[1:end-1]), extrapolation_bc=Line())
    ionpro = exp.(ion_itp(altitude))
    replace!(ionpro, NaN=>0)

    return ionpro
end

function lookupconvert(denold1, ionold1, dennew1, alt)
    if count(>(0), ionold1) < 2 || sum(ionold1) == 1
        return ionold1
    end
    length(dennew1) == length(alt) || throw(ArgumentError("`dennew1` must be length $(length(alt))"))

    # TODO: we can precompute some of the values in this function

    # calculate ionization rate at half grid cells
    alt2 = (alt .+ (alt .+ 1))./2

    do_itp = LinearInterpolation(alt, log.(denold1), extrapolation_bc=Line())
    io_itp = LinearInterpolation(alt, log.(ionold1), extrapolation_bc=Line())
    dn_itp = LinearInterpolation(alt, log.(dennew1), extrapolation_bc=Line())
    denold = exp.(do_itp(alt2))
    ionold = exp.(io_itp(alt2))
    dennew = exp.(dn_itp(alt2))

    # sort the ionization production and mass density
    # from 500 km to the lowest altitude of ionzation production
    reverse!(denold)
    reverse!(dennew)
    reverse!(ionold)

    # find out the lowest altitude of energy deposition in reference and new atmosphere
    # `ionold1` last element is `0`.
    # After interpolation, `ionold` ends with two `0`s and after `reverse!`, the first two
    # entries are zeros. Therefore, we find the next zero after element 2
    minaltold = findnext(isequal(0), ionold, 3) - 1
    denold[minaltold+1:end] .= 0

    sum_denold = sum(denold)
    minaltnew = findfirst(>=(sum_denold), cumsum(dennew))
    if sum(view(dennew, 1:minaltnew)) - sum_denold > sum_denold - sum(view(dennew, 1:minaltnew-1))
        minaltnew -= 1
    end
    dennew[minaltnew+1:end] .= 0

    # cumulatively sum mass density and ionization rate
    denold_csum = round.(cumsum(denold)/sum_denold, digits=15)  # there's a problem if 0.9999999999999999
    ionold_csum = cumsum(ionold)
    dennew_csum = round.(cumsum(dennew)/sum(dennew), digits=15)
    ionnew_csum = zeros(length(dennew))
    ionnew = zeros(length(dennew))

    # interpolate the cumulative sum of ionization rate from 500 km to the lowest altitude
    intind = findfirst(>(denold_csum[3]), dennew_csum)
    if isnothing(intind)
        # pass
    elseif intind > 1
        min_itp = LinearInterpolation(log.(view(denold_csum,3:minaltold)),
                                      log.(view(ionold_csum,3:minaltold)), extrapolation_bc=Line())
        ionnew_csum[intind:minaltnew] = exp.(min_itp(log.(view(dennew_csum, intind:minaltnew))))

        uniqueval = unique(ionold_csum)
        uniind = findfirst(isequal(uniqueval[end]), ionold_csum)
        uni_itp = LinearInterpolation(view(denold_csum,3:uniind),
                                      view(ionold_csum,3:uniind), extrapolation_bc=Line())
        ionnew_csum[1:intind] = uni_itp(view(dennew_csum,1:intind))
    elseif intind == 1
        min_itp = LinearInterpolation(log.(view(denold_csum,3:minaltold)),
                                      log.(view(ionold_csum,3:minaltold)), extrapolation_bc=Line())
        ionnew_csum[intind:minaltnew] = exp.(min_itp(log.(view(dennew_csum, intind:minaltnew))))
    end

    # differentiate the interpolation results and calculate ionization rate at each altitude
    @views ionnew[2:minaltnew] .= ionnew_csum[2:minaltnew] .- ionnew_csum[1:minaltnew-1]
    replace!(x->x < 0 ? zero(x) : x, ionnew)

    return reverse!(ionnew)  # back to alt order
end

end # module
