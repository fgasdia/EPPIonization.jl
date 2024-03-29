"""
    EPPIonization

This package is based upon

    Xu, Wei, Marshall, Robert A, Tyssoy, Hilde Nesse, & Fang, Xiaohua. (2020).
    A Generalized Method for Calculating Atmospheric Ionization by Energetic Electron
    Precipitation. https://doi.org/10.5281/zenodo.3945306

The original Matlab source code is available at https://zenodo.org/record/3945306 under
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode)
license.
"""
module EPPIonization

using Dates
using Interpolations, HDF5, TypedTables
using SpaceIndices, SatelliteToolboxBase, SatelliteToolboxAtmosphericModels
using GPI, FaradayInternationalReferenceIonosphere, LMPTools
import FaradayInternationalReferenceIonosphere as FIRI

export ionizationprofile, neutralprofiles, chargeprofiles, EnergeticElectrons

const SPECIES = GPI.SPECIES
const INITIALIZED = Ref{Bool}(false)

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
    EnergeticElectrons

Describe the energy and pitch angle distribution of energetic precipitating electrons.

# Example
```julia
energy = 90e3:1e4:2.2e6  # 90 keV to 2.2 MeV every 10 keV
energydis = exp.(-energy/2e5)  # f(E) = exp(-E/β) where β ranges from 100 to 300 keV (Whittaker, 2013)
pitchangle = 0:90
pitchangledis = ones(length(pitchangle))

ee = EnergeticElectrons(energy, energydis, pitchangle, pitchangledis)
```
"""
struct EnergeticElectrons{T1,T2}
    energy::T1
    energydis::Vector{Float64}
    pitchangle::T2
    pitchangledis::Vector{Float64}
end

"""
    set_initialized!(::Bool)

Force the global variable `INITIALIZED`, corresponding to whether or not `SpaceIndices.init`
has been run, to `v::Bool`.
"""
function set_initialized!(v::Bool)
    INITIALIZED[] = v
end

"""
    ionizationprofile(ee, massdensity, altitude) → ionizationprofile
    ionizationprofile(energy, energydis, pitchangle, pitchdis, massdensity, altitude) → ionizationprofile

Compute ionization rate profile in pairs/el/cm as a function of `altitude` in km.

To compute the ionization in pairs/cm³/s, multiply the output of this function by the
source precipitating electron flux in el/cm²/s.

# Arguments

    - `ee`: `EnergeticElectrons`
    - `altitude` [km]: altitude vector
    - `energy` [eV]: electron energy between ~3 keV and ~33 MeV
    - `energydis`: energy distribution of precipitating electrons
    - `pitchangle` [deg]: pitch angle between 0 and 90°
    - `pitchdis`: pitch angle distribution
    - `massdensity` [g/cm³/1000]: mass density in new atmosphere in g/cm³ divided by 1000
        defined at kilometer increments from 0 to the maximum of `altitude`

# Example

```julia
alt = 0:100;
energy = 90e3:1e4:2.2e6  # 90 keV to 2.2 MeV every 10 keV
energydis = exp.(-energy/1e5);  # e.g.: f(E) ∝ exp(-E/100 keV)
pitchangle = 0:90;
pitchdis = ones(length(pitchangle));
ion = ionizationprofile(energy, energydis, pitchangle, pitchdis, massdensity, alt);
```
"""
function ionizationprofile(ee, massdensity, altitude)
    ionizationprofile(ee.energy, ee.energydis, ee.pitchangle, ee.pitchangledis, massdensity, altitude)
end

function ionizationprofile(energy, energydis, pitchangle, pitchdis, massdensity, altitude)
    maximum(altitude) > 500 && throw(ArgumentError("`altitude` should not exceed 500 km"))

    en, pa, ion, masden, alt = readlut()

    # Mask values up to maximum of `altitude`
    mask = alt .<= maximum(altitude)
    alt = alt[mask]
    
    # calculate energy bins
    edge = Vector{Float64}(undef, length(en)+1)
    edge[1] = en[1] - (en[2] - en[1])/2
    @views edge[2:length(en)] .= (en[1:end-1] .+ en[2:end])./2
    edge[end] = en[end] + (en[end] - en[end-1])/2
    den = diff(edge)

    # initialize for lookupconvert
    alt2 = alt .+ 0.5
    @assert alt2 == (alt .+ (alt .+ 1))./2
    do_itp = linear_interpolation(alt, log.(view(masden, mask)), extrapolation_bc=Line())
    denold = exp.(do_itp(alt2))
    reverse!(denold)

    # convert the lookup table to a new atmosphere
    ionint = zeros(length(en), length(pa), length(massdensity))
    for j = 1:length(pa) - 1
        for i = 1:length(en)
            lookupconvert!(view(ionint,i,j,:), denold, view(ion,i,j,mask), massdensity, alt, alt2)
        end
    end

    # interpolate the input energy and pitch angle distributions
    enind = findall(x->minimum(energy) <= x <= maximum(energy), en)

    # interpolate in energy
    en_itp = linear_interpolation(energy, energydis, extrapolation_bc=Line())
    enint = en_itp(view(en,enind))

    # interpolate in pitch angle
    pa_itp = linear_interpolation(pitchangle, pitchdis)
    paint = pa_itp(pa)

    # normalize in energy 
    endis = enint./sum(enint.*view(den, enind))

    # normalize in pitch angle
    padis = paint./sum(paint)

    # calculate the ionization profile, sum of ionization contribution from each energy and pitch angle component
    ionen = zeros(length(enind), length(alt))
    ionnew = zeros(length(alt))
    @views for i = 1:length(enind)
        for h = 1:length(pa)
            ionen[i,:] .+= padis[h].*ionint[enind[i],h,:]
        end
        ionnew .+= (endis[i]*den[enind[i]]).*ionen[i,:]
    end

    @views ion_itp = linear_interpolation(alt[1:end-1], log.(ionnew[1:end-1]), extrapolation_bc=Line())
    ionpro = exp.(ion_itp(altitude))
    replace!(ionpro, NaN=>0)

    return ionpro
end

function lookupconvert!(ionnew, denold, ionold1, dennew1, alt, alt2)
    if count(>(0), ionold1) < 2 || sum(ionold1) == 1
        return copy(ionold1)
    end
    length(dennew1) == length(alt) || throw(ArgumentError("`dennew1` must be length $(length(alt))"))

    denold = copy(denold)  # it's mutated below

    # calculate ionization rate at half grid cells
    io_itp = linear_interpolation(alt, log.(ionold1), extrapolation_bc=Line())
    dn_itp = linear_interpolation(alt, log.(dennew1), extrapolation_bc=Line())
    ionold = exp.(io_itp(alt2))
    dennew = exp.(dn_itp(alt2))

    # sort the ionization production and mass density
    # from 500 km to the lowest altitude of ionzation production
    reverse!(ionold)
    reverse!(dennew)

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
    # as of Julia v1.6.0, cumsum!(x, x) is safe and does what would be expected of cumsum!(x)
    cumsum!(denold, denold)
    denold ./= sum_denold  # denold_csum = round.(cumsum(denold)./sum_denold, digits=15)
    cumsum!(ionold, ionold)  # ionold_cumsum = cumsum(ionold)
    sum_dennew = sum(dennew) # dennew_csum = round.(cumsum(dennew)./sum(dennew), digits=15)
    cumsum!(dennew, dennew)
    dennew ./= sum_dennew

    ionnew_csum = zeros(length(dennew))

    # interpolate the cumulative sum of ionization rate from 500 km to the lowest altitude
    intind = findfirst(>(denold[3]), dennew)
    if !isnothing(intind)
        # intind == 1
        min_itp = linear_interpolation(log.(view(denold, 3:minaltold)),
                                      log.(view(ionold, 3:minaltold)),
                                      extrapolation_bc=Line())
        ionnew_csum[intind:minaltnew] .= exp.(min_itp(log.(view(dennew, intind:minaltnew))))

        if intind > 1
            uniqueval = unique(ionold)
            uniind = findfirst(isequal(uniqueval[end]), ionold)
            uni_itp = linear_interpolation(view(denold, 3:uniind),
                                          view(ionold, 3:uniind),
                                          extrapolation_bc=Line())
            
            ionnew_csum[1:intind] .= uni_itp(view(dennew, 1:intind))
        end
    end

    # differentiate the interpolation results and calculate ionization rate at each altitude
    ionnew[2:minaltnew] .= view(ionnew_csum, 2:minaltnew) .- view(ionnew_csum, 1:minaltnew-1)
    replace!(x->x < 0 ? zero(x) : x, ionnew)

    return reverse!(ionnew)  # back to alt order
end

"""
    massdensity(p, z)
    massdensity(p::GPI.Profiles, z)

Return total neutral atmosphere mass density in g/cm³ at altitude `z` in km.

If `p` is not a `GPI.Profiles`, it will be computed.

`z` must be kilometer intervals.
"""
function massdensity(p::GPI.Profiles, z)
    # mass = number density [#/m³] * molecular mass [g/mol] * [1 mol/#] * [m³/cm³]
    NA = 6.0221e23
    mN2 = GPI.getspecies(p, :N2, z)*28.014
    mO2 = GPI.getspecies(p, :O2, z)*31.999
    mO = GPI.getspecies(p, :O, z)*15.999
    m = (mN2 + mO2 + mO)/NA*1e-6
    
    return m
end

function massdensity(p, z)
    prof = GPI.Profiles(p)
    massdensity(prof, z)
end

"""
    neutralprofiles(lat, lon, z, dt::DateTime) → profiles_table

Return `Table` of neutral atmosphere profiles and a boolean value if the `lat`, `lon` (deg)
and UTC time `dt` is daytime. The profiles will be evaluated at heights `z` in km.
"""
function neutralprofiles(lat, lon, z, dt::DateTime)
    if !INITIALIZED[]
        SpaceIndices.init()
        INITIALIZED[] = true
    end

    lon = mod(lon, 360)

    jd = date_to_jd(dt)
    sza = zenithangle(lat, lon, dt)
    
    g_lat = deg2rad(lat)
    g_lon = deg2rad(lon)
    
    h = z.*1000  # m
    p = AtmosphericModels.nrlmsise00.(jd, h, g_lat, g_lon)  # #/m³
    
    f107 = space_index(Val(:F10adj), jd)
    nearest_f107 = FIRI.values(:f10_7)[argmin(abs.(f107 .- FIRI.values(:f10_7)))]

    Ne0 = firi(sza, lat; f10_7=nearest_f107, month=month(dt))
    Ne = FIRI.extrapolate(Ne0, h)
    # With specified f10_7 and month, there is only 1 profile after interpolation

    # Comprehension is more type-stable than `getfield.(p, :T_alt)` (if `:h` isn't Float64)
    df = Table(
        h  = collect(z*1.0),  # convert to Float64 for maximum type stability of Table
        Ne = Ne,
        Tn = [p[i].temperature for i in eachindex(p)],
        O  = [p[i].O_number_density for i in eachindex(p)],
        N2 = [p[i].N2_number_density for i in eachindex(p)],
        O2 = [p[i].O2_number_density for i in eachindex(p)]
    )

    # Sometimes a single (or maybe more?) height has NaN in MSIS. Not sure why...
    mask = isnan.(df.Tn)
    if any(mask)
        itp1 = interpolate(df.h[.!mask], df.Tn[.!mask], FritschButlandMonotonicInterpolation())
        df.Tn .= itp1.(df.h)
    end

    mask .= isnan.(df.O)
    if any(mask)
        itp2 = interpolate(df.h[.!mask], df.O[.!mask], FritschButlandMonotonicInterpolation())
        df.O .= itp2.(df.h)
    end

    mask .= isnan.(df.O2)
    if any(mask)
        itp3 = interpolate(df.h[.!mask], df.O2[.!mask], FritschButlandMonotonicInterpolation())
        df.O2 .= itp3.(df.h)
    end

    mask .= isnan.(df.N2)
    if any(mask)
        itp4 = interpolate(df.h[.!mask], df.N2[.!mask], FritschButlandMonotonicInterpolation())
        df.N2 .= itp4.(df.h)
    end
    
    return df
end

"""
    chargeprofiles(flux, lat, lon, ee, z, dt::DateTime; t=1e7) → (background_profiles, perturbed_profiles)
    chargeprofiles(flux, ee, neutraltable, z, daytime::Bool; t=1e7)

Compute GPI background and EPP-perturbed profiles for precipitating electron `flux` in
el/cm²/s, `lat` and `lon` in degrees, heights `z` in kilometers, and time `dt` for
`EnergeticElectrons` specified by `ee`.
"""
function chargeprofiles(flux, ee, neutraltable, z, daytime::Bool; t=1e7)
    if iszero(flux)
        Nspec0 = chargeprofiles(neutraltable, z, daytime; t)
        return Nspec0, copy(Nspec0)
    end

    # `md` must be defined at kilometer intervals, but the ionization profile can be at
    # finer `z` steps
    zstepped = first(z):last(z)
    p = GPI.Profiles(neutraltable)
    md = massdensity.((p,), zstepped)  # g/cm³
    S = ionizationprofile(ee, md/1000, z)  # pairs/el/cm
    S *= flux  # pairs/cm³/s

    Nspec0, Nspec = gpi(neutraltable, z, daytime, S*1e6; t)  # convert ionization to pairs/m³/s

    return Nspec0, Nspec
end

function chargeprofiles(flux, lat, lon, ee, z, dt::DateTime; t=1e7)
    neutraltable = neutralprofiles(lat, lon, z, dt)
    chargeprofiles(flux, ee, neutraltable, z, isday(zenithangle(lat, lon, dt)); t)
end

"""
    chargeprofiles(lat, lon, z, dt::DateTime; t=1e7) → background_profile
    chargeprofiles(neutraltable, z, daytime::Bool; t=1e7)

Return the unperturbed (zero flux) GPI background profiles only.
"""
function chargeprofiles(neutraltable, z, daytime::Bool; t=1e7)
    Nspec0, _ = equilibrium(neutraltable, z, daytime; t)

    return Nspec0
end

function chargeprofiles(lat, lon, z, dt::DateTime; t=1e7)
    neutraltable = neutralprofiles(lat, lon, z, dt)
    chargeprofiles(neutraltable, z, isday(zenithangle(lat, lon, dt)); t)
end
    
end # module
