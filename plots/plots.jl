using DelimitedFiles, Dates
using EPPIonization
using Plots

using EPPIonization: datadir, massdensity, SPECIES

mask(x) = x < 1 ? NaN : x

function neutralprofileplot()
    z = 0:110
    dt = DateTime(2020, 1, 1, 2, 30)
    lat, lon = 60, 258
    
    np, _ = neutralprofiles(lat, lon, z, dt)

    plot(xlabel="Number density (m⁻³)", ylabel="Altitude (km)",
         xscale=:log10, ylims=(0, 115))
    plot!(mask.(np.Ne), np.h, label="Ne")
    plot!(mask.(np.O2), np.h, label="O2")
    plot!(mask.(np.N2), np.h, label="N2")
    plot!(mask.(np.O), np.h, label="O")
end

function massdensityplot()
    z = 0:110
    dt = DateTime(2020, 1, 1, 2, 30)
    lat, lon = 60, 258
    
    np, _ = neutralprofiles(lat, lon, z, dt)
    md = massdensity.((np,), z)  # g/cm³

    plot(md/1000, z, label="MSIS",
         xlabel="Neutral mass density (g/cm³/1000)", ylabel="Altitude (km)",
         xscale=:log10, ylims=(0, 115), xlims=(10^-14, 10^-5), xticks=exp10.(-12:3:-5))
    exatm = readdlm(datadir("example_atmosphere.txt"))
    plot!(exatm[:,2], exatm[:,1], label="EPP example")
end

function ionizationrateplot()
    z = 0:200
    dt = DateTime(2020, 1, 1, 2, 30)
    lat, lon = 60, 258

    energy = 90e3:1e4:2.2e6  # eV; 90 keV to 2.2 MeV
    energydis = exp.(-energy/2e5)  # f(E) ∝ exp(-E/β) where β ranges from 100 to 300 keV
    pitchangle = 0:90
    pitchdis = ones(length(pitchangle))
    flux = 1e5

    np, _ = neutralprofiles(lat, lon, z, dt)
    zstepped = first(z):last(z)
    md = massdensity.((np,), zstepped)  # g/cm³

    S = ionizationprofile(z, energy, energydis, pitchangle, pitchdis, md/1000)*1e6
    S *= flux
    
    plot(S, z,
         xlabel="Ionization rate (pairs/m³/s)", ylabel="Altitude (km)",
         xscale=:log10, ylims=(0, 205), xlims=(1e3, 1e9), legend=false)
end

function chargeprofileplot(daytime::Bool)
    z = 0:110
    lat, lon = 60, 258
    flux = 1e5

    if daytime
        dt = DateTime(2020, 1, 1, 18, 30)
    else
        dt = DateTime(2020, 1, 1, 2, 30)
    end

    np, daytime = neutralprofiles(lat, lon, z, dt)
    cp0, cp = chargeprofiles(flux, np, z, daytime)

    lws = fill(1.2, (1, length(SPECIES)))
    lws[1] = 2  # Ne

    pal = palette(:tab10)[1:5]
    plot(mask.(cp), z, xscale=:log10,
         labels=permutedims([SPECIES...]), color_palette=pal, linestyle=:solid,
         xlabel="Density (m⁻³)", ylabel="Altitude (km)", 
         yticks=0:20:120, xticks=exp10.([0, 3, 6, 9, 12]),
         legend=:outerright, xlims=(10^0, 10^12), linewidth=lws)
    plot!(mask.(cp0), z, linestyle=:dash, linewidth=lws,
          color_palette=pal, labels=false)
    plot!(mask.(np.Ne), z, color="black", linewidth=2, linestyle=:dash, label="FIRI")
end

function fineprofileplot(daytime::Bool)
    z = 0:110
    zfine = 0:0.25:110
    lat, lon = 60, 258
    flux = 1e5

    if daytime
        dt = DateTime(2020, 1, 1, 18, 30)
    else
        dt = DateTime(2020, 1, 1, 2, 30)
    end

    npf, daytime = neutralprofiles(lat, lon, zfine, dt)
    cp0, cp = chargeprofiles(flux, lat, lon, z, dt)
    cp0f, cpf = chargeprofiles(flux, lat, lon, zfine, dt)

    lws = fill(1.2, (1, length(SPECIES)))
    lws[1] = 2  # Ne

    pal = palette(:tab10)[1:5]
    plot(mask.(cp), z, xscale=:log10,
         labels=permutedims([SPECIES...]), color_palette=pal, linestyle=:solid,
         xlabel="Density (m⁻³)", ylabel="Altitude (km)", 
         yticks=0:20:120, xticks=exp10.([0, 3, 6, 9, 12]),
         legend=:outerright, xlims=(10^0, 10^12), linewidth=lws)
    plot!(mask.(cpf), zfine, linestyle=:dash, linewidth=lws,
          color_palette=pal, labels=false)
end

function relaxationplot(daytime::Bool)
    z = 0:110
    zfine = 0:0.25:110
    lat, lon = 60, 258
    flux = 1e5

    if daytime
        dt = DateTime(2020, 1, 1, 18, 30)
    else
        dt = DateTime(2020, 1, 1, 2, 30)
    end

    ts = exp10.(1:8)

    Nes = Matrix{Float64}(undef, length(z), length(ts))
    Ne0s = similar(Nes)
    for i in eachindex(ts)
        Nspec0, Nspec = chargeprofiles(flux, lat, lon, z, dt; t=ts[i])
        Nes[:,i] = Nspec[:,1]
        Ne0s[:,i] = Nspec0[:,1]
    end
    
    # Turns out t doesn't matter much in terms of runtime
    plot(mask.(Nes), z, xscale=:log10,
         xlabel="Density (m⁻³)", ylabel="Altitude (km)", legendtitle="10^",
         yticks=0:20:120, xticks=exp10.([0, 3, 6, 9, 12]), xlims=(10^0, 10^12),
         linewidth=1.2, labels=permutedims(log10.(ts)), legend=:topleft)
    plot!(mask.(Ne0s), z, linestyle=:dash, linewidth=1.2,
          labels=false)
end
