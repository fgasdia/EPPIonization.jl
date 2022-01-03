using Test, EPPIonization
using EPPIonization: datadir

using DelimitedFiles, Dates
using LMPTools

function comparematlab()
    # electron energy
    energy = exp10.(4:0.01:7)

    # energy distribution, f(E) ∝ exp(-E/100 keV)
    energydis = exp.(-energy/1e5)

    # pitch angle
    pitchangle = 0:90

    # pitch angle distribution, isotropic
    pitchdis = ones(length(pitchangle))

    # example atmosphere, first column: alt 0:500 km, second column: mass density (g/cm²)
    atmo = readdlm(datadir("example_atmosphere.txt"))

    ion = ionizationprofile(energy, energydis, pitchangle, pitchdis, atmo[:,2], atmo[:,1])

    mion = readdlm(datadir("example_ionpro.txt"), ',')
    @test mion[:,2] ≈ ion

    @test_throws ArgumentError ionizationprofile(energy, energydis, pitchangle, pitchdis, atmo[1:10,2], atmo[:,1])
end

function profiles()
    z = 0:110
    dt = DateTime(2020, 1, 1, 2, 30)
    lat, lon = 60, 258
    flux = 1e5

    daytime = isday(zenithangle(lat, lon, dt))  # == false

    np = neutralprofiles(lat, lon, z, dt)

    energy = 90e3:1e4:2.2e6  # 90 keV to 2.2 MeV every 10 keV
    energydis = exp.(-energy/1e5)  # e.g.: f(E) ∝ exp(-E/100 keV)
    pitchangle = 0:90
    pitchdis = ones(length(pitchangle))
    ee = EnergeticElectrons(energy, energydis, pitchangle, pitchdis)

    background1, perturbed1 = chargeprofiles(flux, lat, lon, ee, z, dt)
    background2, perturbed2 = chargeprofiles(flux, ee, np, z, daytime)

    @test background1 == background2
    @test perturbed1 == perturbed2

    # Zero-flux profiles
    background5a = chargeprofiles(lat, lon, z, dt)
    background5b = chargeprofiles(np, z, daytime)
    background6, _ = chargeprofiles(0, lat, lon, ee, z, dt)

    @test background5a == background5b
    @test background5a == background6
    
    # Non-integer steps
    zfine = 0:0.25:110
    np = neutralprofiles(lat, lon, zfine, dt)
    background3, perturbed3 = chargeprofiles(flux, lat, lon, ee, zfine, dt)
    background4, perturbed4 = chargeprofiles(flux, ee, np, zfine, daytime)
    
    @test background3 == background4
    @test perturbed3 == perturbed4

    @test background1 ≈ background3[1:4:end,:] rtol=0.01
    @test background2 ≈ background4[1:4:end,:] rtol=0.01
end

function extraarguments()
    z = 0:110
    dt = DateTime(2020, 1, 1, 2, 30)
    lat, lon = 60, 258
    flux = 1e5

    daytime = isday(zenithangle(lat, lon, dt))  # == false

    energy = 90e3:1e4:2.2e6  # 90 keV to 2.2 MeV every 10 keV
    energydis = exp.(-energy/1e5)  # e.g.: f(E) ∝ exp(-E/100 keV)
    pitchangle = 0:90
    pitchdis = ones(length(pitchangle))
    ee = EnergeticElectrons(energy, energydis, pitchangle, pitchdis)

    # Accuracy check
    np = neutralprofiles(lat, lon, z, dt)
    npf = neutralprofiles(lat, lon, z, dt; datafilepath="wdc")
    @test npf == np

    # Smoke tests
    @test chargeprofiles(flux, lat, lon, ee, z, dt; datafilepath="wdc") isa Tuple{<:Matrix,<:Matrix}
    @test chargeprofiles(lat, lon, z, dt; datafilepath="wdc") isa Matrix
end

@testset "EPPIonization" begin
    comparematlab()
    profiles()
    extraarguments()
end
