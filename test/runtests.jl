using Test, EPPIonization
using EPPIonization: datadir

using DelimitedFiles, Dates


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
    atom = readdlm(datadir("example_atmosphere.txt"))

    ion = ionizationprofile(atom[:,1], energy, energydis, pitchangle, pitchdis, atom[:,2])

    mion = readdlm(datadir("example_ionpro.txt"), ',')
    @test mion[:,2] ≈ ion

    @test_throws ArgumentError ionizationprofile(atom[:,1], energy, energydis, pitchangle, pitchdis, atom[1:10,2])
end

function profiles()
    z = 0:110
    dt = DateTime(2020, 1, 1, 2, 30)
    lat, lon = 60, 258
    flux = 1e5

    np, daytime = neutralprofiles(lat, lon, z, dt)

    @test daytime == false

    background1, perturbed1 = chargeprofiles(flux, lat, lon, z, dt)
    background2, perturbed2 = chargeprofiles(flux, np, z, daytime)

    @test background1 == background2
    @test perturbed1 == perturbed2

    # Zero-flux profiles
    background5a, _ = chargeprofiles(lat, lon, z, dt)
    background5b, _ = chargeprofiles(np, z, daytime)
    background6, _ = chargeprofiles(0, lat, lon, z, dt)

    @test background5a == background5b
    @test background5a ≈ background6
    
    # Non-integer steps
    zfine = 0:0.25:110
    np, daytime = neutralprofiles(lat, lon, zfine, dt)
    background3, perturbed3 = chargeprofiles(flux, lat, lon, zfine, dt)
    background4, perturbed4 = chargeprofiles(flux, np, zfine, daytime)
    
    @test background3 == background4
    @test perturbed3 == perturbed4

    @test background1 ≈ background3[1:4:end,:] rtol=0.01
    @test background2 ≈ background4[1:4:end,:] rtol=0.01
end

@testset "EPPIonization" begin
    comparematlab()
    profiles()
end
