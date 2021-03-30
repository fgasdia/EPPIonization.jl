using Test, EPPIonization
using DelimitedFiles, Dates

using EPPIonization: datadir


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
end

@testset "EPPIonization" begin
    comparematlab()
    profiles()
end
