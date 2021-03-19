using Test, EPPIonization
using DelimitedFiles

@testset "EPPIonization" begin
    ## Compare to Matlab
    # electron energy
    energy = exp10.(4:0.01:7)

    # energy distribution, f(E) ∝ exp(-E/100 keV)
    energydis = exp.(-energy/1e5)

    # pitch angle
    pitchangle = 0:90

    # pitch angle distribution, isotropic
    pitchdis = ones(length(pitchangle))

    # example atmosphere, first column: alt 0:500 km, second column: mass density (g/cm²)
    atom = readdlm("example_atmosphere.txt")

    ion = ionprofile(atom[:,1], energy, energydis, pitchangle, pitchdis, atom[:,2])

    mion = readdlm("example_ionpro.txt", ',')
    @test mion[:,2] ≈ ion

    @test_throws ArgumentError ionprofile(atom[:,1], energy, energydis, pitchangle, pitchdis, atom[1:10,2])
end
