using Test
using NL2sol
using Glob

print_steps = false
verbose = false
include("util.jl")


for fn in glob("./optests/*.jl")
    include(fn)
end

@testset "NL2sol tests" begin
    @test rosenbrock(scale=3)
    @test helix(scale=3)
    @test singular(scale=3)
    @test woods(scale=3)
    @test zangwill(scale=3)

    @test engvall(scale=3)
    @test branin(scale=3)
    @test beale(scale=2)
    @test cragg_levy(scale=2)
    @test box(scale=2)
    #@test davidon(scale 1, verbose=true)  # false convergence
    @test freudenstein_roth(scale=3)
    @test watson6(scale=1)
    @test watson9(scale=1)
    @test watson12(scale=1)
    #@test watson20(scale=1, verbose=true) # false convergence
    @test chebyquad(scale=2)
    @test brown_dennis(scale=3)
    @test bard(scale=3)
    @test jennrich_sampson(scale=3)
    @test kowalik_osborne(scale=3)
    @test osborne1(scale=3)
    @test osborne2(scale=2)
    @test madsen(scale=3)
    #@test meyer(verbose=true)             # false convergence
    @test brown5()
    @test brown10()
    @test brown30()
    @test brown40()
end
