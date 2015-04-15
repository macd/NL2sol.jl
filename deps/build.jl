# OK, there is no real build.jl. This is just a sham until the
# time I receive BinDeps enlightenment

# using BinDeps
# @BinDeps.setup

# libnl2sol = library_dependency("libnl2sol", aliases=["libnl2sol"])

# dpsrc = joinpath(Pkg.dir(), "Nl2sol/deps/src")
# dpbld = joinpath(dpsrc, "build")

# println("dpsrc ", dpsrc)
# println("dpbld ", dpbld)

# @build_steps begin
#     ChangeDirectory(dpsrc)
#     @build_steps begin
#         CreateDirectory(dpbld)
#         @build_steps begin
#             ChangeDirectory(dpbld)
#             @build_steps begin
#                 `cmake ..`
#                 `make`
#             end
#         end
#     end
# end

# @BinDeps.install Dict(:libnl2sol => :libnl2sol)
