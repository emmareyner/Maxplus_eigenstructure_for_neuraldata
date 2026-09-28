using Pkg

root = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(root)

@assert VERSION >= v"1.12.0" && VERSION < v"1.13.0" "This project targets Julia 1.12.x; got $(VERSION)."

Pkg.instantiate()
Pkg.precompile()

mp = only(info for info in values(Pkg.dependencies()) if info.name == "MaxPlus")
@assert v"0.5.0" <= mp.version < v"0.6.0" "Expected MaxPlus 0.5.x, got $(mp.version)."

println("Environment ready.")
println("Julia:   ", VERSION)
println("MaxPlus: ", mp.version)
println("Project: ", root)
