using Pkg

root = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(root)

@assert VERSION >= v"1.12.0" && VERSION < v"1.13.0" "Expected Julia 1.12.x, got $(VERSION)."

mp = only(info for info in values(Pkg.dependencies()) if info.name == "MaxPlus")
@assert v"0.5.0" <= mp.version < v"0.6.0" "Expected MaxPlus 0.5.x, got $(mp.version)."

println("OK: Julia ", VERSION, " / MaxPlus ", mp.version)
