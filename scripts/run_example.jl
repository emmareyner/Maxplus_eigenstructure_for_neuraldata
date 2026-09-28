using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using Maxplus_eigenvectors

root = normpath(joinpath(@__DIR__, ".."))
input = joinpath(root, "data", "examples", "weighted_network.csv")

example = load_adjacency_csv(input)
analysis = analyse_network(
    example.matrix;
    labels=example.labels,
    name="synthetic example",
)

println("\nRight projective spectrum")
display(finite_right_spectrum_summary(analysis.right))

println("\nRight eigendirection supports")
display(finite_right_eigenvector_support(analysis.right, analysis.labels))
