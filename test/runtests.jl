using Test
using SparseArrays
using MaxPlus: plustimes
using Maxplus_eigenvectors

const ROOT = normpath(joinpath(@__DIR__, ".."))

@testset "Maxplus_eigenvectors" begin
    @testset "Basic matrices" begin
        A = I(3)
        @test size(A) == (3, 3)
        @test all(A[i, i] == one(MP) for i in 1:3)
        @test all(A[i, j] == zero(MP) for i in 1:3 for j in 1:3 if i != j)
        @test size(nB(2, 3)) == (2, 3)
        @test size(nT(2, 3)) == (2, 3)
    end

    @testset "MaxPlus 0.5 compatibility" begin
        A = MP([0.0 1.0; 2.0 0.0])
        @test eltype(A) <: MP
        @test plustimes.(mmp_mpower(A, 2)) ≈ [3.0 1.0; 2.0 3.0]
        @test plustimes(karp_max_cycle_mean_mp(A)) ≈ 1.5

        U = ufnf₀(A; verbose=false)
        @test plustimes(U.λ) ≈ 1.5
    end

    @testset "Adjacency CSV loader" begin
        path = joinpath(ROOT, "data", "examples", "weighted_network.csv")
        loaded = load_adjacency_csv(path)

        @test loaded.labels == ["A", "B", "C", "D", "E", "F"]
        @test size(loaded.matrix) == (6, 6)
        @test loaded.matrix[2, 1] == 2.0  # A -> B
        @test loaded.matrix[1, 2] == 1.0  # B -> A
        @test loaded.matrix[5, 4] == -0.5 # D -> E
    end

    @testset "Edge-list CSV loader" begin
        path = joinpath(ROOT, "data", "examples", "weighted_edges.csv")
        loaded = load_edgelist_csv(path)

        @test loaded.labels == ["A", "B", "C", "D", "E", "F"]
        @test size(loaded.matrix) == (6, 6)
        @test loaded.matrix[2, 1] == 2.0
        @test loaded.matrix[5, 4] == -0.5
    end

    @testset "Numeric to max-plus conversion" begin
        C = sparse([1, 2], [2, 1], [1.0, 3.0], 2, 2)
        A = to_maxplus(C)

        @test A[1, 2] == MP(1.0)
        @test A[2, 1] == MP(3.0)
        @test A[1, 1] == zero(MP)

        dense = [0.0 1.0; 2.0 0.0]
        B = to_maxplus(dense)
        @test B[1, 1] == zero(MP)
        @test B[1, 2] == MP(1.0)
    end

    @testset "Generic network analysis" begin
        C = sparse([1, 2], [2, 1], [1.0, 2.0], 2, 2)
        result = analyse_network(
            C;
            labels=["x", "y"],
            verbose=false,
        )

        @test result.labels == ["x", "y"]
        @test result.left !== nothing
        @test check_eigenpairs(result.matrix, result.right; side=:right)
        @test size(finite_right_spectrum_summary(result.right), 1) > 0
        @test size(finite_right_eigenvector_support(result.right, result.labels), 1) > 0
    end

    @testset "Network summary" begin
        A = sparse([1, 2], [2, 1], [1.0, 2.0], 2, 2)
        summary = network_summary((example=A,))

        @test summary.Nodes[1] == 2
        @test summary.Edges[1] == 2
        @test summary.Symmetric[1] == false
    end
    @testset "C. elegans helper smoke test" begin
        include(joinpath(ROOT, "src", "celegans_helpers.jl"))
        @test isdefined(Main, :CElegansHelpers)
        @test CElegansHelpers.DEFAULT_COOK_FILE ==
            "SI 5 Connectome adjacency matrices, corrected July 2020 (1).xlsx"
    end

end
