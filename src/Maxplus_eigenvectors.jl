module Maxplus_eigenvectors

import CSV
import Graphs
import MaxPlus

using Graphs: SimpleDiGraph, strongly_connected_components, topological_sort_by_dfs, weakly_connected_components, bfs_tree
using MaxPlus: MP, MI, ε, mptop, plustimes
using LinearAlgebra: diagm, issymmetric
using SparseArrays: SparseMatrixCSC, sparse, findnz, nnz
using DataFrames: DataFrame, nrow, ncol, eachrow, names

include("primitives.jl")
include("ufnf.jl")
include("closures.jl")
include("eigenvectors.jl")
include("validation.jl")
include("critical_structure.jl")
include("scc_structure.jl")
include("network.jl")

export MP, MI, ε
export configure_tropical_display!, I, nB, nT, top
export mmp_mpower_raw, mmp_mpower, rclosure
export Cycle, cycle2pairs, weight
export UFNF, UFNF0, UFNF1, UFNF2, UFNF3
export ufnf₀, ufnf₁_raw, ufnf₁, ufnf₂_raw, ufnf₂, ufnf₃, ufnf
export karp_max_cycle_mean_mp
export tclosure, trclosure
export SemiEigen0, SemiEigen1, SemiEigen2
export eigen₀, eigen₁, eigen₂, eigen₃, eigvals, mpeigen, leigen₃, lmpeigen
export check_mp_matrix_equality, check_eigenpairs, check_right_semieigen, check_left_semieigen
export projectively_equivalent, projective_eigenvector_groups
export finite_right_spectrum_summary, finite_right_eigenvector_support
export load_adjacency_csv, load_edgelist_csv, to_maxplus, network_summary, analyse_network
export eigen_direction
export saturation_edges
export critical_cycle_mean
export critical_structure
export critical_structure_summary
export visualize_critical_structure
export analyse_scc_structure
export visualize_scc_structure

end # module Maxplus_eigenvectors
