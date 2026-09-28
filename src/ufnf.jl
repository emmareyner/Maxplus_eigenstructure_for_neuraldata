# Upper Frobenius Normal Form (UFNF) decompositions.

"""
    UFNF

Abstract type for Upper Frobenius Normal Form decompositions.
"""
abstract type UFNF{T <: Number}
end


"""
Type alias for a cycle represented by matrix-node indices.
"""
const Cycle = Vector{Int}
# Helper functions.


"""
    cycle2pairs(cc) → iterable of index pairs

Convert a cycle into the ordered matrix-index pairs corresponding to its edges.
"""
cycle2pairs(cc) = length(cc) == 1 ? (cc[1], cc[1]) : zip(cc, [cc[2:end]; cc[1]])

"""
    weight(A::Matrix{T}, cc::Cycle) → T

Return the semiring product of the edge weights along cycle `cc`.
"""
function weight(A::Matrix{T}, cc::Cycle)::T where {T}
    return reduce(*, A[CartesianIndex.(cycle2pairs(cc))]; init=one(T))
end


"""
    UFNF0{T <: Number} <: UFNF{T}

UFNF representation of an irreducible matrix, storing its maximum cycle mean and critical cycles.
"""
@kwdef struct UFNF0{T <: Number} <: UFNF{T}
    λ::T  # Critical cycle mean.

    nodes::Array{Int}  # Matrix node indices.
    cycles::Array{Cycle}  # Critical cycles.
end


"""
    UFNF1{T <: Number} <: UFNF{T}

UFNF representation of a reducible weakly connected matrix.

Stores the strongly connected components, their condensation digraph, and the UFNF0 decomposition of each irreducible component.
"""
@kwdef struct UFNF1{T <: Number} <: UFNF{T}
    sccs::Array{Array{Int}}  # Strongly connected components.
    cdg::SimpleDiGraph{Int}
    ufnf0s::Array{UFNF0{T}}  # UFNF0 decomposition of each irreducible component.
    UFNF1{T}(sccs,cdg,u) where {T<:Number} = new(sccs,cdg,u)
end


"""
    ufnf₁_raw(Aᵣ::Matrix{T}, n::Int) where T <: Number → UFNF1{T}

Construct the UFNF1 representation of a non-empty square matrix without empty rows or columns.
"""
function ufnf₁_raw(Aᵣ::Matrix{T}, n::Int) where  T <: Number
    dg = Graphs.DiGraph(Aᵣ .!= zero(T))
    sccs = strongly_connected_components(dg)  # Strongly connected components.
    cdg = Graphs.condensation(dg, sccs)  # Condensation digraph.

    # Build the UFNF1 decomposition.
    nsccs = length(sccs)  # Number of strongly connected components.
    if nsccs != 1
        order = topological_sort_by_dfs(cdg)
        sccs = sccs[order]  # Topological order yields the UFNF1 node order.

    end

    ufnf0s = Array{UFNF0{T}}(undef, nsccs)
    for s in 1:nsccs
        nodes = sccs[s]  # Nodes are local to this component.
        ufnf0s[s] = ufnf₀(Aᵣ[nodes,nodes])
    end
    return UFNF1{T}(sccs, cdg, ufnf0s)
end

"""
    ufnf₁(Aᵣ::Matrix{T}) where T <: Number → UFNF1{T}

Validate the matrix dimensions and construct its UFNF1 representation.
"""
function ufnf₁(Aᵣ::Matrix{T}) where T <: Number
    (n,m) = size(Aᵣ)
    @assert n == m && n != 0 "Can only do UFNF1 for a non-empty square matrix."
    return ufnf₁_raw(Aᵣ,n)
end


"""
    UFNF2{T <: Number} <: UFNF{T}

UFNF representation of a block-diagonal reducible matrix.

Stores the weakly connected components, their condensation digraphs, and the corresponding UFNF1 decompositions.
"""
@kwdef struct UFNF2{T <: Number}  <: UFNF{T}
    wccs::Array{Array{Array{Int}}}  # Weakly connected components.
    cdg::Array{SimpleDiGraph{Int}}  # Condensation digraphs.
    ufnf1s::Array{UFNF1{T}}  # UFNF1 decompositions of the weakly connected components.
    UFNF2{T}(wccs,cdg,u) where {T<:Number} = new(wccs,cdg,u)
end


"""
    ufnf₂_raw(Ad::Matrix{T}) where T <: Number → UFNF2{T}

Construct the UFNF2 representation from the weakly connected components of `Ad`.
"""
function ufnf₂_raw(Ad::Matrix{T}) where T <: Number
    # Find weakly connected components.
    dg = Graphs.DiGraph(Ad .!= zero(T))
    blocks = weakly_connected_components(dg)

    # Order components by size for a stable UFNF representation.
    sort!(blocks; by=length)

    # Build the block-diagonal UFNF2 representation.
    nwccs = length(blocks)  # Number of weakly connected components.

    wccs = Array{Array{Array{Int}}}(undef, nwccs)  # UFNF-ordered weakly connected components.
    cdg = Array{SimpleDiGraph{Int}}(undef,nwccs)  # Condensation digraphs.
    ufnf1s = Array{UFNF1{T}}(undef,nwccs)
    for w in 1:nwccs
        block = blocks[w]
        # Local node indices are relative to the current weakly connected component.
        ufnf1s[w] = ufnf₁_raw(Ad[block,block], length(block))
        wccs[w] = map(dSccs -> block[dSccs], ufnf1s[w].sccs)
        cdg[w] = ufnf1s[w].cdg

    end
    return UFNF2{T}(wccs,cdg,ufnf1s)
end

"""
    ufnf₂(Ad::Matrix{T}) where T <: Number → UFNF2{T}

Validate the matrix dimensions and construct its UFNF2 representation.
"""
function ufnf₂(Ad::Matrix{T}) where T <: Number
    (n,m) = size(Ad)
    n == m && n != 0 || error("Can only do UFNF2 for a non-empty square matrix.")
    return ufnf₂_raw(Ad)
end


"""
    UFNF3{T <: Number} <: UFNF{T}

Generic UFNF representation based on the four index classes `iota`, `alpha`, `beta`, and `omega`.

The non-empty `beta` core is represented by a UFNF2 decomposition when present.
"""
@kwdef struct UFNF3{T <: Number}  <: UFNF{T}
    iota::Vector{Int}  # Empty rows and columns.
    alpha::Vector{Int}  # Empty columns with non-empty rows.
    beta::Vector{Int}  # Non-empty rows and columns.
    omega::Vector{Int}  # Empty rows with non-empty columns.
    # The beta block contains the non-empty UFNF2 core.

    ufnf2::Union{UFNF2{T}, Nothing}  # No UFNF2 core exists for an empty beta block.
    UFNF3{T}(i,a,b,o,u) where {T <: Number} = new(i,a,b,o,u)
end
# Constructor for matrices with an empty beta block.
function UFNF3{T}(i,a,b,o) where {T <: Number}
    return UFNF3{T}(i,a,b,o,nothing)
end


"""
    ufnf₃(A::Matrix{T}) where T <: Number → UFNF3{T}

Construct the generic UFNF3 representation of a square matrix.

Rows and columns are partitioned into:
- `iota`: empty rows and empty columns,
- `alpha`: non-empty rows and empty columns,
- `beta`: non-empty rows and non-empty columns,
- `omega`: empty rows and non-empty columns.

The `beta` block is decomposed recursively as UFNF2.
"""
function ufnf₃(A::Matrix{T}) where T <: Number
    (n,m) = size(A)
    n == m || error("UFNF3 only applicable to square matrices")
    allZeroCols = collect(1:n)[all(A .== zero(T), dims=1)']
    # Empty matrix: only the empty-row/column class remains.
    length(allZeroCols) == n && return UFNF3{T}(allZeroCols,[],[],[])
    # Identify the remaining UFNF3 index classes.
    allZeroRows = collect(1:n)[all(A .== zero(T), dims=2)]
    # Decompose indices according to the UFNF3 classes.
    iota = allZeroCols ∩ allZeroRows  # Empty rows and columns.
    beta = setdiff(1:n, allZeroCols ∪ allZeroRows)  # Non-empty rows and columns.
    # Decompose the non-empty beta core as UFNF2.
    ufnf2 = ufnf₂_raw(A[beta,beta])


    betareordered =
        reduce(vcat,
            mapreduce(
                distalCcs -> map(distalSccs -> beta[distalSccs], distalCcs),
                vcat,
                ufnf2.wccs
            )
        )
    @assert Set(beta) == Set(betareordered)
    return UFNF3{T}(; iota,
                        alpha = setdiff(allZeroCols, iota),  # Empty columns with non-empty rows.
                        beta=betareordered,
                        omega = setdiff(allZeroRows, iota),  # Empty rows with non-empty columns.
                        ufnf2)
end


"""
    ufnf(A)

Return the generic UFNF representation of `A`.
"""
ufnf(A) = ufnf₃(A)


# Karp-based UFNF0 specialization for max-plus matrices.

# UFNF0 implementation using Karp's maximum cycle mean algorithm.
# The eigenvalue is computed without enumerating simple cycles. Critical nodes
# are recovered from the normalized closure for use by the downstream routines.

function karp_max_cycle_mean_mp(A::Matrix{T}) where {T<:MP}
    n, m = size(A)
    @assert n == m && n != 0 "Karp only applies to non-empty square matrices."

    z = zero(T)

    # Convert max-plus weights to conventional arithmetic.
    W = plustimes.(A)

    # Finite directed edges.
    edges = Tuple{Int, Int, Float64}[]

    for i in 1:n, j in 1:n
        if A[i, j] != z
            push!(edges, (i, j, W[i, j]))
        end
    end

    @assert !isempty(edges) "No finite edges found in this UFNF0 block."

    # D[k+1, j] stores the maximum weight of a length-k path ending at j.
    D = fill(-Inf, n + 1, n)
    D[1, :] .= 0.0

    for k in 1:n
        for (i, j, wij) in edges
            candidate = D[k, i] + wij
            if candidate > D[k + 1, j]
                D[k + 1, j] = candidate
            end
        end
    end

    # Karp's formula, with D[n+1, j] representing paths of length n.
    node_scores = fill(-Inf, n)

    for j in 1:n
        if isfinite(D[n + 1, j])
            candidates = Float64[]

            for k in 0:n-1
                if isfinite(D[k + 1, j])
                    push!(candidates, (D[n + 1, j] - D[k + 1, j]) / (n - k))
                end
            end

            if !isempty(candidates)
                node_scores[j] = minimum(candidates)
            end
        end
    end

    λ_float = maximum(node_scores)

    @assert isfinite(λ_float) "Karp failed to find a finite maximum cycle mean."

    return T(λ_float)
end


"""
    ufnf₀(A::Matrix{T}; verbose::Bool = true) where T<:MP → UFNF0{T}

Compute the UFNF0 representation of a non-empty square max-plus matrix using
Karp's maximum cycle mean algorithm. A singleton SCC without a self-loop is
kept in the UFNF structure with bottom eigenvalue and no critical cycle.
"""
function ufnf₀(A::Matrix{T}; verbose::Bool = true) where {T<:MP}
    n, m = size(A)
    @assert n == m && n != 0 "Can only do UFNF0 for non-empty square matrices."

    # A singleton SCC without a self-loop belongs to the condensation
    # structure but contains no finite directed cycle.
    if n == 1 && A[1, 1] == zero(T)
        verbose && println("ufnf₀: acyclic singleton SCC; no finite cycle.")
        return UFNF0{T}(zero(T), [1], Cycle[])
    end

    # Maximum cycle mean.
    λ = karp_max_cycle_mean_mp(A)

    # Normalize so that critical cycles have mean one(MP).
    Anorm = A / λ

    # Critical nodes have unit diagonal entries in the normalized closure.
    Anormplus = Anorm * (I(n) + Anorm)^(n - 1)

    critical_nodes = Int[]

    for i in 1:n
        if (Anormplus[i, i] == one(T)) || (Anormplus[i, i] ≈ one(T))
            push!(critical_nodes, i)
        end
    end

    # Numerical fallback for small floating-point deviations.
    if isempty(critical_nodes)
        verbose && println("ufnf₀/Karp warning: no critical nodes detected exactly.")
        verbose && println("Using maximal diagonal entries of normalized closure as fallback.")

        diag_values = [plustimes(Anormplus[i, i]) for i in 1:n]
        max_diag = maximum(diag_values)

        for i in 1:n
            if abs(diag_values[i] - max_diag) < 1e-8
                push!(critical_nodes, i)
            end
        end
    end

    @assert !isempty(critical_nodes) "No critical nodes detected after Karp normalization."

    # Downstream routines use these entries to select critical closure columns.
    critical_cycles = Cycle[[i] for i in critical_nodes]

    verbose && println("ufnf₀: using Karp maximum cycle mean.")
    verbose && println("  λ = ", λ)
    verbose && println("  plustimes(λ) = ", plustimes(λ))
    verbose && println("  number of representative critical nodes = ", length(critical_nodes))

    return UFNF0{T}(λ, collect(1:n), critical_cycles)
end