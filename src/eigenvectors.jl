# Max-plus semi-eigenvalue and eigenvector routines.

"""
    SemiEigen0{T <: Number}

Eigendecomposition result for an irreducible UFNF0 matrix.

- `values`: eigenvalue repeated for each critical cycle.
- `vectors`: right eigenvectors associated with each critical cycle.
- `cycles`: critical cycles in local node indices.
"""
@kwdef struct SemiEigen0{T <: Number}
    values::Array{T}  # Eigenvalue repeated for each critical cycle.
    vectors::Array{Matrix{T}}  # Right eigenvectors, coindexed with critical cycles.
    cycles::Array{Cycle}  # Critical cycles in local node indices.
end

"""
    eigen₀(A::Matrix{T}, ufnf0::UFNF0{T}; tol=1e-10) where T <: Number → SemiEigen0{T}

Compute the right eigendecomposition of an irreducible matrix from its UFNF0 representation.
"""
function eigen₀(
    A::Matrix{T},
    ufnf0::UFNF0{T};
    tol=1e-10,
) where T <: Number

    (; λ, nodes, cycles) = ufnf0

    # An acyclic SCC has no critical cycle and therefore generates
    # no local finite eigendirection.
    if isempty(cycles)
        return SemiEigen0{T}(
            T[],
            Matrix{T}[],
            Cycle[],
        )
    end

    Anormplus = tclosure(
        A / λ,
        UFNF0{T}(one(T), nodes, cycles),
    )

    return SemiEigen0{T}(
        map(c -> λ, cycles),
        map(c -> Anormplus[:, c], cycles),
        cycles,
    )
end

"""
    SemiEigen1{T <: Number}

Eigendecomposition result for a UFNF1 matrix.

Stores eigenvalues, right eigenvectors, associated critical nodes, and strongly connected components.
"""
@kwdef struct SemiEigen1{T <: Number}
    values::Array{T}  # Eigenvalues with multiplicities induced by critical cycles.
    vectors::Matrix{T}  # Right eigenvectors stored as columns.
    nodes::Array{Int}  # Critical-node indices coindexed with eigenvalues and eigenvectors.
    sccs::Vector{Vector{Int}}  # Strongly connected components.
end

"""
    eigen₁(Ar::Matrix{T}, ufnf1::UFNF1{T}; tol=1e-10) where T <: Number → SemiEigen1{T}

Compute the right eigendecomposition of a UFNF1 matrix.
"""
function eigen₁(Ar::Matrix{T}, ufnf1::UFNF1{T}; tol=1e-10) where T <: Number


    (; sccs, cdg, ufnf0s) = ufnf1

    m = size(sccs,1)


    lambdas = T[]
    nodes = Int[]
    vectors = Vector{Vector{T}}()

    for (i,scc) in enumerate(sccs)
        sp0 = eigen₀(Ar[scc,scc],ufnf0s[i];tol)

        if isempty(sp0.cycles)
            continue
        end

        # Retain nodes that belong to critical cycles.
        critical_nodes = reduce(vcat, [scc[cy] for cy in sp0.cycles])

        # Local eigenvalue of the current strongly connected component.
        ρ = sp0.values[1]

        # Normalize by the local eigenvalue.
        Arnorm = Ar/ρ
        Arnormp = tclosure(Arnorm,ufnf₁(Arnorm);tol)

        # Candidate columns correspond to critical nodes.
        candidate_vectors = Arnormp[:,critical_nodes]

        # Retain only columns satisfying the right-eigenvector equation.
        #     Ar * v = ρ * v

        for j in 1:size(candidate_vectors,2)
            vj = candidate_vectors[:,j]

            lhs = Ar * vj
            rhs = vj * ρ

            ok = true
            for idx in eachindex(lhs)
                if !((lhs[idx] == rhs[idx]) || (lhs[idx] ≈ rhs[idx]))
                    ok = false
                    break
                end
            end

            if ok
                push!(lambdas, ρ)
                push!(nodes, critical_nodes[j])
                push!(vectors, vj)
            end
        end

    end

    @assert length(lambdas) == length(vectors)
    @assert length(nodes) == length(vectors)

    vector_matrix = isempty(vectors) ? zeros(T,size(Ar,1),0) : reduce(hcat,vectors)


    return SemiEigen1{T}(
        lambdas,
        vector_matrix,
        nodes,
        sccs
    )
end

"""
    eigvals(s1::SemiEigen1)

Return the eigenvalues of `s1` in descending order.
"""
function eigvals(s1::SemiEigen1{T}) where T <: Number
    return(sort(s1.values, rev=true))
end


"""
    SemiEigen2{T <: Number}

Eigendecomposition result shared by UFNF2 and UFNF3 matrices.

Results are stored by weakly connected component, together with the corresponding node sets.
"""
@kwdef struct SemiEigen2{T <: Number}
    values::Array{Array{T}}  # Eigenvalues by weakly connected component.
    vectors::Array{Matrix{T}}  # Eigenvectors for each non-empty component.
    nodes::Array{Array{Int}}  # Node indices for each weakly connected component.
    wccs::Array{Array{Array{Int}}}  # Weakly connected components.
end

"""
    eigen₂(Ad::Matrix{T}, ufnf2::UFNF2{T}; tol=1e-10) where T <: Number → SemiEigen2{T}

Compute the right eigendecomposition of a matrix from its UFNF2 representation.
"""
function eigen₂(Ad::Matrix{T}, ufnf2::UFNF2{T}; tol=1e-10) where T <: Number
    # Unpack the UFNF2 decomposition.
    (; wccs, cdg, ufnf1s) = ufnf2
    n = size(Ad,1)

    # Allocate result containers.
    nWccs = length(wccs)
    values = Array{Array{T}}(undef, nWccs)
    vectors = Array{Matrix{T}}(undef, nWccs)
    nodes = Array{Array{Int}}(undef, nWccs)  # Node sets, not cycles.
    # Process each weakly connected component.
    for (k,wcc) in enumerate(wccs)

        # Reconstruct the original node order used when ufnf1s[k] was built.
        # Do NOT use reduce(union,wcc), because that returns SCC-order and
        # can destroy the local indexing expected by ufnf1s[k].
        original_wcc = Vector{Int}(undef, sum(length, wcc))

        for (local_scc, global_scc) in zip(ufnf1s[k].sccs, wcc)
            original_wcc[local_scc] = global_scc
        end

        wcc = original_wcc

        # Compute the local UFNF1 spectrum.
        sp1 = eigen₁(Ad[wcc,wcc],ufnf1s[k]; tol=tol)
        values[k] = sp1.values
        vectors[k] = sp1.vectors

        @assert size(sp1.vectors,1) == length(wcc)

        # The rows of sp1.vectors are in the same order as Ad[wcc,wcc],
        # so nodes[k] must store that same row order.
        nodes[k] = wcc

    end
    # Assemble the eigendecomposition.
    return SemiEigen2{T}(;
        values,
        vectors,
        nodes,
        wccs)

end

"""
    eigen₃(A::Matrix{T}; tol=1e-10) where T <: Number → SemiEigen2{T}

Compute the right eigendecomposition of a generic matrix.

The UFNF2 eigendecomposition of the non-empty `beta` block is lifted back to the full UFNF3 space.
"""
function eigen₃(A::Matrix{T}; tol=1e-10) where T <: Number
    ufnf3 = ufnf₃(A)
    (; iota, alpha, beta, omega, ufnf2) = ufnf3

    n = size(A,1)
    allNodes = collect(1:n)

    values = Array{T}[]
    vectors = Matrix{T}[]
    nodes = Array{Int}[]
    wccs = Array{Array{Int}}[]

    # 1) Right eigenvectors associated with empty columns
    # If column j is empty, then A * e_j = zero(T) * e_j.
    # Empty columns are iota ∪ alpha.
    allZeroCols = sort(vcat(iota, alpha))

    if !isempty(allZeroCols)
        Z = [(i == j ? one(T) : zero(T)) for i in 1:n, j in allZeroCols]

        push!(values, fill(zero(T), length(allZeroCols)))
        push!(vectors, Z)
        push!(nodes, allNodes)
        push!(wccs, [[j] for j in allZeroCols])
    end

    # 2) Proper eigenvectors from the beta UFNF2 core
    if !isempty(beta)

        # ufnf₃ currently returns beta as betareordered, but ufnf2 was built
        # on the original beta order. Recover the original beta order.
        flat_local_beta = reduce(vcat, reduce(vcat, ufnf2.wccs))

        beta_original = Vector{Int}(undef, length(beta))
        beta_original[flat_local_beta] = beta

        sp2 = eigen₂(A[beta_original,beta_original], ufnf2; tol=tol)

        for k in eachindex(sp2.values)
            vals = sp2.values[k]
            Vlocal = sp2.vectors[k]
            localNodes = sp2.nodes[k]

            Vfull = zeros(T, n, size(Vlocal,2))

            for j in 1:size(Vlocal,2)
                ρ = vals[j]

                # Vector on the beta block.
                vβ = fill(zero(T), length(beta_original))
                vβ[localNodes] = Vlocal[:,j]

                # Lift to full space.
                v = fill(zero(T), n)
                v[beta_original] = vβ

                # alpha has empty columns but non-empty rows.
                # Therefore alpha coordinates are forced by:
                #
                #     A[alpha,beta] * vβ = ρ * v_alpha
                #
                # hence:
                #
                #     v_alpha = (A[alpha,beta] * vβ) / ρ
                #
                if !isempty(alpha)
                    if ρ == zero(T)
                        tmp = A[alpha,beta_original] * vβ
                        @assert all(tmp .== zero(T))
                        v[alpha] = fill(zero(T), length(alpha))
                    else
                        v[alpha] = (A[alpha,beta_original] * vβ) / ρ
                    end
                end

                # iota and omega remain zero(T) for proper right eigenvectors.
                Vfull[:,j] = v
            end

            # Verify lifted right eigenvector equations.
            lhs = A * Vfull
            rhs = Vfull * diagm(vals)

            bad = CartesianIndex{2}[]

            for idx in CartesianIndices(lhs)
                if !((lhs[idx] == rhs[idx]) || (lhs[idx] ≈ rhs[idx]))
                    push!(bad, idx)
                end
            end

            if !isempty(bad)
                println("Failed lifted UFNF3 block k = ", k)
                println("beta returned by ufnf₃ = ", beta)
                println("beta_original recovered = ", beta_original)
                println("localNodes = ", localNodes)
                println("global beta nodes = ", beta_original[localNodes])
                println("bad positions = ", bad)
                println("lhs[bad] = ", lhs[bad])
                println("rhs[bad] = ", rhs[bad])
                error("UFNF3 right eigenvector equation failed")
            end

            push!(values, vals)
            push!(vectors, Vfull)
            push!(nodes, allNodes)
            push!(wccs, map(scc -> beta_original[scc], sp2.wccs[k]))
        end
    end

    return SemiEigen2{T}(;
        values,
        vectors,
        nodes,
        wccs
    )
end

function mpeigen(A::Matrix{T}; tol=1e-10) where T <: Number
    return eigen₃(A; tol=tol)
end

"""
    leigen₃(A::Matrix{T}; tol=1e-10) where T <: Number → SemiEigen2{T}

Compute left eigenvectors of `A` as right eigenvectors of `transpose(A)`.

Returned eigenvectors are stored as columns.
"""
function leigen₃(A::Matrix{T}; tol=1e-10) where T <: Number
    return eigen₃(Matrix(transpose(A)); tol=tol)
end

function lmpeigen(A::Matrix{T}; tol=1e-10) where T <: Number
    return leigen₃(A; tol=tol)
end

