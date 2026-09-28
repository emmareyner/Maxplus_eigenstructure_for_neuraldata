# Transitive and reflexive-transitive closures.

"""
    tclosure(A::Matrix{T}, ufnf0::UFNF0{T}; tol=1e-10) where T <: Number → Matrix{T}

Return the transitive closure of an irreducible matrix with a precomputed UFNF0 representation.
"""
function tclosure(A::Matrix{T}, ufnf0::UFNF0{T}; tol=1.0e-10) where T <: Number
    n = size(A,1)
    # A superunit cycle mean leads to saturation.
    if ufnf0.λ > one(T) + tol
        return(nT(n))
    else  # Otherwise the finite closure is obtained directly.
        return A * (I(n) + A)^(n-1)
    end
end

"""
    tclosure(A::Matrix{T}, ufnf1::UFNF1{T}; tol=1e-10) where T <: Number → Matrix{T}

Return the transitive closure of a reducible matrix with a precomputed UFNF1 representation.
"""
function tclosure(A::Matrix{T}, ufnf1::UFNF1{T}; tol=1.0e-10) where T <: Number
    n = size(A,1)
    Aplus = nB(n)  # Allocate the result.

    # Access the component decomposition.
    (; sccs, cdg, ufnf0s) = ufnf1
    m = length(sccs)
    Astars = Array{Matrix{T}}(undef, m)  # Reflexive closures of component blocks.
    lambdas = map(u -> u.λ, ufnf0s)

    # Build the closure recursively over the condensation digraph.

    for (k,scc) in Iterators.reverse(enumerate(sccs))
        # Superunit cycle means saturate the corresponding block.
        if (lambdas[k] > one(T) + tol)

            Astars[k] = nT(length(scc))
            Aplus[scc,scc] = Astars[k]


        else        # Finite component closure.
            Aplus[scc,scc] = tclosure(A[scc,scc], ufnf0s[k]; tol=tol)
            Astars[k] = rclosure(Aplus[scc,scc])
        end
        successors = reduce(union,bfs_tree(cdg,k; dir=:in).fadjlist)
        if !isempty(successors)
            succnodes = reduce(vcat,sccs[successors])
            Aplus[scc,succnodes] =
                Astars[k] * A[scc,succnodes] * rclosure(Aplus[succnodes,succnodes])
        end
    end


    return Aplus
end

"""
    tclosure(A::Matrix{T}, ufnf2::UFNF2{T}; tol=1e-10) where T <: Number → Matrix{T}

Return the transitive closure of a UFNF2 matrix by processing its weakly connected UFNF1 blocks.
"""
function tclosure(A::Matrix{T}, ufnf2::UFNF2{T}; tol=1.0e-10) where T <: Number
    n = size(A,1)
    (; wccs, ufnf1s) = ufnf2
    A2plus = nB(n)
    for (w,wcc) in enumerate(wccs)
        theNodes = reduce(vcat,wcc)
        A2plus[theNodes,theNodes] = tclosure(A[theNodes,theNodes], ufnf1s[w]; tol=tol)
    end
    return A2plus
end

"""
    tclosure(A::Matrix{T}, ufnf3::UFNF3{T}; tol=1e-10) where T <: Number → Matrix{T}

Return the transitive closure of a matrix with a precomputed UFNF3 representation.
"""
function tclosure(A::Matrix{T}, ufnf3::UFNF3{T}; tol=1.0e-10) where T <: Number
    (; alpha, beta, omega) = ufnf3
    # Return immediately when no cycles are present.
    if isempty(beta) return nB(size(A,1)) end

    A1plus = nB(size(A,1))  # Allocate the result.


    # Compute the closure on the non-empty beta block.
    A1plus[beta,beta] = tclosure(A[beta,beta], ufnf₂(A[beta,beta]))


    if !isempty(omega)  # Fill beta-to-omega entries.
        A1plus[beta,omega] = rclosure(A1plus[beta,beta]) * A[beta,omega]
    end
    # Iota rows and columns remain empty.
    if !isempty(alpha)  # Fill alpha rows.
        A1plus[alpha,beta] = A[alpha,beta] * rclosure(A1plus[beta,beta])
        if !isempty(omega)  # Fill alpha-to-omega entries.
            A1plus[alpha,omega] = A1plus[alpha,beta] * A[beta,omega] + A[alpha,omega]
        end
    end
     return A1plus
end

"""
    tclosure(A::Matrix{<:MP}) → Matrix{<:MP}

Return the max-plus transitive closure, or Kleene plus, `A^+`.
"""
function tclosure(A::Matrix{T}) where {T<:MP}
    m,n = size(A)
    m == n || error("Can only find tclosure on a square matrix.")
    n == 1 && return (A[n,n] > 0 ? nT(n) : A)
    return tclosure(A,ufnf(A))
end

"""
    trclosure(A) → Matrix

Return the reflexive transitive closure `A^* = I + A^+`.
"""
trclosure(A) = rclosure(tclosure(A))

