# Basic max-plus utilities and display configuration.

"""
    configure_tropical_display!()

Configure compact tropical-number display and plain-text rendering inside notebooks.
"""
function configure_tropical_display!()
    if isdefined(MaxPlus, :set_tropical_display)
        MaxPlus.set_tropical_display(1)
    end
    return nothing
end
Base.show(io::IO, ::MIME"text/latex", x::MP) = show(io, MIME"text/plain", x)
Base.show(io::IO, ::MIME"text/latex", A::AbstractArray{<:MP}) = show(io, MIME"text/plain", A)

"""
Max-plus identity matrix.
"""
function I(n::Int, m::Int)
    return MP([(i != j ? -Inf : 0.0) for i in 1:n, j in 1:m])  # Max-plus identity matrix.
end
I(n) = I(n,n)

"""
Return an `n × m` max-plus bottom matrix. With one argument, return a square matrix.
"""
function nB(n::Int, m::Int)
    return fill(zero(MP), (n,m))  # Concrete Float64 max-plus bottom matrix (MaxPlus.jl >= 0.5).

end
nB(n::Int) = nB(n,n)

"""
Return an `n × m` max-plus top matrix. With one argument, return a square matrix.
"""
function nT(n::Int, m::Int)
    return fill(mptop,(n,m))  # Max-plus top element.
end
nT(n::Int) = nT(n,n)

const MP_ISAPPROX_DIST = 1e-10

Base.isapprox(a::MP, b::MP) =
    isapprox(plustimes(a) .+ MP_ISAPPROX_DIST, plustimes(b) .+ MP_ISAPPROX_DIST)

Base.isapprox(A::AbstractMatrix{<:MP}, B::AbstractMatrix{<:MP}) =
    plustimes(A) .+ MP_ISAPPROX_DIST ≈ plustimes(B) .+ MP_ISAPPROX_DIST

"""
    mmp_mpower_raw(A::Matrix{<:MP}, k::Int64; n = size(A,1)) → Matrix{<:MP}

Compute the `k`-th max-plus matrix power by repeated squaring.

Preconditions: `A` is square and `k ≥ 0`.
"""
function mmp_mpower_raw(A::Matrix{T}, k::Int64; n = size(A,1)) where {T<:MP}
    # Preconditions: A is square and k is non-negative.
    k == 1 && return A  # Early return.
    k == 0 && return [T(i != j ? -Inf : 0.0) for i in 1:n, j in 1:n]  # Max-plus identity matrix.
    Y = [T(i != j ? -Inf : 0.0) for i in 1:n, j in 1:n]

    # Accumulator starts as the identity matrix.
    while k > 0
        if isodd(k)
            Y = Y * A
        end
        A = A * A
        k >>= 1  # Divide k by two.
    end
    return Y  # Return the accumulated power.
end

"""
    mmp_mpower(A::Matrix{<:MP}, k::Int64) → Matrix{<:MP}

Compute the `k`-th max-plus power of a square matrix using repeated squaring.

- `A`: square matrix.
- `k`: non-negative exponent.

Returns `A^k`.
"""
function mmp_mpower(A::Matrix{T}, k::Int64) where {T<:MP}
    m,n = size(A)
    m == n || error("Can only find k-th power on a square matrix.")
    k >= 0 || error("Can only do positive powers by this method.")
    return mmp_mpower_raw(A, k; n = n)
end

"""
    rclosure(A::Matrix{T}) where T <: Number → Matrix{T}

Return the reflexive closure of `A`.
"""
function rclosure(A::Matrix{T}) where T <: Number
    (m,n) = size(A)
    for i in 1:min(n,m)
        A[i,i] = max(A[i,i], one(T))  # Ensure the diagonal contains the multiplicative identity.
    end
    return A
end

"""Return the top element of the requested tropical algebra."""
top(::Type{MP}) = MP(Inf)
top(::Type{MI}) = MI(-Inf)
