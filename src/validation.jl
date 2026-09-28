# Validation and reporting helpers for eigendecompositions.

# Verify the right and left eigenvector equations.

function check_mp_matrix_equality(lhs, rhs; label="")
    bad = CartesianIndex{2}[]

    for idx in CartesianIndices(lhs)
        if !((lhs[idx] == rhs[idx]) || (lhs[idx] ≈ rhs[idx]))
            push!(bad, idx)
        end
    end

    if !isempty(bad)
        println("Failed check: ", label)
        println("bad positions = ", bad)
        println("lhs[bad] = ", lhs[bad])
        println("rhs[bad] = ", rhs[bad])
        error("Matrix equality check failed: $label")
    else
        println("Passed check: ", label)
    end
end


function check_right_semieigen(A::AbstractMatrix{<:MP}, sp; tol=1e-10)
    for k in eachindex(sp.values)
        vals = sp.values[k]
        V = sp.vectors[k]

        @assert length(vals) == size(V, 2)

        lhs = A * V
        rhs = V * diagm(vals)

        bad = CartesianIndex{2}[]

        for idx in CartesianIndices(lhs)
            if !((lhs[idx] == rhs[idx]) || (lhs[idx] ≈ rhs[idx]))
                push!(bad, idx)
            end
        end

        if isempty(bad)
            println("Passed right eigenvector check for block $k")
        else
            println("Failed right eigenvector check for block $k")
            println("Number of bad positions: ", length(bad))
            println("First bad positions:")
            display(bad[1:min(end, 10)])
            println("lhs values:")
            display(lhs[bad[1:min(end, 10)]])
            println("rhs values:")
            display(rhs[bad[1:min(end, 10)]])
            error("Right eigenvector check failed for block $k")
        end
    end
end

"""
    projectively_equivalent(v, w; tol=1e-10)

Return `true` when two max-plus vectors represent the same projective
direction: they have the same support and differ only by an additive
constant on that support.
"""
function projectively_equivalent(v, w; tol=1e-10)
    support_v = findall(x -> x != zero(MP), v)
    support_w = findall(x -> x != zero(MP), w)

    support_v == support_w || return false
    isempty(support_v) && return true

    a = plustimes.(v[support_v])
    b = plustimes.(w[support_w])

    differences = a .- b
    offset = differences[1]

    return all(isapprox.(differences, offset; atol=tol, rtol=0.0))
end


"""
    projective_eigenvector_groups(values, V; tol=1e-10)

Group eigenvector columns representing the same max-plus projective
eigendirection.
"""
function projective_eigenvector_groups(values, V; tol=1e-10)
    groups = Vector{Vector{Int}}()
    used = falses(length(values))

    for j in eachindex(values)
        used[j] && continue

        λ = values[j]

        if λ == zero(MP)
            continue
        end

        group = [j]
        used[j] = true

        for q in (j + 1):length(values)
            used[q] && continue

            same_lambda =
                values[q] == λ ||
                isapprox(
                    plustimes(values[q]),
                    plustimes(λ);
                    atol=tol,
                    rtol=0.0,
                )

            if same_lambda &&
               projectively_equivalent(V[:, j], V[:, q]; tol=tol)

                push!(group, q)
                used[q] = true
            end
        end

        push!(groups, group)
    end

    return groups
end

function finite_right_spectrum_summary(sp; tol=1e-10)
    rows = DataFrame(
        Block = Int[],
        Direction = Int[],
        LambdaMP = MP[],
        Lambda = Float64[],
        Multiplicity = Int[],
        SupportSize = Int[]
    )

    for k in eachindex(sp.values)
        vals = sp.values[k]
        V = sp.vectors[k]

        groups = projective_eigenvector_groups(vals, V; tol=tol)

        for (direction, group) in enumerate(groups)
            j = first(group)

            λ = vals[j]
            v = V[:, j]
            support = findall(x -> x != zero(MP), v)

            push!(rows, (
                k,
                direction,
                λ,
                plustimes(λ),
                length(group),
                length(support)
            ))
        end
    end

    return rows
end

function finite_right_eigenvector_support(sp, labels; tol=1e-10)
    rows = DataFrame(
        Block = Int[],
        Direction = Int[],
        Lambda = Float64[],
        Multiplicity = Int[],
        SupportSize = Int[],
        SupportLabels = String[]
    )

    for k in eachindex(sp.values)
        vals = sp.values[k]
        V = sp.vectors[k]

        groups = projective_eigenvector_groups(vals, V; tol=tol)

        for (direction, group) in enumerate(groups)
            j = first(group)

            λ = vals[j]
            v = V[:, j]
            support = findall(x -> x != zero(MP), v)

            push!(rows, (
                k,
                direction,
                plustimes(λ),
                length(group),
                length(support),
                join(labels[support], ", ")
            ))
        end
    end

    return rows
end


function check_left_semieigen(A::AbstractMatrix{<:MP}, sp)
    for k in eachindex(sp.values)
        vals = sp.values[k]
        U = sp.vectors[k]

        @assert length(vals) == size(U, 2)

        lhs = Matrix(transpose(U)) * A
        rhs = diagm(vals) * Matrix(transpose(U))

        bad = CartesianIndex{2}[]

        for idx in CartesianIndices(lhs)
            if !((lhs[idx] == rhs[idx]) || (lhs[idx] ≈ rhs[idx]))
                push!(bad, idx)
            end
        end

        if isempty(bad)
            println("Passed left eigenvector check for block $k")
        else
            println("Failed left eigenvector check for block $k")
            println("Number of bad positions: ", length(bad))
            println("First bad positions:")
            display(bad[1:min(end, 10)])
            println("lhs values:")
            display(lhs[bad[1:min(end, 10)]])
            println("rhs values:")
            display(rhs[bad[1:min(end, 10)]])
            error("Left eigenvector check failed for block $k")
        end
    end
end

