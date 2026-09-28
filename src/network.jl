# Generic network I/O and max-plus analysis helpers.
#
# Matrix convention throughout the package:
#     A[target, source]
# so an edge source -> target is stored in A[target, source].

"""
    load_adjacency_csv(path; label_column=1, require_matching_labels=true)

Load a square weighted adjacency matrix from CSV.

By default, the first column contains row/node labels and the remaining column
headers contain the same labels. Columns are automatically reordered to the row
label order when necessary.

Returns a named tuple `(matrix, labels, index_by_label)`.
"""
function load_adjacency_csv(
    path::AbstractString;
    label_column::Union{Nothing,Int}=1,
    require_matching_labels::Bool=true,
)
    isfile(path) || throw(ArgumentError("File not found: $path"))

    df = CSV.read(path, DataFrame; normalizenames=false)

    if label_column === nothing
        nrow(df) == ncol(df) || throw(
            ArgumentError(
                "Adjacency CSV without a label column must be square; got " *
                "$(nrow(df)) rows and $(ncol(df)) columns."
            )
        )

        labels = String.(names(df))
        matrix_columns = collect(1:ncol(df))
    else
        1 <= label_column <= ncol(df) || throw(
            ArgumentError("label_column=$label_column is outside the CSV columns.")
        )

        labels = string.(df[:, label_column])
        matrix_columns = [j for j in 1:ncol(df) if j != label_column]

        length(labels) == length(matrix_columns) || throw(
            ArgumentError(
                "Adjacency CSV must contain one label column plus an n×n matrix; " *
                "got $(length(labels)) rows and $(length(matrix_columns)) matrix columns."
            )
        )

        column_labels = String.(names(df)[matrix_columns])

        if require_matching_labels
            length(unique(labels)) == length(labels) || throw(
                ArgumentError("Row labels must be unique.")
            )
            Set(labels) == Set(column_labels) || throw(
                ArgumentError(
                    "Row labels and matrix-column labels differ. " *
                    "Use the same node labels on rows and columns."
                )
            )

            position = Dict(label => j for (j, label) in enumerate(column_labels))
            matrix_columns = matrix_columns[[position[label] for label in labels]]
        end
    end

    n = length(labels)
    A = Matrix{Float64}(undef, n, n)

    for i in 1:n, j in 1:n
        value = df[i, matrix_columns[j]]
        ismissing(value) && throw(
            ArgumentError("Missing value at matrix row $i, column $j in $path")
        )

        if value isa Number
            A[i, j] = Float64(value)
        else
            parsed = tryparse(Float64, strip(string(value)))
            parsed === nothing && throw(
                ArgumentError("Non-numeric matrix value '$value' at row $i, column $j.")
            )
            A[i, j] = parsed
        end
    end

    index_by_label = Dict(label => i for (i, label) in enumerate(labels))
    return (matrix=A, labels=labels, index_by_label=index_by_label)
end


"""
    load_edgelist_csv(path; source_col=:source, target_col=:target,
                      weight_col=:weight, directed=true, labels=nothing)

Load a weighted edge-list CSV and convert it to a sparse adjacency matrix using
the package `A[target, source]` convention.

If `weight_col=nothing`, every edge receives weight `1.0`. Duplicate edges are
summed. If `labels` is omitted, node labels are inferred in first-appearance
order from the source and target columns.

Returns a named tuple `(matrix, labels, index_by_label)`.
"""
function load_edgelist_csv(
    path::AbstractString;
    source_col=:source,
    target_col=:target,
    weight_col=:weight,
    directed::Bool=true,
    labels=nothing,
)
    isfile(path) || throw(ArgumentError("File not found: $path"))

    df = CSV.read(path, DataFrame; normalizenames=false)
    available = Set(Symbol.(names(df)))
    s_col = source_col isa Symbol ? source_col : Symbol(source_col)
    t_col = target_col isa Symbol ? target_col : Symbol(target_col)

    s_col in available || throw(
        ArgumentError("Source column '$source_col' not found. Available columns: $(names(df))")
    )
    t_col in available || throw(
        ArgumentError("Target column '$target_col' not found. Available columns: $(names(df))")
    )

    w_col = weight_col === nothing ? nothing : (weight_col isa Symbol ? weight_col : Symbol(weight_col))
    if w_col !== nothing && !(w_col in available)
        throw(
            ArgumentError("Weight column '$weight_col' not found. Available columns: $(names(df))")
        )
    end

    sources_raw = df[!, s_col]
    targets_raw = df[!, t_col]

    any(ismissing, sources_raw) && throw(ArgumentError("Source column contains missing values."))
    any(ismissing, targets_raw) && throw(ArgumentError("Target column contains missing values."))

    sources = string.(sources_raw)
    targets = string.(targets_raw)

    node_labels = if labels === nothing
        unique(vcat(sources, targets))
    else
        string.(collect(labels))
    end

    length(unique(node_labels)) == length(node_labels) || throw(
        ArgumentError("Node labels must be unique.")
    )

    index_by_label = Dict(label => i for (i, label) in enumerate(node_labels))

    unknown = setdiff(Set(vcat(sources, targets)), Set(node_labels))
    isempty(unknown) || throw(
        ArgumentError("Edge list contains labels not present in `labels`: $(collect(unknown))")
    )

    weights = Vector{Float64}(undef, nrow(df))
    if w_col === nothing
        fill!(weights, 1.0)
    else
        for i in 1:nrow(df)
            value = df[i, w_col]
            ismissing(value) && throw(ArgumentError("Weight column contains a missing value at row $i."))

            if value isa Number
                weights[i] = Float64(value)
            else
                parsed = tryparse(Float64, strip(string(value)))
                parsed === nothing && throw(
                    ArgumentError("Non-numeric weight '$value' at row $i.")
                )
                weights[i] = parsed
            end
        end
    end

    rows = Int[]
    cols = Int[]
    vals = Float64[]

    for i in eachindex(sources)
        source = index_by_label[sources[i]]
        target = index_by_label[targets[i]]

        push!(rows, target)
        push!(cols, source)
        push!(vals, weights[i])

        if !directed && source != target
            push!(rows, source)
            push!(cols, target)
            push!(vals, weights[i])
        end
    end

    n = length(node_labels)
    A = sparse(rows, cols, vals, n, n, +)

    return (matrix=A, labels=node_labels, index_by_label=index_by_label)
end


"""
    to_maxplus(A; zero_is_absent=true)

Convert a real-valued adjacency matrix to a dense max-plus matrix.

When `zero_is_absent=true` (default), numeric zeros denote missing edges and are
mapped to max-plus bottom (`-Inf`). Set it to `false` only for a dense matrix in
which zero is a genuine edge weight and absence is already encoded as `-Inf`.
"""
function to_maxplus(
    A::AbstractMatrix{<:Real};
    zero_is_absent::Bool=true,
)
    n, m = size(A)
    n == m || throw(ArgumentError("A must be square; got $(n)×$(m)."))

    B = fill(zero(MP), n, m)

    for i in 1:n, j in 1:m
        x = Float64(A[i, j])
        isnan(x) && throw(ArgumentError("NaN found at A[$i,$j]."))

        if x == -Inf || (zero_is_absent && iszero(x))
            continue
        end

        B[i, j] = MP(x)
    end

    return B
end

function to_maxplus(
    A::SparseMatrixCSC{<:Real};
    zero_is_absent::Bool=false,
)
    n, m = size(A)
    n == m || throw(ArgumentError("A must be square; got $(n)×$(m)."))

    B = fill(zero(MP), n, m)
    Iidx, Jidx, V = findnz(A)

    for k in eachindex(V)
        x = Float64(V[k])
        isnan(x) && throw(ArgumentError("NaN found at A[$(Iidx[k]),$(Jidx[k])]."))

        if x == -Inf || (zero_is_absent && iszero(x))
            continue
        end

        B[Iidx[k], Jidx[k]] = MP(x)
    end

    return B
end

# Already-max-plus matrices are copied without reinterpretation.
to_maxplus(A::AbstractMatrix{<:MP}; zero_is_absent::Bool=true) = Matrix(A)


"""
    network_summary(layers)

Return a compact DataFrame describing one or more adjacency matrices. `layers`
may be a NamedTuple, dictionary, or any object supported by `pairs`.
"""
function network_summary(layers)
    rows = DataFrame(
        Layer=String[],
        Nodes=Int[],
        Edges=Int[],
        Symmetric=Bool[],
        Density=Float64[],
    )

    for (name, A) in pairs(layers)
        n, m = size(A)
        n == m || throw(ArgumentError("Layer '$name' is not square: $(n)×$(m)."))
        edges = if A isa SparseMatrixCSC
            nnz(A)
        elseif eltype(A) <: MP
            count(x -> x != zero(MP), A)
        elseif any(x -> x == -Inf, A)
            count(x -> x != -Inf, A)
        else
            count(!iszero, A)
        end

        push!(rows, (
            string(name),
            n,
            edges,
            issymmetric(A),
            edges / (n * n),
        ))
    end

    return rows
end


"""
    check_eigenpairs(A, spectrum; side=:right)

Verify every eigenvector equation stored in a spectrum returned by `mpeigen` or
`lmpeigen`. For a left spectrum, pass `transpose(A)` as `A` because left
eigenvectors are represented internally as right eigenvectors of the transpose.
"""
function check_eigenpairs(
    A::AbstractMatrix{<:MP},
    spectrum;
    side::Symbol=:right,
)
    for k in eachindex(spectrum.values)
        λ = spectrum.values[k]
        V = spectrum.vectors[k]

        length(λ) == size(V, 2) || error(
            "$side spectrum block $k has $(length(λ)) eigenvalues but $(size(V,2)) vectors."
        )

        lhs = A * V
        rhs = V * diagm(λ)

        valid = all(
            idx -> lhs[idx] == rhs[idx] || lhs[idx] ≈ rhs[idx],
            eachindex(lhs),
        )

        valid || error("$side eigenvector check failed in block $k")
    end

    return true
end


"""
    analyse_network(A; labels=nothing, name="network", tol=1e-10,
                    zero_is_absent=nothing, compute_left=true,
                    validate=true, verbose=true)

Run the complete max-plus eigendecomposition for an arbitrary weighted network.
The input follows the convention `A[target, source]`.

Returns a named tuple containing the max-plus matrix, labels, right spectrum,
and (optionally) left spectrum.
"""
function analyse_network(
    A::AbstractMatrix;
    labels=nothing,
    name::AbstractString="network",
    tol::Float64=1e-10,
    zero_is_absent::Union{Nothing,Bool}=nothing,
    compute_left::Bool=true,
    validate::Bool=true,
    verbose::Bool=true,
)
    n, m = size(A)
    n == m || throw(ArgumentError("A must be square; got $(n)×$(m)."))

    node_labels = labels === nothing ? ["node_$i" for i in 1:n] : string.(collect(labels))
    length(node_labels) == n || throw(
        ArgumentError("labels has length $(length(node_labels)), but A is $(n)×$(n).")
    )
    length(unique(node_labels)) == n || throw(ArgumentError("labels must be unique."))

    zero_absent = zero_is_absent === nothing ? !(A isa SparseMatrixCSC) : zero_is_absent
    A_mp = to_maxplus(A; zero_is_absent=zero_absent)
    right = mpeigen(A_mp; tol=tol)

    symmetric = issymmetric(A_mp)
    left = if compute_left
        symmetric ? right : lmpeigen(A_mp; tol=tol)
    else
        nothing
    end

    if validate
        check_eigenpairs(A_mp, right; side=:right)
        if compute_left
            check_eigenpairs(Matrix(transpose(A_mp)), left; side=:left)
        end
    end

    if verbose
        println(name)
        println("  nodes: ", n)
        println("  symmetric: ", symmetric)
        println("  right spectrum blocks: ", length(right.values))
        compute_left && println("  left spectrum blocks: ", length(left.values))
        validate && println("  eigenvector checks: passed")
    end

    return (
        name=String(name),
        labels=node_labels,
        matrix=A_mp,
        right=right,
        left=left,
        symmetric=symmetric,
    )
end
