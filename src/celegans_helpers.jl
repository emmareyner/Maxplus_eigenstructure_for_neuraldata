module CElegansHelpers

using CSV
using DataFrames
using Graphs
using SparseArrays
using XLSX

using Maxplus_eigenvectors

export DEFAULT_COOK_FILE,
       DEFAULT_MONOAMINE_FILE,
       DEFAULT_NEUROPEPTIDE_FILE,
       DEFAULT_FUNCTIONAL_FILE,
       load_celegans_structural,
       load_celegans_functional,
       load_celegans_thesis_layers,
       scc_statistics,
       critical_classes_table,
       thesis_reproducibility_checks

const DEFAULT_COOK_FILE =
    "SI 5 Connectome adjacency matrices, corrected July 2020 (1).xlsx"
const DEFAULT_MONOAMINE_FILE =
    "08062023_monoamine_connectome.csv"
const DEFAULT_NEUROPEPTIDE_FILE =
    "01022024_neuropeptide_connectome_short_range_model.csv"
const DEFAULT_FUNCTIONAL_FILE =
    joinpath("functional", "functional_connectivity_all.csv")

# ---------------------------------------------------------------------------
# Input normalisation
# ---------------------------------------------------------------------------

function _norm_neuron_name(x)
    (ismissing(x) || x === nothing) && return ""

    s = uppercase(strip(string(x)))
    isempty(s) && return ""

    # Normalise numbered motor-neuron names (e.g. DA01 -> DA1).
    m = match(r"^(DA|DB|DD|VA|VB|VD|AS|VC)0*([0-9]+)$", s)
    m === nothing && return s

    return m.captures[1] * string(parse(Int, m.captures[2]))
end

function _tofloat(x)
    (ismissing(x) || x === nothing) && return 0.0
    x isa Number && return Float64(x)

    s = strip(string(x))
    isempty(s) && return 0.0

    value = tryparse(Float64, s)
    return value === nothing ? 0.0 : value
end

function _read_square_csv(path)
    isfile(path) || throw(ArgumentError("File not found: $path"))

    df = CSV.read(path, DataFrame; normalizenames=false)
    row_ids = _norm_neuron_name.(df[:, 1])
    col_ids = _norm_neuron_name.(names(df)[2:end])

    nrow(df) == length(col_ids) || throw(
        ArgumentError("Expected a square adjacency CSV in $path")
    )
    Set(row_ids) == Set(col_ids) || throw(
        ArgumentError("Row and column neuron labels differ in $path")
    )

    col_index = Dict(label => j for (j, label) in enumerate(col_ids))
    perm = [col_index[label] for label in row_ids]

    raw = Matrix(df[:, 2:end])[:, perm]
    W = zeros(Float64, length(row_ids), length(row_ids))

    for i in axes(raw, 1), j in axes(raw, 2)
        W[i, j] = _tofloat(raw[i, j])
    end

    return row_ids, W
end

function _read_sheet(path, sheet_name)
    isfile(path) || throw(ArgumentError("File not found: $path"))

    XLSX.openxlsx(path) do xf
        sheet_name in XLSX.sheetnames(xf) || throw(
            ArgumentError("Sheet '$sheet_name' not found in $path")
        )
        return xf[sheet_name][:]
    end
end

"""
    _cook_layer(path, sheet_name, neuron_ids)

Extract the neuron-to-neuron block from a Cook et al. worksheet. The source
file stores source/presynaptic neurons in rows and target/postsynaptic neurons
in columns, so the returned matrix is initially `W[source, target]`.
"""
function _cook_layer(path, sheet_name, neuron_ids)
    raw = _read_sheet(path, sheet_name)
    valid = Set(neuron_ids)

    row_scores = [
        count(x -> _norm_neuron_name(x) in valid, raw[r, :])
        for r in axes(raw, 1)
    ]
    col_scores = [
        count(x -> _norm_neuron_name(x) in valid, raw[:, c])
        for c in axes(raw, 2)
    ]

    header_row = argmax(row_scores)
    label_col = argmax(col_scores)

    row_scores[header_row] > 0 || error(
        "Could not detect neuron header row in '$sheet_name'."
    )
    col_scores[label_col] > 0 || error(
        "Could not detect neuron-label column in '$sheet_name'."
    )

    row_positions = Dict{String,Int}()
    for r in axes(raw, 1)
        label = _norm_neuron_name(raw[r, label_col])
        label in valid && (row_positions[label] = r)
    end

    col_positions = Dict{String,Int}()
    for c in axes(raw, 2)
        label = _norm_neuron_name(raw[header_row, c])
        label in valid && (col_positions[label] = c)
    end

    idx = Dict(label => i for (i, label) in enumerate(neuron_ids))
    W = zeros(Float64, length(neuron_ids), length(neuron_ids))

    for (source, r) in row_positions
        i = idx[source]
        for (target, c) in col_positions
            j = idx[target]
            W[i, j] = _tofloat(raw[r, c])
        end
    end

    metadata = (
        header_row=header_row,
        label_col=label_col,
        source_neurons_found=length(row_positions),
        target_neurons_found=length(col_positions),
    )

    return W, metadata
end

# ---------------------------------------------------------------------------
# Thesis data loaders
# ---------------------------------------------------------------------------

"""
    load_celegans_structural(connectome_dir; kwargs...)

Load the four structural/extrasynaptic connectomes used in the thesis and
harmonise them to a common 302-neuron order.

All returned matrices follow the package convention `A[target, source]`.
"""
function load_celegans_structural(
    connectome_dir;
    cook_file=DEFAULT_COOK_FILE,
    monoamine_file=DEFAULT_MONOAMINE_FILE,
    neuropeptide_file=DEFAULT_NEUROPEPTIDE_FILE,
    chemical_sheet="hermaphrodite chemical",
    gap_sheet="hermaphrodite gap jn symmetric",
)
    cook_path = joinpath(connectome_dir, cook_file)
    monoamine_path = joinpath(connectome_dir, monoamine_file)
    neuropeptide_path = joinpath(connectome_dir, neuropeptide_file)

    for path in (cook_path, monoamine_path, neuropeptide_path)
        isfile(path) || throw(ArgumentError(
            "Required C. elegans data file not found: $path"
        ))
    end

    # The neuropeptide file supplies the canonical 302-neuron ordering used in
    # the original thesis-analysis workflow.
    labels, neuropeptide_source_target = _read_square_csv(neuropeptide_path)

    monoamine_labels, monoamine_source_target = _read_square_csv(monoamine_path)
    Set(monoamine_labels) == Set(labels) || error(
        "Monoamine and neuropeptide files do not contain the same neuron set."
    )

    ma_index = Dict(label => i for (i, label) in enumerate(monoamine_labels))
    ma_perm = [ma_index[label] for label in labels]
    monoamine_source_target = monoamine_source_target[ma_perm, ma_perm]

    synaptic_source_target, syn_meta = _cook_layer(
        cook_path, chemical_sheet, labels
    )
    gap_source_target, gap_meta = _cook_layer(
        cook_path, gap_sheet, labels
    )

    # Source datasets use rows=source, columns=target. The mathematical package
    # uses A[target, source], hence the transpose here.
    layers = (
        synaptic=sparse(transpose(synaptic_source_target)),
        gap_junction=sparse(transpose(gap_source_target)),
        monoamine=sparse(transpose(monoamine_source_target)),
        neuropeptide=sparse(transpose(neuropeptide_source_target)),
    )

    return (
        labels=labels,
        index_by_label=Dict(label => i for (i, label) in enumerate(labels)),
        layers=layers,
        metadata=(
            cook_file=cook_path,
            monoamine_file=monoamine_path,
            neuropeptide_file=neuropeptide_path,
            synaptic_detection=syn_meta,
            gap_detection=gap_meta,
        ),
    )
end

"""
    load_celegans_functional(path)

Load the functional signal-propagation matrix used in the thesis. Missing/NaN
entries remain absent connections. The CSV is already oriented according to
`A[target, source]`, matching the original analysis notebook.
"""
function load_celegans_functional(path)
    isfile(path) || throw(ArgumentError(
        "Required C. elegans functional file not found: $path"
    ))

    df = CSV.read(path, DataFrame; normalizenames=false)
    labels = string.(df[:, 1])
    column_labels = string.(names(df)[2:end])

    length(labels) == length(column_labels) || throw(
        ArgumentError("Functional connectivity CSV must be square.")
    )
    Set(labels) == Set(column_labels) || throw(
        ArgumentError("Functional row and column neuron labels differ.")
    )

    col_index = Dict(label => j for (j, label) in enumerate(column_labels))
    perm = [col_index[label] for label in labels]
    raw = Matrix(df[:, 2:end])[:, perm]

    rows = Int[]
    cols = Int[]
    vals = Float64[]

    # Preserve the original notebook representation: every finite entry is
    # stored in the sparse matrix, including explicit zeros. The subsequent
    # Max-plus conversion is called with zero_is_absent=true, so zero-valued
    # entries do not become graph edges.
    for i in axes(raw, 1), j in axes(raw, 2)
        value = raw[i, j]
        if ismissing(value)
            continue
        end

        x = value isa Number ? Float64(value) : tryparse(Float64, strip(string(value)))
        x === nothing && continue
        isfinite(x) || continue

        push!(rows, i)
        push!(cols, j)
        push!(vals, x)
    end

    A = sparse(rows, cols, vals, length(labels), length(labels))

    return (
        matrix=A,
        labels=labels,
        index_by_label=Dict(label => i for (i, label) in enumerate(labels)),
    )
end

"""
    load_celegans_thesis_layers(connectome_dir; functional_file=DEFAULT_FUNCTIONAL_FILE)

Load the five network layers analysed in the thesis. Each entry contains a
matrix in `A[target, source]` orientation and its matching labels.
"""
function load_celegans_thesis_layers(
    connectome_dir;
    functional_file=DEFAULT_FUNCTIONAL_FILE,
)
    structural = load_celegans_structural(connectome_dir)
    functional = load_celegans_functional(joinpath(connectome_dir, functional_file))

    layers = Dict{Symbol,Any}(
        :synaptic => (matrix=structural.layers.synaptic, labels=structural.labels),
        :gap_junction => (matrix=structural.layers.gap_junction, labels=structural.labels),
        :monoamine => (matrix=structural.layers.monoamine, labels=structural.labels),
        :neuropeptide => (matrix=structural.layers.neuropeptide, labels=structural.labels),
        :functional => (matrix=functional.matrix, labels=functional.labels),
    )

    return (
        layers=layers,
        structural=structural,
        functional=functional,
    )
end

# ---------------------------------------------------------------------------
# Compact reproducibility helpers
# ---------------------------------------------------------------------------

function _simple_digraph(A)
    n, m = size(A)
    n == m || throw(ArgumentError("A must be square."))
    g = SimpleDiGraph(n)

    if A isa SparseMatrixCSC
        I, J, V = findnz(A)
        for k in eachindex(V)
            iszero(V[k]) && continue
            # A[target, source] -> source -> target
            add_edge!(g, J[k], I[k])
        end
    else
        for target in 1:n, source in 1:n
            value = A[target, source]
            iszero(value) && continue
            add_edge!(g, source, target)
        end
    end

    return g
end

"""
    scc_statistics(A)

Return the SCC count, largest SCC size, number of non-singleton SCCs, and the
sorted SCC-size vector for a network in `A[target, source]` orientation.
"""
function scc_statistics(A)
    sccs = strongly_connected_components(_simple_digraph(A))
    sizes = sort(length.(sccs); rev=true)

    return (
        count=length(sccs),
        largest=isempty(sizes) ? 0 : first(sizes),
        non_singleton=count(>(1), sizes),
        sizes=sizes,
        components=sccs,
    )
end

"""
    critical_classes_table(analysis)

Return one row per critical class across all finite right projective
eigendirections in an `analyse_network` result.
"""
function critical_classes_table(analysis)
    rows = DataFrame(
        Block=Int[],
        Direction=Int[],
        Lambda=Float64[],
        CriticalClass=Int[],
        CriticalLabels=String[],
        RepresentativeCycle=String[],
        CycleMean=Float64[],
    )

    spectrum = finite_right_spectrum_summary(analysis.right)

    for spectrum_row in eachrow(spectrum)
        block = Int(spectrum_row.Block)
        direction = Int(spectrum_row.Direction)
        cs = critical_structure(
            analysis.matrix,
            analysis.right;
            block=block,
            direction=direction,
        )
        summary = critical_structure_summary(cs, analysis.labels)

        for row in eachrow(summary)
            push!(rows, (
                block,
                direction,
                Float64(row.Lambda),
                Int(row.CriticalClass),
                String(row.CriticalLabels),
                String(row.RepresentativeCycle),
                Float64(row.CycleMean),
            ))
        end
    end

    return rows
end

"""
    thesis_reproducibility_checks(analyses, layer_inputs)

Run compact assertions for the headline structures reported in the thesis.
The checks intentionally target stable biological/mathematical results rather
than every printed table or plotting detail.
"""
function thesis_reproducibility_checks(analyses, layer_inputs)
    expected_scc = Dict(
        :synaptic => (count=11, largest=275),
        :gap_junction => (count=6, largest=296),
        :monoamine => (count=286, largest=17),
        :neuropeptide => (count=17, largest=286),
        :functional => (count=42, largest=147),
    )

    for layer in keys(expected_scc)
        stats = scc_statistics(layer_inputs[layer].matrix)
        expected = expected_scc[layer]
        @assert(stats.count == expected.count,
            "$layer: expected $(expected.count) SCCs, found $(stats.count).")
        @assert(stats.largest == expected.largest,
            "$layer: expected largest SCC $(expected.largest), found $(stats.largest).")
    end

    expected_lambda = Dict(
        :synaptic => [15.0, 37.0],
        :gap_junction => [1.0, 401.0],
        :monoamine => [2.5],
        :neuropeptide => [14.0],
        :functional => [0.921569],
    )

    for layer in keys(expected_lambda)
        found = sort(unique([
            round(x; digits=6)
            for x in finite_right_spectrum_summary(analyses[layer].right).Lambda
        ]))
        expected = sort(expected_lambda[layer])
        @assert(found == expected,
            "$layer: expected finite right eigenvalues $expected, found $found.")
    end

    critical = Dict(
        layer => critical_classes_table(analyses[layer])
        for layer in keys(analyses)
    )

    @assert any(occursin.("RIAR", critical[:synaptic].CriticalLabels))
    @assert any(occursin.("SMDVL", critical[:synaptic].CriticalLabels))
    @assert any(==(14.0), critical[:neuropeptide].Lambda)
    @assert any(==(2.5), critical[:monoamine].Lambda)
    @assert any(occursin.("PQR", critical[:neuropeptide].CriticalLabels))
    @assert any(occursin.("RIH", critical[:functional].CriticalLabels))
    @assert any(occursin.("SIBVR", critical[:functional].CriticalLabels))

    println("✓ Thesis headline SCC, eigenvalue, and critical-class checks passed.")
    return true
end

end # module CElegansHelpers
