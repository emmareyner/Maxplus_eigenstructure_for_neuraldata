# Critical max-plus structure extraction and Graphviz visualisation.
#
# Include this file INSIDE module Maxplus_eigenvectors after validation.jl.
#
# It contains both:
#   1) mathematical extraction of recurrent/critical structure;
#   2) Graphviz SVG visualisation.
#
# Package matrix convention:
#     A[receiving, sending]
# so a network arrow source -> target corresponds to A[target, source].

const _CRITICAL_CLASS_COLOURS = [
    "#c62828",
    "#1565c0",
    "#2e7d32",
    "#6a1b9a",
    "#ef6c00",
    "#00838f",
    "#ad1457",
    "#455a64",
]

_dot_escape(s) = replace(
    string(s),
    "\\" => "\\\\",
    "\"" => "\\\"",
    "\n" => "\\n",
)

_graphviz_colour(cid::Int) =
    _CRITICAL_CLASS_COLOURS[mod1(cid, length(_CRITICAL_CLASS_COLOURS))]


# ---------------------------------------------------------------------------
# Select the same projective eigendirection used by validation/reporting
# ---------------------------------------------------------------------------

"""
    eigen_direction(sp, block, direction; tol=1e-10)

Select the projective eigendirection reported by
`finite_right_spectrum_summary` / `finite_right_eigenvector_support`.

Returns:
- eigenvalue;
- representative eigenvector;
- support;
- projectively equivalent candidate columns;
- projective multiplicity.

`multiplicity` is the number of candidate columns collapsed into this
projective direction. It is not automatically the number of critical
cycles or critical classes.
"""
function eigen_direction(sp, block::Int, direction::Int; tol=1e-10)
    1 <= block <= length(sp.values) ||
        throw(BoundsError(sp.values, block))

    vals = sp.values[block]
    V = sp.vectors[block]
    groups = projective_eigenvector_groups(vals, V; tol=tol)

    1 <= direction <= length(groups) ||
        throw(BoundsError(groups, direction))

    group = groups[direction]
    j = first(group)

    λ = vals[j]
    v = V[:, j]

    # Partition the eigendirection coordinates.
    vf = plustimes.(v)

    finite_support = findall(isfinite, vf)
    top_support = findall(==(Inf), vf)
    bottom_support = findall(==(-Inf), vf)

    # Preserve the previous algebraic support definition:
    # everything except bottom (-Inf).
    support = sort(vcat(finite_support, top_support))

    return (
        block = block,
        direction = direction,
        representative_column = j,
        group = group,
        multiplicity = length(group),
        lambda_mp = λ,
        lambda = plustimes(λ),
        vector = v,

        # Previous non-bottom support.
        support = support,

        # Explicit coordinate classes.
        finite_support = finite_support,
        top_support = top_support,
        bottom_support = bottom_support,
    )
end


# ---------------------------------------------------------------------------
# Graph helpers
# ---------------------------------------------------------------------------

function _finite_edges(A::AbstractMatrix{<:MP}, nodes)
    allowed = Set(Int.(collect(nodes)))
    edges = Tuple{Int,Int,Float64}[]

    # A[target, source] -> source -> target
    for target in allowed, source in allowed
        a = A[target, source]

        if a != zero(MP)
            push!(edges, (source, target, plustimes(a)))
        end
    end

    sort!(edges; by=e -> (e[1], e[2]))
    return edges
end


function _graph_from_edges(n::Int, edges)
    g = SimpleDiGraph(n)

    for (source, target, _) in edges
        Graphs.add_edge!(g, source, target)
    end

    return g
end


function _cyclic_classes(g::SimpleDiGraph, allowed_nodes)
    allowed = Set(Int.(collect(allowed_nodes)))
    classes = Vector{Vector{Int}}()

    for scc in strongly_connected_components(g)
        cls = sort([u for u in scc if u in allowed])
        isempty(cls) && continue

        if length(cls) > 1 || Graphs.has_edge(g, cls[1], cls[1])
            push!(classes, cls)
        end
    end

    sort!(classes; by=cls -> (minimum(cls), length(cls)))
    return classes
end


function _class_map(classes)
    mapping = Dict{Int,Int}()

    for (cid, cls) in enumerate(classes), node in cls
        mapping[node] = cid
    end

    return mapping
end


function _edges_inside_classes(edges, classes)
    cmap = _class_map(classes)

    return [
        e for e in edges
        if haskey(cmap, e[1]) &&
           haskey(cmap, e[2]) &&
           cmap[e[1]] == cmap[e[2]]
    ]
end


# ---------------------------------------------------------------------------
# Saturation graph and critical graph
# ---------------------------------------------------------------------------

"""
    saturation_edges(A, λ, v; nodes, tol=1e-10)

Return tight edges of the right max-plus eigenvector equation

    A[target, source] + v[source] = λ + v[target].

Edges are returned in network direction:
    (source, target, weight)
"""
function saturation_edges(
    A::AbstractMatrix{<:MP},
    λ::MP,
    v::AbstractVector{<:MP};
    nodes=findall(i -> isfinite(plustimes(v[i])), eachindex(v)),
    tol=1e-10,
)
    allowed = Set(Int.(collect(nodes)))
    λf = plustimes(λ)

    isfinite(λf) ||
        throw(ArgumentError("saturation_edges requires a finite eigenvalue."))

    edges = Tuple{Int,Int,Float64}[]

    for target in allowed, source in allowed
        a = A[target, source]

        a == zero(MP) && continue

        af = plustimes(a)
        vs = plustimes(v[source])
        vt = plustimes(v[target])

        # Saturation through ordinary finite potentials is only meaningful
        # when edge weight and both vector coordinates are finite.
        isfinite(af) || continue
        isfinite(vs) || continue
        isfinite(vt) || continue

        lhs = af + vs
        rhs = λf + vt

        if isapprox(lhs, rhs; atol=tol, rtol=0.0)
            push!(edges, (source, target, af))
        end
    end

    sort!(edges; by=e -> (e[1], e[2]))
    return edges
end

function _restricted_path(
    g::SimpleDiGraph,
    start::Int,
    goal::Int,
    allowed::Set{Int},
)
    start == goal && return [start]

    queue = Int[start]
    head = 1
    parent = Dict{Int,Int}(start => 0)

    while head <= length(queue)
        u = queue[head]
        head += 1

        for v in Graphs.outneighbors(g, u)
            v in allowed || continue
            haskey(parent, v) && continue

            parent[v] = u

            if v == goal
                path = Int[goal]
                current = goal

                while parent[current] != 0
                    current = parent[current]
                    push!(path, current)
                end

                reverse!(path)
                return path
            end

            push!(queue, v)
        end
    end

    return Int[]
end


function _representative_cycle(g::SimpleDiGraph, cls::Vector{Int})
    allowed = Set(cls)

    # Prefer a self-loop if the class contains one.
    for u in cls
        if Graphs.has_edge(g, u, u)
            return [u]
        end
    end

    # In an SCC, any internal edge u -> v has a return path v -> ... -> u.
    for u in cls
        for v in Graphs.outneighbors(g, u)
            v in allowed || continue
            u == v && continue

            back = _restricted_path(g, v, u, allowed)
            isempty(back) && continue

            # u -> v -> ... -> u, omitting repeated terminal u.
            return vcat([u], back[1:end-1])
        end
    end

    error("Could not recover a representative cycle from a cyclic SCC.")
end


"""
    critical_cycle_mean(A, cycle)

Arithmetic mean of the edge weights around a directed cycle.

The cycle is specified in network direction:
    source1 -> source2 -> ... -> source1

The matrix is read using the package convention A[target, source].
"""
function critical_cycle_mean(A::AbstractMatrix{<:MP}, cycle)
    c = Int.(collect(cycle))
    isempty(c) && throw(ArgumentError("Cycle cannot be empty."))

    weights = Float64[]

    for k in eachindex(c)
        source = c[k]
        target = c[mod1(k + 1, length(c))]
        a = A[target, source]

        a == zero(MP) &&
            throw(ArgumentError("Missing edge in cycle: $source -> $target"))

        push!(weights, plustimes(a))
    end

    return sum(weights) / length(weights)
end


"""
    critical_structure(A, sp; block, direction, tol=1e-10)

Recover recurrent and critical structure for one reported right max-plus
eigendirection.

The routine does not enumerate all simple cycles.

It returns:
- `support_edges`: all finite edges inside the eigenvector support;
- `cycle_classes`: SCCs of the support that contain cycles;
- `cycle_edges`: union of all support edges that belong to at least one cycle;
- `saturation_edges`: all tight edges of A*v = λ*v;
- `critical_classes`: cyclic SCCs of the saturation graph;
- `critical_edges`: union of all edges in the critical graph;
- `representative_cycles`: one explicit critical cycle per critical class;
- `representative_cycle_means`: their cycle means.

Every directed cycle contained in the critical graph has mean λ, up to
numerical tolerance, because the eigenvector potentials telescope around
a closed path.
"""
function critical_structure(
    A::AbstractMatrix{<:MP},
    sp;
    block::Int,
    direction::Int,
    tol=1e-10,
)
    ed = eigen_direction(sp, block, direction; tol=tol)

    # For recurrent/critical reconstruction, use only coordinates
    # with finite max-plus potential.
    finite_nodes = ed.finite_support

    isempty(finite_nodes) &&
        throw(
            ArgumentError(
                "The selected eigendirection has no finite coordinates."
            )
        )

    # All recurrent structure inside the finite support.
    support_edges = _finite_edges(A, finite_nodes)
    support_graph = _graph_from_edges(size(A, 1), support_edges)

    cycle_classes = _cyclic_classes(support_graph, finite_nodes)
    cycle_edges = _edges_inside_classes(support_edges, cycle_classes)

    # Tight/saturated structure for this eigendirection.
    sat_edges = saturation_edges(
        A,
        ed.lambda_mp,
        ed.vector;
        nodes=finite_nodes,
        tol=tol,
    )

    sat_graph = _graph_from_edges(size(A, 1), sat_edges)

    critical_classes = _cyclic_classes(sat_graph, finite_nodes)
    critical_edges = _edges_inside_classes(sat_edges, critical_classes)

    critical_graph = _graph_from_edges(size(A, 1), critical_edges)

    representative_cycles = [
        _representative_cycle(critical_graph, cls)
        for cls in critical_classes
    ]

    representative_cycle_means = [
        critical_cycle_mean(A, cycle)
        for cycle in representative_cycles
    ]

    for μ in representative_cycle_means
        isapprox(μ, ed.lambda; atol=max(tol, 1e-9), rtol=0.0) ||
            error(
                "Recovered critical cycle mean $μ differs from λ=$(ed.lambda)."
            )
    end

    return merge(
        ed,
        (
            support_edges = support_edges,
            cycle_classes = cycle_classes,
            cycle_edges = cycle_edges,
            saturation_edges = sat_edges,
            critical_classes = critical_classes,
            critical_edges = critical_edges,
            representative_cycles = representative_cycles,
            representative_cycle_means = representative_cycle_means,
        ),
    )
end


"""
    critical_structure_summary(cs, labels)

Return one row per critical class, with one representative critical cycle
and its verified cycle mean.
"""
function critical_structure_summary(cs, labels)
    rows = DataFrame(
        CriticalClass = Int[],
        Lambda = Float64[],
        ClassSize = Int[],
        CriticalLabels = String[],
        RepresentativeCycle = String[],
        CycleLength = Int[],
        CycleMean = Float64[],
    )

    for cid in eachindex(cs.critical_classes)
        cls = cs.critical_classes[cid]
        cycle = cs.representative_cycles[cid]

        cycle_string =
            join(labels[cycle], " → ") *
            " → " *
            labels[first(cycle)]

        push!(rows, (
            cid,
            cs.lambda,
            length(cls),
            join(labels[cls], ", "),
            cycle_string,
            length(cycle),
            cs.representative_cycle_means[cid],
        ))
    end

    return rows
end


# ---------------------------------------------------------------------------
# Graphviz visualisation
# ---------------------------------------------------------------------------

"""
    visualize_critical_structure(
        A,
        sp,
        labels;
        block,
        direction,
        filename,
        tol=1e-10,
        engine=:dot,
        rankdir="TB",
        show_weights=true,
        show_noncritical_weights=true,
        weight_digits=3,
        keep_dot=true,
    )

Create a Graphviz SVG showing:

- all recurrent edges in the selected eigenvector support in grey;
- the critical graph in thicker coloured edges;
- each critical class in a different colour;
- edge weights, optionally for both critical and non-critical recurrent edges.

The layout is purely graphical. Distances between nodes have no network or
max-plus meaning.

Supported Graphviz engines:
    :dot, :neato, :fdp, :sfdp

Returns:
    (
        svgfile,
        dotfile,
        structure,
        summary,
    )
"""
function visualize_critical_structure(
    A::AbstractMatrix{<:MP},
    sp,
    labels;
    block::Int,
    direction::Int,
    filename::AbstractString,
    tol=1e-10,
    engine::Symbol=:dot,
    rankdir::AbstractString="TB",
    show_weights::Bool=true,
    show_noncritical_weights::Bool=true,
    weight_digits::Int=3,
    keep_dot::Bool=true,
)
    engine in (:dot, :neato, :fdp, :sfdp) ||
        throw(
            ArgumentError(
                "engine must be :dot, :neato, :fdp, or :sfdp"
            )
        )

    cs = critical_structure(
        A,
        sp;
        block=block,
        direction=direction,
        tol=tol,
    )

    recurrent_edges = cs.cycle_edges

    isempty(recurrent_edges) &&
        throw(
            ArgumentError(
                "No recurrent edges found in this eigendirection support."
            )
        )

    recurrent_nodes = sort(
        unique(
            vcat(
                [e[1] for e in recurrent_edges],
                [e[2] for e in recurrent_edges],
            )
        )
    )

    critical_edge_set = Set(
        (source, target)
        for (source, target, _) in cs.critical_edges
    )

    critical_class_by_node = Dict{Int,Int}()

    for (cid, cls) in enumerate(cs.critical_classes), node in cls
        critical_class_by_node[node] = cid
    end

    # Normalise output names.
    base, ext = splitext(String(filename))

    if lowercase(ext) in (".dot", ".svg")
        stem = base
    else
        stem = String(filename)
    end

    dotfile = stem * ".dot"
    svgfile = stem * ".svg"

    mkpath(dirname(abspath(dotfile)))

    open(dotfile, "w") do io
        println(io, "digraph G {")

        if engine == :dot
            println(
                io,
                """  graph [
                    layout=dot,
                    rankdir="$rankdir",
                    overlap=false,
                    splines=true,
                    outputorder=edgesfirst,
                    forcelabels=true,
                    nodesep=0.50,
                    ranksep=0.85,
                    pad=0.25,
                    bgcolor="white"
                ];"""
            )
        else
            println(
                io,
                """  graph [
                    layout=$(String(engine)),
                    overlap=prism,
                    splines=true,
                    outputorder=edgesfirst,
                    forcelabels=true,
                    sep="+20",
                    pad=0.25,
                    bgcolor="white"
                ];"""
            )
        end

        println(
            io,
            """  node [
                shape=ellipse,
                fontsize=11,
                fontname="Helvetica",
                margin="0.08,0.05",
                style="filled",
                fillcolor="white",
                color="#8a8f94",
                penwidth=1.0
            ];"""
        )

        println(
            io,
            """  edge [
                fontsize=8,
                fontname="Helvetica",
                arrowsize=0.65
            ];"""
        )

        # Nodes
        for node in recurrent_nodes
            name = _dot_escape(labels[node])

            if haskey(critical_class_by_node, node)
                cid = critical_class_by_node[node]
                colour = _graphviz_colour(cid)

                println(
                    io,
                    """  n$node [
                        label="$name",
                        color="$colour",
                        penwidth=2.6
                    ];"""
                )
            else
                println(io, """  n$node [label="$name"];""")
            end
        end

        # Recurrent edges
        for (source, target, weight) in recurrent_edges
            iscritical = (source, target) in critical_edge_set

            if iscritical
                cid = critical_class_by_node[source]
                colour = _graphviz_colour(cid)

                attrs = String[
                    "color=\"$colour\"",
                    "fontcolor=\"$colour\"",
                    "penwidth=3.0",
                ]

                if show_weights
                    push!(
                        attrs,
                        "label=\"$(round(weight, digits=weight_digits))\""
                    )
                end

                println(
                    io,
                    "  n$source -> n$target [" *
                    join(attrs, ", ") *
                    "];"
                )
            else
                attrs = String[
                    "color=\"#b7b9bc\"",
                    "fontcolor=\"#73777b\"",
                    "penwidth=0.8",
                ]

                if show_weights && show_noncritical_weights
                    push!(
                        attrs,
                        "label=\"$(round(weight, digits=weight_digits))\""
                    )
                end

                println(
                    io,
                    "  n$source -> n$target [" *
                    join(attrs, ", ") *
                    "];"
                )
            end
        end

        println(io, "}")
    end

    engine_cmd = String(engine)
    Sys.which(engine_cmd) === nothing && throw(
        ArgumentError("Graphviz executable '$engine_cmd' was not found on PATH.")
    )
    run(Cmd([engine_cmd, "-Tsvg", dotfile, "-o", svgfile]))

    summary = critical_structure_summary(cs, labels)

    println("λ = ", cs.lambda)
    println("Projective multiplicity = ", cs.multiplicity)
    println("Non-bottom support size = ", length(cs.support))
    println("Finite support size = ", length(cs.finite_support))
    println("+Inf coordinates = ", length(cs.top_support))
    println("-Inf coordinates = ", length(cs.bottom_support))
    println("Recurrent nodes shown = ", length(recurrent_nodes))
    println("Recurrent edges shown = ", length(recurrent_edges))
    println("Critical classes = ", length(cs.critical_classes))
    println("Critical edges = ", length(cs.critical_edges))
    println("SVG = ", svgfile)

    if !keep_dot
        rm(dotfile; force=true)
    end

    return (
        svgfile = svgfile,
        dotfile = keep_dot ? dotfile : nothing,
        structure = cs,
        summary = summary,
    )
end
