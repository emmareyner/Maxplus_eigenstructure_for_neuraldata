# Integrated SCC + critical-class analysis and condensation-graph visualisation.
#
# Include this file INSIDE module Maxplus_eigenvectors after:
#   validation.jl
#   critical_structure.jl
#
# Matrix convention:
#   A[target, source]
# so graph arrows are drawn source -> target.

"""
    analyse_scc_structure(A, sp, labels; block, direction, tol=1e-10)

Compute the SCC decomposition and condensation DAG for the full network,
and automatically locate the selected eigendirection support and critical
classes inside those SCCs.

Returns a named tuple containing the graph, SCCs, condensation graph,
critical SCCs, and a summary DataFrame.
"""
function analyse_scc_structure(
    A::AbstractMatrix{<:MP},
    sp,
    labels;
    block::Int,
    direction::Int,
    tol=1e-10,
)
    n = size(A, 1)
    size(A, 2) == n || throw(ArgumentError("A must be square."))
    length(labels) == n || throw(
        ArgumentError("labels has length $(length(labels)), but A is $(n)×$(n).")
    )

    ed = eigen_direction(
        sp,
        block,
        direction;
        tol=tol,
    )

    cs = critical_structure(
        A,
        sp;
        block=block,
        direction=direction,
        tol=tol,
    )

    # Build graph in network direction: source -> target.
    g = SimpleDiGraph(n)

    for target in 1:n, source in 1:n
        if A[target, source] != zero(MP)
            Graphs.add_edge!(g, source, target)
        end
    end

    sccs = strongly_connected_components(g)
    cdg = Graphs.condensation(g, sccs)

    # Node -> SCC map.
    scc_of_node = zeros(Int, n)

    for (sid, nodes) in enumerate(sccs)
        for node in nodes
            scc_of_node[node] = sid
        end
    end

    support_set = Set(ed.support)

    # Critical class -> containing SCC.
    critical_class_scc = Dict{Int,Int}()

    for (cid, cls) in enumerate(cs.critical_classes)
        sids = unique(scc_of_node[node] for node in cls)

        length(sids) == 1 || error(
            "Critical class $cid spans more than one SCC: $sids"
        )

        critical_class_scc[cid] = first(sids)
    end

    critical_sccs = sort(unique(collect(values(critical_class_scc))))

    rows = DataFrame(
        SCC = Int[],
        Size = Int[],
        SupportSize = Int[],
        SupportFraction = Float64[],
        InDegree = Int[],
        OutDegree = Int[],
        Critical = Bool[],
        CriticalClasses = String[],
        Labels = String[],
    )

    for (sid, nodes) in enumerate(sccs)
        support_nodes = [node for node in nodes if node in support_set]

        class_ids = sort([
            cid
            for (cid, mapped_sid) in critical_class_scc
            if mapped_sid == sid
        ])

        push!(
            rows,
            (
                sid,
                length(nodes),
                length(support_nodes),
                length(support_nodes) / length(nodes),
                Graphs.indegree(cdg, sid),
                Graphs.outdegree(cdg, sid),
                !isempty(class_ids),
                isempty(class_ids) ? "" : join(class_ids, ", "),
                join(labels[nodes], ", "),
            )
        )
    end

    sort!(rows, :Size, rev=true)

    return (
        block = block,
        direction = direction,
        lambda = ed.lambda,
        multiplicity = ed.multiplicity,
        support = ed.support,
        graph = g,
        sccs = sccs,
        scc_of_node = scc_of_node,
        condensation = cdg,
        critical_structure = cs,
        critical_class_scc = critical_class_scc,
        critical_sccs = critical_sccs,
        summary = rows,
    )
end


"""
    visualize_scc_structure(A, sp, labels; block, direction, filename, ...)

Create a Graphviz SVG of the condensation DAG.

Visual encoding:
- white node: SCC outside the selected eigenvector support;
- light-grey node: SCC contains at least one support node;
- red border: SCC contains a critical class;
- label: SCC id, total size, and number of support nodes.

No critical nodes need to be supplied manually.
"""
function visualize_scc_structure(
    A::AbstractMatrix{<:MP},
    sp,
    labels;
    block::Int,
    direction::Int,
    filename::AbstractString,
    tol=1e-10,
    engine::Symbol=:dot,
    rankdir::AbstractString="LR",
    keep_dot::Bool=true,
)
    engine in (:dot, :neato, :fdp, :sfdp) ||
        throw(ArgumentError("engine must be :dot, :neato, :fdp, or :sfdp"))

    result = analyse_scc_structure(
        A,
        sp,
        labels;
        block=block,
        direction=direction,
        tol=tol,
    )

    cdg = result.condensation
    rows = result.summary

    size_by_scc = Dict(row.SCC => row.Size for row in eachrow(rows))
    support_by_scc = Dict(row.SCC => row.SupportSize for row in eachrow(rows))
    critical_by_scc = Dict(row.SCC => row.Critical for row in eachrow(rows))
    classes_by_scc = Dict(row.SCC => row.CriticalClasses for row in eachrow(rows))
    labels_by_scc = Dict(row.SCC => row.Labels for row in eachrow(rows))

    base, ext = splitext(String(filename))
    stem = lowercase(ext) in (".dot", ".svg") ? base : String(filename)

    dotfile = stem * ".dot"
    svgfile = stem * ".svg"

    mkpath(dirname(abspath(dotfile)))

    open(dotfile, "w") do io
        println(io, "digraph G {")
        println(
            io,
            "  graph [layout=$(String(engine)), rankdir=\"$rankdir\", " *
            "overlap=false, splines=true, nodesep=0.55, ranksep=0.85, " *
            "pad=0.25, bgcolor=\"white\"];"
        )
        println(
            io,
            "  node [shape=ellipse, fontname=\"Helvetica\", fontsize=11, " *
            "style=\"filled\", fillcolor=\"white\", color=\"#777777\", penwidth=1.1];"
        )
        println(
            io,
            "  edge [color=\"#8a8a8a\", arrowsize=0.7, penwidth=1.0];"
        )

        for sid in Graphs.vertices(cdg)
            n_scc = size_by_scc[sid]
            n_support = support_by_scc[sid]
            iscritical = critical_by_scc[sid]
            classes = classes_by_scc[sid]

            label = "SCC $sid\\nn=$n_scc\\nsupport=$n_support"

            if iscritical
                label *= "\\ncritical class $classes"
            end

            fill = n_support > 0 ? "#eeeeee" : "white"
            colour = iscritical ? "#c62828" : "#777777"
            penwidth = iscritical ? 3.0 : 1.1

            tooltip = replace(
                labels_by_scc[sid],
                "\\" => "\\\\",
                "\"" => "\\\"",
            )

            println(
                io,
                "  s$sid [label=\"$label\", fillcolor=\"$fill\", " *
                "color=\"$colour\", penwidth=$penwidth, tooltip=\"$tooltip\"];"
            )
        end

        for e in Graphs.edges(cdg)
            println(io, "  s$(Graphs.src(e)) -> s$(Graphs.dst(e));")
        end

        println(io, "}")
    end

    engine_cmd = String(engine)
    Sys.which(engine_cmd) === nothing && throw(
        ArgumentError("Graphviz executable '$engine_cmd' was not found on PATH.")
    )
    run(Cmd([engine_cmd, "-Tsvg", dotfile, "-o", svgfile]))

    println("λ = ", result.lambda)
    println("Projective multiplicity = ", result.multiplicity)
    println("Number of SCCs = ", length(result.sccs))
    println("Critical SCCs = ", result.critical_sccs)
    println("SVG = ", svgfile)

    if !keep_dot
        rm(dotfile; force=true)
    end

    return (
        svgfile = svgfile,
        dotfile = keep_dot ? dotfile : nothing,
        analysis = result,
        summary = result.summary,
    )
end
