# Data formats

The package accepts arbitrary square numeric matrices directly. The CSV helpers below are convenience functions, not requirements.

## 1. Labeled adjacency matrix

Use `load_adjacency_csv(path)` for a CSV in which the first column contains row labels and the remaining column headers contain the same node labels.

Example:

```csv
node,A,B,C
A,0,1,0
B,2,0,0
C,0.5,0,3
```

The package convention is `A[target, source]`. In this example, the value `2` in row `B`, column `A` represents the edge `A -> B` with weight `2`.

The loader checks that row and column labels match and reorders matrix columns to the row-label order if necessary.

For dense matrices, numeric zero means "no edge" by default when the matrix is converted to max-plus form. If zero is a valid edge weight, use `-Inf` for absent edges and pass `zero_is_absent=false` to `to_maxplus` or `analyse_network`. Edge-list inputs are stored sparsely, so explicitly stored zero-weight edges can be preserved.

## 2. Weighted edge list

Use `load_edgelist_csv(path)` for:

```csv
source,target,weight
A,B,2.0
B,A,1.0
B,C,0.5
```

The loader creates an adjacency matrix using `A[target, source]`. Duplicate edges are summed.

Useful options:

```julia
load_edgelist_csv(
    path;
    source_col=:from,
    target_col=:to,
    weight_col=:strength,
    directed=true,
)
```

For unweighted data:

```julia
load_edgelist_csv(path; weight_col=nothing)
```

For an undirected edge list:

```julia
load_edgelist_csv(path; directed=false)
```

## 3. Other file types

For XLSX, JLD2, Arrow, HDF5, database queries or domain-specific formats, load the data with the appropriate Julia package and provide the resulting square matrix directly:

```julia
result = analyse_network(A; labels=labels)
```

Keeping generic input handling separate from study-specific preprocessing makes the eigenanalysis reusable across domains. The optional `src/celegans_helpers.jl` file contains the thesis-specific harmonisation used only by the C. elegans reproducibility notebook.
