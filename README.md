# Maxplus_eigenvectors

Reusable Julia tools for max-plus eigenanalysis of weighted directed networks, including reducible matrices, together with a reproducibility workflow for the *C. elegans* analyses presented in the accompanying thesis.

The mathematical package computes max-plus eigenvalues and complete right/left eigendirections, strongly and weakly connected decompositions, transitive closures, critical classes and representative critical cycles. Generic loading helpers accept square adjacency matrices and edge lists. The *C. elegans* data handling is kept in one isolated helper file so that the mathematical core remains reusable with unrelated datasets.

## Requirements

- Julia 1.12.x
- MaxPlus.jl 0.5.x
- Graphviz is optional and is required only for SVG graph visualisations (`dot`, `neato`, `fdp` or `sfdp`).

## Installation

After cloning or downloading the repository, instantiate its Julia environment:

```bash
cd maxplus_eigenvectors_28sept
julia --project=. scripts/setup.jl
```

Run the test suite:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Check the active Julia and MaxPlus versions with:

```bash
julia --project=. scripts/verify_environment.jl
```

## Matrix convention

All network matrices use

```text
A[target, source]
```

so an edge `source -> target` with weight `w` is stored as `A[target, source] = w`.

For dense numeric adjacency matrices, `0` is interpreted as an absent edge by default and is converted to max-plus bottom (`-Inf`). If zero is a genuine edge weight, encode absent edges as `-Inf` and call `to_maxplus(A; zero_is_absent=false)`. For the thesis datasets, zero denotes absence and the reproducibility notebook therefore converts them with `zero_is_absent=true` explicitly.

## Notebooks

The repository deliberately contains two different workflows.

### `notebooks/maxplus_eigenstructure_playground.ipynb`

This is the **dataset-agnostic playground**. It uses the small synthetic adjacency matrix committed under `data/examples/` and demonstrates the normal workflow for a new dataset:

```julia
input_file = joinpath(repo_root, "data", "examples", "weighted_network.csv")
data = load_adjacency_csv(input_file)
A = data.matrix
labels = data.labels

analysis = analyse_network(
    A;
    labels=labels,
    name="synthetic example",
)
```

Replace the synthetic input with any compatible matrix or edge list to use the package on another network.

### `notebooks/maxplus_eigenstructure_celegans.ipynb`

This is the **thesis reproducibility notebook**. It reconstructs the five network layers analysed in the thesis:

- chemical synapses;
- electrical gap junctions;
- monoaminergic signalling;
- neuropeptidergic signalling;
- functional signal propagation.

It explicitly separates the original numeric adjacency matrix `C` from its max-plus representation `A = to_maxplus(C)`, computes and validates the complete eigensystems, reports the SCC organisation and critical classes, reproduces the gap-junction sensitivity control, exports compact result tables, and finishes with automated checks for the main results reported in the thesis.

### Data required for the *C. elegans* notebook

To reproduce the thesis results, the notebook expects the following files:

```text
data/
└── connectomes/
    ├── SI 5 Connectome adjacency matrices, corrected July 2020 (1).xlsx
    ├── 08062023_monoamine_connectome.csv
    ├── 01022024_neuropeptide_connectome_short_range_model.csv
    └── functional/
        └── functional_connectivity_all.csv
```

These files correspond to the specific connectome releases used in the thesis:

- **Chemical synapses and electrical gap junctions.** The anatomical connectomes are taken from Cook et al. (2019), *Whole-animal connectomes of both Caenorhabditis elegans sexes*, Nature 571, 63–71, DOI: `10.1038/s41586-019-1352-7`. The notebook uses the corrected July 2020 release of Supplementary Information 5 distributed through WormWiring: `SI 5 Connectome adjacency matrices, corrected July 2020.xlsx`.

- **Monoaminergic signalling.** The monoaminergic network is based on Bentley et al. (2016), *The Multilayer Connectome of Caenorhabditis elegans*, PLOS Computational Biology 12(12): e1005283, DOI: `10.1371/journal.pcbi.1005283`. The file `08062023_monoamine_connectome.csv` is the harmonised version used in the thesis analysis.

- **Neuropeptidergic signalling.** The neuropeptidergic network is taken from Ripoll-Sánchez et al. (2023), *The neuropeptidergic connectome of C. elegans*, Neuron 111, 3570–3589.e5. The thesis uses the short-range model distributed as `01022024_neuropeptide_connectome_short_range_model.csv`, in which connections are restricted according to anatomical proximity of neuronal processes.

- **Functional signal propagation.** The functional network is based on Randi et al. (2023), *Neural signal propagation atlas of Caenorhabditis elegans*, Nature 623, 406–414, DOI: `10.1038/s41586-023-06683-4`. The file `functional_connectivity_all.csv` contains the functional signal-propagation matrix used in the thesis analysis, derived from the authors' WormNeuroAtlas data.

`src/celegans_helpers.jl` contains the dataset-specific loading and harmonisation functions used by this notebook. In particular, it aligns the structural layers to a common 302-neuron ordering and converts their source/target orientation to the package convention `A[target, source]`.

The versions listed above are the datasets used to obtain the results reported in the thesis and should therefore be used when **reproducing those results exactly**. They are not, however, a requirement of the max-plus methodology itself. Future releases or updated reconstructions of the *C. elegans* synaptic, electrical, neuromodulatory, or functional connectomes can be analysed using the same workflow, provided that the corresponding input loaders and neuron labelling are adapted when necessary. Results obtained with newer connectome versions should be regarded as new analyses rather than exact reproductions of the thesis results.


## Quick start: adjacency matrix

The generic example uses `data/examples/weighted_network.csv`, whose first column and header contain matching node labels.

```julia
using Maxplus_eigenvectors

example = load_adjacency_csv("data/examples/weighted_network.csv")
C = example.matrix
labels = example.labels

# Optional but useful when the conversion should be visible explicitly.
A = to_maxplus(C)

analysis = analyse_network(
    A;
    labels=labels,
    name="example network",
)

finite_right_spectrum_summary(analysis.right)
finite_right_eigenvector_support(analysis.right, labels)
```

## Quick start: edge list

A weighted edge list needs `source`, `target` and `weight` columns by default:

```julia
using Maxplus_eigenvectors

example = load_edgelist_csv("data/examples/weighted_edges.csv")
analysis = analyse_network(
    example.matrix;
    labels=example.labels,
    name="example network",
)
```

Column names can be changed with `source_col`, `target_col` and `weight_col`. Set `weight_col=nothing` for an unweighted graph or `directed=false` for an undirected graph.

## Direct matrix use

Data loading is optional. Any square dense or sparse Julia matrix can be analysed directly:

```julia
using Maxplus_eigenvectors
using SparseArrays

C = sparse([
    0.0  1.0  0.0
    2.0  0.0  0.0
    0.5  0.0  3.0
])

A = to_maxplus(C; zero_is_absent=true)
result = analyse_network(A; verbose=false)
```

If a matrix is already represented with `MaxPlus.MP` entries, `analyse_network` uses it directly.

## Reducible matrices and eigendirections

The implementation uses a hierarchy of UFNF decompositions:

- `UFNF0`: irreducible strongly connected block;
- `UFNF1`: condensation-DAG decomposition;
- `UFNF2`: weakly connected components;
- `UFNF3`: generic square matrices, including empty rows/columns.

The public entry points are:

```julia
spR = mpeigen(A_mp)   # right eigendirections
spL = lmpeigen(A_mp)  # left eigendirections
```

Projectively equivalent columns can be collapsed with `projective_eigenvector_groups`, and `finite_right_spectrum_summary` reports the finite projective directions.

## Critical and SCC structure

For one selected projective eigendirection:

```julia
summary = finite_right_spectrum_summary(analysis.right)
block = summary.Block[1]
direction = summary.Direction[1]

cs = critical_structure(
    analysis.matrix,
    analysis.right;
    block=block,
    direction=direction,
)

critical_structure_summary(cs, analysis.labels)

scc = analyse_scc_structure(
    analysis.matrix,
    analysis.right,
    analysis.labels;
    block=block,
    direction=direction,
)
```

If Graphviz is installed, SVG files can be generated with `visualize_critical_structure` and `visualize_scc_structure`. Generated files belong in `figures/`, which is ignored by Git apart from its placeholder.

## `src/` versus `scripts/`

The distinction is intentional:

- **`src/` contains reusable implementation code.** These files define functions and algorithms that are imported or included by analyses. The mathematical package lives here (`eigenvectors.jl`, `closures.jl`, `critical_structure.jl`, etc.). `celegans_helpers.jl` is also here because it defines reusable helper functions for the reproducibility notebook, but it is kept separate from the main `Maxplus_eigenvectors` module so the core package is not conceptually tied to *C. elegans*.
- **`scripts/` contains executable workflows.** They are meant to be run from the command line and perform a task using the code in `src/`: setting up the environment, checking versions, or running the small example. They are not part of the package API.

A useful rule is: if another Julia file should `call` a function repeatedly, that function belongs in `src/`; if the file is something a user runs from start to finish to perform a task, it belongs in `scripts/`.

## Repository structure

```text
.
├── data/
│   ├── README.md
│   └── examples/                         # small synthetic inputs only
├── docs/                                 # API and input-format notes
├── figures/                              # generated figures (ignored)
├── notebooks/
│   ├── maxplus_eigenstructure_playground.ipynb
│   └── maxplus_eigenstructure_celegans.ipynb
├── results/                              # generated tables/results (ignored)
├── scripts/                              # executable setup/example workflows
├── src/
│   ├── celegans_helpers.jl               # thesis-data loading/reproducibility helpers
│   └── ...                               # reusable max-plus implementation
└── test/                                 # automated tests
```

## Reproducibility and dependency files

`Project.toml` declares the packages and compatibility ranges required by the repository. `Manifest.toml` is intentionally ignored: cloning the repository and running `Pkg.instantiate()` resolves compatible versions for the declared environment instead of publishing a machine-specific lockfile.

`MIGRATION_NOTES.md` records the technical changes made when the implementation was migrated to MaxPlus.jl 0.5.x. It is developer documentation rather than a file needed to run an analysis.

`CITATION.cff` is machine-readable citation metadata. GitHub recognises this standard and can display a **Cite this repository** button with the software title, author, version and licence. It does not affect execution of the code.

`.github/workflows/ci.yml` defines the GitHub Actions continuous-integration job. On pushes and pull requests, GitHub creates a clean Julia 1.12 environment, instantiates the project and runs the tests. This is useful for detecting breakage that may not appear on the developer's machine.

`.gitattributes` and `.gitignore` control repository behaviour: `.gitattributes` provides Git-level file handling rules, while `.gitignore` prevents local environments, private/large datasets and generated results from being committed accidentally.

## Citation

GitHub can generate citation metadata directly from [`CITATION.cff`](CITATION.cff).

## License

MIT. See [`LICENSE`](LICENSE).
