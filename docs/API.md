# Public API overview

## Matrix primitives

- `I`, `nB`, `nT`, `top`
- `mmp_mpower`, `rclosure`
- `tclosure`, `trclosure`

## UFNF decomposition

- `ufnf₀`, `ufnf₁`, `ufnf₂`, `ufnf₃`, `ufnf`
- `karp_max_cycle_mean_mp`

## Eigenanalysis

- `mpeigen`: right max-plus eigendirections
- `lmpeigen`: left max-plus eigendirections
- `projective_eigenvector_groups`
- `finite_right_spectrum_summary`
- `finite_right_eigenvector_support`

## Generic data/network helpers

- `load_adjacency_csv`
- `load_edgelist_csv`
- `to_maxplus`
- `network_summary`
- `analyse_network`
- `check_eigenpairs`

## Critical structure

- `eigen_direction`
- `saturation_edges`
- `critical_cycle_mean`
- `critical_structure`
- `critical_structure_summary`
- `visualize_critical_structure`

## SCC structure

- `analyse_scc_structure`
- `visualize_scc_structure`

Graphviz visualisation functions require a Graphviz executable (`dot`, `neato`, `fdp` or `sfdp`) available on the system PATH.
