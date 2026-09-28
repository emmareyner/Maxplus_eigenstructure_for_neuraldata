# Data directory

Only synthetic example files are version-controlled in `data/examples/`.

Local datasets can be placed elsewhere under `data/`; `.gitignore` excludes them automatically. This prevents large, licensed, private or study-specific datasets from being committed accidentally.

For the thesis reproducibility notebook, use:

```text
data/connectomes/
├── SI 5 Connectome adjacency matrices, corrected July 2020 (1).xlsx
├── 08062023_monoamine_connectome.csv
├── 01022024_neuropeptide_connectome_short_range_model.csv
└── functional/
    └── functional_connectivity_all.csv
```

See `docs/DATA_FORMATS.md` for the generic CSV layouts used by the playground.
