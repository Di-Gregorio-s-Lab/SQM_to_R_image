# Lean SqueezeMeta plotting pipeline

`sqm_plots_lean.R` generates functional, taxonomic, flow, pie, enzyme, KEGG taxa stacked bar, and KEGG Pathview outputs from a SqueezeMeta project loaded through SQMtools.

The script can run a complete analysis or selected output modes. It validates the command line and the SQM object before rendering, records task errors in a TSV file, and leaves the input project unchanged.

## Requirements

- R with `Rscript` available on `PATH`
- A completed SqueezeMeta project that `SQMtools::loadSQM()` can read
- Internet access on the first run, or an existing KEGG pathway catalog cache
- The R packages required by the selected mode

Install the CRAN packages with:

```r
install.packages(c("SQMtools", "ggplot2", "ggalluvial", "plotly", "htmlwidgets"))
```

Install Pathview from Bioconductor with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("pathview")
```

`SQMtools` and `ggplot2` are always checked. FLOW also needs `ggalluvial`; HTML FLOW output needs `plotly` and `htmlwidgets`. Every non-`--plan_only` run currently checks for `pathview`, even when the selected mode does not render pathway maps.

## Basic use

Run the normal profile:

```sh
Rscript sqm_plots_lean.R \
  --project_dir /path/to/squeezemeta-project \
  --output_dir out/analysis \
  --mode normal
```

PowerShell uses backticks for line continuation:

```powershell
Rscript sqm_plots_lean.R `
  --project_dir C:/path/to/squeezemeta-project `
  --output_dir out/analysis `
  --mode normal
```

Use `Rscript sqm_plots_lean.R --help` for the short CLI reference. [GUIDE.md](GUIDE.md) documents every option, calculation, and output path.

## Modes

| Mode | Output |
|---|---|
| `huge` | FUNZ, enzyme, FLOW, TAXON, PIE, KEGG taxa, and PATHVIEW outputs |
| `normal` | FUNZ, enzyme, FLOW, TAXON, KEGG taxa, and PATHVIEW for the defined pathways |
| `funz` | KO bar plots and enzyme plots |
| `flow` | Taxon to KO alluvial PNG and Sankey HTML outputs |
| `taxon` | Global and pathway-specific taxonomy plots |
| `pie` | Taxonomic composition pies for each sample and KO |
| `kegg_taxa` | Stacked TPM bars by pathway, sample, and taxonomy rank, with a per-image TPM legend |
| `pathview` | KEGG Pathview exports |

Elementary modes can be combined as a comma-separated list, for example `--mode funz,flow`. The `huge` and `normal` profiles must be used alone.

The values `funz`, `definiti`, `enzimi`, `insieme`, and `separato` are retained in options or output paths for compatibility.

## Plan a run

`--plan_only` validates the project, ranks pathways for each analysis context, writes `top20.tsv`, and stops before rendering:

```sh
Rscript sqm_plots_lean.R \
  --project_dir /path/to/squeezemeta-project \
  --output_dir out/plan \
  --mode huge \
  --plan_only
```

## Outputs

The selected modes write under the requested output directory:

```text
funz/                 KO and enzyme plots, data, and manifest
flowplot/             alluvial and Sankey outputs, data, and manifest
taxonomy_global/      global SQMtools taxonomy plots
taxonomy_by_pathway/  pathway-specific SQMtools taxonomy plots
pie/                  taxonomic pie plots, data, and manifest
kegg_taxa/            pathway-by-sample stacked taxon TPM plots and TSV data
pathview/             SQMtools Pathview exports
top20.tsv              contextual pathway ranking
errors.tsv             task errors for the completed run
<run_id>.log           arguments, warnings, errors, and final status
_cache/kegg/           reusable KEGG pathway catalog cache
```

Existing output directories are not cleared. Files with the same deterministic name may be overwritten, while unrelated files remain in place.

KEGG taxa plots select up to 10 classified taxa per pathway and sample bar, each contributing at least 1% of that bar's TPM by default. Taxa selected in any sample of a pathway remain visible in all its samples; remaining classified taxa are grouped as `Other`, and `Unclassified` stays separate. Set `--top_n_kegg_taxa` and `--min_kegg_taxon_percent` to change these limits independently of FLOW and PIE. Each requested dimension produces PNG pages for every taxonomy rank and pathway selection, with a color and TPM legend by pathway and sample on every page. TSV data are saved alongside the PNGs.

## Documentation

- [GUIDE.md](GUIDE.md) describes the CLI, data flow, calculations, and generated files.

## Known limitations

- The first run needs access to the KEGG REST API unless the pathway catalog cache already exists.
- `--refresh_kegg` refreshes the pathway catalog used by the current pipeline. KGML cache helpers exist in the script but do not determine pathway membership.
- TAXON and PATHVIEW do not produce section manifests.
- TAXON uses the native `SQMtools::plotTaxonomy()` behavior and produces PNG files without taxonomy sidecar TSV files.
- The pipeline records independent task failures and exits with status `1` when `errors.tsv` is not empty. Partial outputs remain available.

## License

Copyright 2026 Giacomo Bernabei.

This repository is available under the [PolyForm Noncommercial License 1.0.0](LICENSE). Noncommercial use is permitted under that license. Commercial use requires a separate agreement with the copyright holder. This is a source-available license, not an OSI-approved open source license.
