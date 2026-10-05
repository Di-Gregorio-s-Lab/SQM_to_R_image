# Guide to `sqm_plots_lean.R`

## Purpose

This guide describes the behavior of `sqm_plots_lean.R`: accepted inputs, pathway selection, TPM calculations, rendering, and generated files.

The file can be run with `Rscript` or sourced to reuse its functions. `main()` runs only during direct script execution.

## Execution flow

```text
parse and validate CLI arguments
check packages required by the selected modes
load the project once with SQMtools::loadSQM()
validate ORF identifiers, samples, taxonomy ranks, and TPM matrices
create one global context or one context per requested taxon
rank pathways and write top20.tsv in each context
resolve defined and top20 pathway tasks
prepare ORF, sample, and KO data for each pathway
render the requested modes
write section manifests, errors.tsv, and the run log
```

When `--taxa` is present, the script creates independent taxon-filtered contexts and does not also render the global context.

## Command line

Arguments accept both `--name value` and `--name=value` forms. `--project_dir`, `--output_dir`, and `--mode` are required for an analysis run.

```sh
Rscript sqm_plots_lean.R \
  --project_dir /path/to/project \
  --output_dir out/analysis \
  --mode normal
```

Unknown options, missing values, invalid enumerations, and invalid positive integers stop the run.

| Option | Default | Meaning |
|---|---|---|
| `--mode` | required | `huge`, `normal`, or a comma-separated subset of `funz,flow,taxon,pie,kegg_taxa,pathview` |
| `--samples` | all samples | Samples to use, in the given order |
| `--tax_mode` | `prokfilter` | `loadSQM()` taxonomic filter: `prokfilter`, `allfilter`, or `nofilter` |
| `--taxa` | none | Taxa to process as separate contexts |
| `--pathways` | ten curated IDs | Pathways used by the `defined` selection |
| `--pathway_selection_modes` | `defined,top20` | Pathway selection strategies |
| `--pathway_top_n` | `20` | Number of ranked pathways |
| `--top_n_ko` | `20` | KOs kept as separate categories in FUNZ and FLOW |
| `--top_n_taxa` | `15` | Taxa kept as separate categories in FLOW and PIE |
| `--top_n_kegg_taxa` | `10` | Maximum classified taxa selected independently in each KEGG pathway and sample bar |
| `--min_kegg_taxon_percent` | `1` | Minimum relative contribution, in percent of the pathway and sample TPM, for a classified taxon to be selected in that bar |
| `--taxonomy_ranks` | `phylum,class,order,family,genus,species` | Taxonomy ranks to render |
| `--taxonomy_counts` | `abund,percent` | Count types passed to `plotTaxonomy()` |
| `--flowplot_formats` | `png,html` | FLOW formats |
| `--pathview_sample_modes` | `insieme,separato` | Combined and per-sample Pathview exports |
| `--enzyme_ecs` | curated list | EC codes included in enzyme plots |
| `--enzyme_plot_types` | `bar,line` | Enzyme plot types |
| `--dimensions` | `12x9,16x9,12x16` | Plot dimensions in inches |
| `--plot_dpi` | `600` | PNG resolution |
| `--workers` | up to 4 physical cores | Parallel rendering workers |
| `--plan_only` | false | Write contextual pathway rankings without rendering |
| `--refresh_kegg` | false | Refresh the KEGG pathway catalog cache |

`normal` forces pathway selection to `defined`. By default, PIE uses only `defined`, including under `huge`; passing `--pathway_selection_modes` explicitly can include `top20` for PIE. KEGG taxa bar plots use the current pathway selection: `defined` under `normal`, or separate `defined` and `top20` groups when both are selected. The two KEGG taxa selection options are independent of `--top_n_taxa` for FLOW and PIE.

The former `all` and `enzimi` modes are not accepted. Enzyme plots are part of `funz`.

## Project loading and validation

The script loads the project as follows:

```r
SQMtools::loadSQM(
  project_path = normalized_project_dir,
  tax_mode = tax_mode,
  trusted_functions_only = FALSE,
  load_sequences = FALSE
)
```

The loaded object must provide `sqm$orfs$table`, `sqm$orfs$tax`, `sqm$orfs$tpm`, and `sqm$functions$KEGG$tpm`. The three ORF tables need unique row names and matching ORF ID sets. TPM matrices must contain finite, nonnegative numeric values without missing entries in the selected samples.

If `--samples` is omitted, the script uses every column of `sqm$orfs$tpm`. Requested samples and taxonomy ranks must exist.

## Taxon contexts

Each `--taxa` value is matched exactly after trimming whitespace and ignoring letter case. A taxon must occur in exactly one taxonomy column. Missing or ambiguous names stop the run.

Matching ORFs are passed to `SQMtools::subsetORFs()` with TPM and copy number rescaling disabled. Results are written under:

```text
<output_dir>/taxon_filter/<rank>/<taxon>/
```

## Pathway selection

### Defined pathways

A requested pathway can be a five-digit KEGG ID or a full KEGG pathway name. The script resolves it against the KEGG `pathway/ko` catalog and requires one match.

The resolved name is passed to:

```r
SQMtools::subsetFun(
  SQM = sqm,
  fun = pathway_name,
  columns = "KEGGPATH",
  ignore_case = FALSE,
  fixed = TRUE,
  allow_empty = FALSE
)
```

The KO row names in the subset's official KEGG TPM matrix define the pathway KOs used by FUNZ, FLOW, PIE, and KEGG taxa bar plots. A pathway that cannot be prepared is added to `errors.tsv`; other tasks continue.

### Top pathways

The ranking is recalculated inside every analysis context with the selected samples. It reads three-level entries from `sqm$misc$KEGG_paths`, keeps the six KEGG PATHWAY roots, and excludes categories marked "not included". Each pathway receives the sum of its member KO values from `sqm$functions$KEGG$tpm`.

A KO can contribute its full TPM to more than one pathway. Ties are ordered by pathway name, root, and category. Duplicate pathway names in different hierarchies are rejected because they cannot be resolved safely.

Every context writes `top20.tsv`. Its columns are `rank`, `pathway_id`, `pathway_name`, `pathway_root`, `pathway_category`, `total_tpm`, `samples`, `filtered_taxon`, and `filtered_taxon_rank`. Rows without a resolvable KEGG ID remain in the ranking but do not become rendering tasks.

## TPM allocation

FUNZ, FLOW, PIE, and KEGG taxa bar plots share one prepared data table. The script extracts every `Kxxxxx` identifier from each ORF's `KEGG ID` value. For each sample, it divides the ORF TPM equally across all KOs annotated on that ORF, then keeps the shares for KOs in the prepared pathway subset.

```text
allocated ORF-KO TPM = ORF TPM / number of KOs annotated on the ORF
```

The script does not rescale these shares. For each positive sample and KO pair, their sum must match the pathway subset value in `sqm$functions$KEGG$tpm` within a relative tolerance of `1e-8`. A mismatch stops preparation for that pathway and records observed and expected values.

Missing taxonomy values become `Unclassified`. KO names and EC codes come from `sqm$misc$KEGG_names`; EC codes are read from the `[EC:...]` block.

## Analytical outputs

### FUNZ and enzymes

FUNZ ranks KOs by their official TPM sum over the selected samples. It keeps `top_n_ko` KOs and combines the rest as `Other`. Percentages use the sample's pathway KO total as the denominator.

```text
funz/pathway/<definiti|top20>/<pathway>/
  barplot_ko_data.tsv
  barplot_ko_<dimension>.png
```

Enzyme plots use the complete official KO TPM matrix for the context, not a single pathway subset. KO values sharing an EC code are added together for each sample.

```text
funz/enzimi/insieme/
funz/enzimi/separato/<EC>/
```

Each enzyme directory contains `enzimi_data.tsv` and the requested bar or line PNG files.

### FLOW

FLOW chooses the leading KOs and taxa from totals across the selected samples, then writes a separate table for each sample. `Unclassified` remains distinct; `Other` contains categories outside the selected top values.

```text
flow_percent  = 100 * edge TPM / sample pathway TPM
taxon_percent = 100 * taxon TPM / sample pathway TPM
ko_percent    = 100 * KO TPM / sample pathway TPM
```

```text
flowplot/<definiti|top20>/<pathway>/<rank>/
  flow_<sample>_data.tsv
  flow_<sample>_<dimension>.png
  flow_<sample>.html
```

The PNG is an alluvial taxon to KO plot. The self-contained HTML file is a Sankey diagram with KO metadata and TPM values in its hover text.

### TAXON

Global and pathway taxonomy plots are delegated to `SQMtools::plotTaxonomy()` with `N = 15`, `rescale = FALSE`, `ignore_unmapped = TRUE`, `ignore_unclassified = TRUE`, and `no_partial_classifications = FALSE`.

```text
taxonomy_global/<abund|percent>/<rank>/taxonomy_<dimension>.png
taxonomy_by_pathway/<definiti|top20>/<pathway>/<abund|percent>/<rank>/taxonomy_<dimension>.png
```

TAXON writes PNG files only. It does not create taxonomy TSV files or a section manifest.

### PIE

For each sample, pathway KO, and rank, PIE expresses allocated taxon TPM against the official KO total, keeps `top_n_taxa` taxa, and combines the rest as `Other`.

```text
pie/<definiti|top20>/<pathway>/<sample>/<KO>/<rank>/
  pie_data.tsv
  pie_<dimension>.png
```

`pct` is a proportion between 0 and 1. Combinations without positive TPM produce no files.

### KEGG taxa stacked bars

The `kegg_taxa` mode renders one stacked bar per pathway and sample, with pathway IDs along the X axis and TPM on the Y axis. Each requested taxonomy rank has its own images. The `defined` and `top20` selections remain separate; all selected samples appear side by side within each pathway.

For each pathway and sample bar, classified taxa are ranked by allocated TPM. Up to `--top_n_kegg_taxa` taxa contributing at least `--min_kegg_taxon_percent` of that bar's total TPM are selected. The union of the taxa selected in any sample of the pathway is displayed in every sample of that pathway, including samples where a selected taxon's contribution is below the threshold. Remaining classified taxa become `Other`; `Unclassified` remains separate. The stacked values sum to each bar's full pathway TPM, including zero-TPM bars.

The same taxon keeps its palette color across pathways and pages of a rank. Every PNG includes a table-like legend with the color, taxon name, and TPM by pathway and sample for the pathways shown on that page. The TSV retains the plotted values for checking the bars and legend.

```text
kegg_taxa/<definiti|top20>/<rank>/
  kegg_taxa_data.tsv
  kegg_taxa_<dimension>_page_<number>.png
```

The requested `--dimensions` are rendered separately. Pathways are split across numbered PNG pages when the bars and full legend cannot fit legibly. With many samples, the legend repeats taxon rows in blocks of sample columns. A page containing one pathway grows vertically if needed so its legend remains complete. `kegg_taxa/manifest_kegg_taxa.tsv` lists the generated TSV and PNG files.

### PATHVIEW

PATHVIEW delegates rendering to `SQMtools::exportPathway()` with `count = "tpm"`, `split_samples = FALSE`, and `log_scale = FALSE`.

```text
pathview/<definiti|top20>/<insieme|separato>/<pathway>/
```

`insieme` exports all selected samples in one call. `separato` creates one sample directory per call. PATHVIEW has no section manifest.

## Cache, errors, and logs

The KEGG pathway catalog is stored at `<output_dir>/_cache/kegg/pathway_catalog.tsv`. A valid cache is reused. An invalid cache is downloaded again, and `--refresh_kegg` forces a new catalog download. The script also contains KGML cache helpers, but the current pipeline does not call them to determine pathway membership.

FUNZ, FLOW, and PIE write `manifest_funz.tsv`, `manifest_flow.tsv`, and `manifest_pie.tsv` in their section roots. Each manifest lists only files produced by that section and uses paths relative to the section root. TAXON and PATHVIEW have no manifests.

Task failures are collected in `<output_dir>/errors.tsv` with `task_id`, `mode`, `pathway_id`, and `error`. The pipeline attempts other independent tasks and exits with status `1` if the final table contains errors.

Every CLI run also writes `<output_dir>/<run_id>.log` when the output directory can be parsed. Existing output directories are not cleared. A regenerated filename can be overwritten, while other files remain untouched.

## Current limits

- The input SqueezeMeta project is never modified.
- Pathway membership follows the `SQMtools::subsetFun()` subset, not KGML ortholog nodes.
- `top20.tsv` is written for every context, including `normal` and `--plan_only` runs.
- Rendering tasks can run in parallel, but completed results and errors are sorted to keep output deterministic.
- Every non-`--plan_only` run checks for `pathview`, regardless of the selected output mode.
- `--refresh_kegg` refreshes the catalog used by the pipeline. It does not cause a KGML download during a normal run.
