# Script SqueezeMeta

Pipeline R per produrre grafici funzionali, tassonomici e KEGG da un progetto
SqueezeMeta.

## Esecuzione

Il punto di ingresso mantenuto è `sqm_plots_lean.R`.

```powershell
Rscript sqm_plots_lean.R `
  --project_dir in/Au_sip `
  --output_dir out/analisi `
  --mode normal
```

Modalità disponibili:

- `huge`: tutti gli output, inclusi PIE e pathway Top 20;
- `normal`: FUNZ, ENZIMI, FLOW, TAXON e PATHVIEW sui pathway definiti;
- `funz`, `flow`, `taxon`, `pie`, `pathview`: esecuzione mirata, anche in una
  lista separata da virgole.

Usare `Rscript sqm_plots_lean.R --help` per le opzioni essenziali. La
descrizione completa di default, calcoli e struttura degli output è in
[`GUIDA_LOGICA_SQM_PLOTS.md`](GUIDA_LOGICA_SQM_PLOTS.md); le regole analitiche
sono in [`REGOLE_SCRIPT_R_SQUEEZEMETA.md`](REGOLE_SCRIPT_R_SQUEEZEMETA.md).

## Pianificazione rapida

Per validare il progetto e scrivere il ranking dei pathway senza generare
grafici:

```powershell
Rscript sqm_plots_lean.R `
  --project_dir in/Au_sip `
  --output_dir out/plan `
  --mode huge `
  --plan_only
```

## Test della versione lean

```powershell
Get-ChildItem tests/test_sqm_plots_lean*.R | ForEach-Object { Rscript $_.FullName }
```

I test lean sono indipendenti dai test storici del vecchio `sqm_plots.R`.

## Struttura

- `sqm_plots_lean.R` — pipeline corrente;
- `in/Au_sip/` — dataset SqueezeMeta di prova, da preservare;
- `out/` — risultati e cache KEGG;
- `tests/test_sqm_plots_lean*.R` — test mirati della pipeline corrente;
- `archive/` — script e report storici.
