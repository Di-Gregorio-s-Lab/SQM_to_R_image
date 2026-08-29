# Script SqueezeMeta

Questa cartella contiene gli script R per le analisi SqueezeMeta e il dataset di prova.

## Punto di ingresso

Usare `sqm_plots.R`: è lo script attualmente mantenuto e coperto dai test in `tests/`.

```powershell
Rscript sqm_plots.R --project_dir in/Au_sip --output_dir out/<analisi> --mode <modalita>
```

Le opzioni disponibili sono descritte con `Rscript sqm_plots.R --help`. Le convenzioni analitiche sono in `REGOLE_SCRIPT_R_SQUEEZEMETA.md`.

## Struttura

- `in/Au_sip/` — dataset SqueezeMeta di prova (circa 3 GB): da preservare.
- `out/` — risultati delle analisi. `out/_smoke/` contiene output di smoke test storici.
- `tests/` — test dello script corrente.
- `archive/legacy-scripts/` — script sostituiti da `sqm_plots.R`, conservati per riferimento.
- `archive/tool-reports/` — report e cache generati dagli strumenti, non necessari per l'esecuzione.

Non sono stati cancellati file: gli elementi storici sono solo stati spostati in archivio.
