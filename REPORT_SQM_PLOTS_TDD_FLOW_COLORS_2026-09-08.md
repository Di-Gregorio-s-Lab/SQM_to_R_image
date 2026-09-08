# Report TDD — ripristino FLOW semplice

Data: 2026-09-08  
Branch: `fix/p3-correctness`

## Esito

La soluzione a palette unificata della Tranche 5 è stata sostituita perché alterava la geometria percepita del FLOW reale CS8 e rendeva l'HTML difficile da leggere.

Il renderer produttivo ora usa una rappresentazione semplice e stabile:

- ribbon colorati soltanto in base al taxon di origine;
- colonne tassonomica e KO neutre (`grey95`);
- nessuna legenda nel PNG;
- geometria di ribbon e colonne calcolata dallo stesso `flow_tbl` wide;
- Sankey HTML con layout Plotly automatico `snap`, nodi neutri e link colorati per taxon;
- etichette semplici, con TPM, KO, EC e percentuali disponibili nell'hover;
- HTML non self-contained con directory di supporto `<nome>_files`.

I dati scientifici, i TSV, il valore di K01563, i manifest e i log non sono stati modificati.

## Causa della regressione

I ribbon erano calcolati direttamente dal `flow_tbl` wide, mentre le colonne venivano ricostruite separatamente con `ggalluvial::to_lodes_form()`. La conversione di `stratum` a carattere perdeva l'ordine dei factor e disponeva i rettangoli in un ordine diverso da quello usato per i ribbon.

Nel PNG reale `flowplot_family_CS8T2_16x9.png` lo scarto verticale massimo fra ribbon e strato tassonomico era `42.06794`. Con il renderer semplice lo scarto è `2.84e-14`, cioè solo rumore numerico.

L'HTML precedente forzava inoltre le coordinate dei nodi. Il ripristino di `arrangement = "snap"` lascia a Plotly il posizionamento coerente del Sankey. La directory `_files` mancava perché il widget era stato reso self-contained; `save_html_widget()` è tornata a usare `selfcontained = FALSE`.

## Checkpoint TDD

1. `21686a1` — RED: `test: reproduce displaced FLOW strata`
2. `c930802` — GREEN: `fix: restore simple aligned FLOW PNG`
3. `9b78709` — RED: `test: reproduce forced FLOW HTML layout`
4. `fcf6b66` — GREEN: `fix: restore automatic FLOW HTML layout`
5. `b8aee25` — RED: `test: require FLOW HTML support directory`
6. `e0d878d` — GREEN: `fix: restore FLOW HTML support directory`
7. `db4a565` — RED: `test: define simple FLOW runtime dependencies`
8. `bc79745` — GREEN: `fix: simplify FLOW renderer dependencies`
9. `2073bbd` — integrazione: `test: verify simple CS8 FLOW rendering`

I commit sono locali e non è stato eseguito alcun push.

## Garanzie verificate

| Garanzia | Evidenza | Risultato |
|---|---|---|
| Ribbon e colonna tassonomica condividono gli stessi intervalli verticali | `tests/test_t5_simple_flow_regression.R` | PASS |
| Il PNG usa colori tassonomici soltanto sui ribbon, colonne neutre e nessuna legenda | `tests/test_t5_flow_colors.R` | PASS |
| L'HTML usa il layout automatico e non forza coordinate `x/y` | `tests/test_t5_simple_html_regression.R` | PASS |
| L'HTML crea e può sovrascrivere la directory `<nome>_files` | `tests/test_t5_simple_html_regression.R` | PASS |
| FLOW non richiede più `ggnewscale`, `rmarkdown` o Pandoc | test preflight e dipendenze | PASS |
| I cinque campioni CS8 mantengono K01563 entro `1e-8` | `tests/test_t5_integration_CS8.R` | PASS |
| PNG, HTML e directory di supporto vengono realmente prodotti | `tests/test_t5_integration_CS8.R` | PASS |
| Il rendering non modifica `flow_tbl` | test sintetici e integrazione CS8 | PASS |
| Oracle SQMtools/pathview/plotTaxonomy | `tests/test_t1_sqmtools_oracles.R` | PASS |
| Parità scientifica della Tranche 3 | `tests/test_t3_integration_CS8.R --oracle=fixture` | PASS |
| Lifecycle degli output della Tranche 4 | `tests/test_t4_integration_CS8.R` | PASS |

## Verifica finale

```powershell
rtk Rscript tests/run_fast_tests.R
rtk Rscript tests/test_t1_sqmtools_oracles.R

$env:SQM_CS8_PROJECT_DIR = "C:\Users\unico\OneDrive - University of Pisa\Documenti\UNIPI\Grani\Progetti\TCE\Shotgun\minion\montescudaio\CS8_All\CS8_All"

rtk Rscript tests/test_t3_integration_CS8.R --oracle=fixture
rtk Rscript tests/test_t4_integration_CS8.R
rtk Rscript tests/test_t5_integration_CS8.R
```

Risultato finale:

- 32 test rapidi GREEN;
- oracle SQMtools 1.7.2 GREEN;
- integrazioni CS8 scientifica, lifecycle e rendering FLOW GREEN;
- nessuna variazione dei valori ufficiali K01563;
- nessun artefatto KEGG aggiunto al repository.
