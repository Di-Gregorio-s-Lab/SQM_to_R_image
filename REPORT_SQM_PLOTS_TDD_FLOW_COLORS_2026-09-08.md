# Report TDD — mapping colori FLOW

Data: 2026-09-08  
Branch: `fix/p3-correctness`

## Obiettivo

Correggere la colonna tassonomica bianca nei FLOW PNG e usare un solo mapping nominato `categoria -> colore` in ogni grafico. Flussi, colonne, due legende PNG, nodi e link HTML devono condividere tale assegnazione senza modificare i dati scientifici.

## Checkpoint TDD

1. Baseline GREEN: 29 test rapidi.
2. RED colonna tassonomica: `da1ab8d` — `test: reproduce blank FLOW taxonomy strata`.
3. GREEN colonna tassonomica: `5473431` — `fix: color FLOW taxonomy strata`.
4. RED mapping condiviso: `b39b3b7` — `test: define shared FLOW category colors`.
5. GREEN mapping condiviso: `ba8ed07` — `fix: unify FLOW category color mapping`.
6. Integrazione CS8: `2561ee9` — `test: verify CS8 FLOW color parity`.

I commit sono locali e non è stato eseguito alcun push.

## Evidenza RED

Il primo test costruiva il grafico con `ggplot2::ggplot_build()` e osservava:

```text
Alpha: flow=#E32636, stratum=WHITE, legend=#E32636
Beta:  flow=#5D8AA8, stratum=WHITE, legend=#5D8AA8
Other: flow=GREY70, stratum=WHITE, legend=GREY70
```

La conversione di `stratum` a testo, da sola, non era sufficiente: `StatStratum` ricalcolava la variabile dopo il mapping. La correzione usa quindi `after_stat(stratum)` per il riempimento delle colonne.

Il secondo RED mostrava due palette indipendenti:

```text
taxonomy=#5D8AA8,#E32636,GREY70
function=#5D8AA8,#E32636,GREY70
expected function=#EFDECD,#FFBF00,GREY70
```

Tassoni e KO distinti ripartivano dagli stessi primi colori.

## Implementazione GREEN

- `build_flow_color_map()` costruisce un unico vettore nominato usando prima i livelli tassonomici e poi quelli KO.
- Categorie mancanti, vuote e duplicate vengono eliminate; `Other` compare una sola volta ed è sempre `grey70`.
- La palette usa `unique(colors_hex)`. Se le categorie superano i colori disponibili, viene generata una palette HCL deterministica; colori mancanti o duplicati causano un errore esplicito.
- Le due scale PNG restano separate per titoli, break e label, ma ricevono lo stesso mapping completo.
- Gli alluvia mantengono il colore del taxon sorgente; entrambe le colonne usano lo strato calcolato da `StatStratum`.
- Il Sankey usa lo stesso mapping per i nodi e deriva i link dal colore del taxon sorgente applicando soltanto la trasparenza.
- Nessuna modifica è stata apportata a CLI, TSV, motore KEGG, manifest, log o nomi degli output.

## Verifiche

| Garanzia | Evidenza | Risultato |
|---|---|---|
| Flusso, colonna tassonomica e legenda coincidono | `tests/test_t5_flow_colors.R` | PASS |
| Tassoni e KO distinti non riutilizzano colori | `tests/test_t5_flow_colors.R` | PASS |
| `Other` rimane `grey70` nelle due sezioni | `tests/test_t5_flow_colors.R` | PASS |
| Il mapping non dipende dall'ordine delle righe | `tests/test_t5_flow_colors.R` | PASS |
| PNG e HTML condividono i colori | `tests/test_t5_flow_colors.R` | PASS |
| L'overflow usa colori HCL non riciclati | `tests/test_t5_flow_colors.R` | PASS |
| I cinque campioni CS8 renderizzano PNG e HTML non vuoti | `tests/test_t5_integration_CS8.R` | PASS |
| Il TPM K01563 CS8 resta invariato entro `1e-8` | `tests/test_t5_integration_CS8.R` | PASS |
| Oracle SQMtools/pathview/plotTaxonomy | `tests/test_t1_sqmtools_oracles.R` | PASS |
| Parità scientifica Tranche 3 fixture | `tests/test_t3_integration_CS8.R --oracle=fixture` | PASS |
| Lifecycle output Tranche 4 | `tests/test_t4_integration_CS8.R` | PASS |

Comandi eseguiti:

```powershell
rtk Rscript tests/run_fast_tests.R
rtk Rscript tests/test_t1_sqmtools_oracles.R

$env:SQM_CS8_PROJECT_DIR = "C:\Users\unico\OneDrive - University of Pisa\Documenti\UNIPI\Grani\Progetti\TCE\Shotgun\minion\montescudaio\CS8_All\CS8_All"
rtk Rscript tests/test_t3_integration_CS8.R --oracle=fixture
rtk Rscript tests/test_t4_integration_CS8.R
rtk Rscript tests/test_t5_integration_CS8.R
```

Risultato finale:

- 30 test rapidi GREEN;
- oracle SQMtools 1.7.2 GREEN;
- integrazioni CS8 scientifica, lifecycle e colori GREEN;
- nessuna variazione dei valori ufficiali K01563;
- nessun artefatto KEGG aggiunto al repository.
