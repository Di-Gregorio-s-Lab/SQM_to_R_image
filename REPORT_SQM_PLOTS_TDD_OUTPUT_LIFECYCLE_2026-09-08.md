# Report TDD — lifecycle degli output e log di esecuzione

Data: 2026-09-08  
Branch: `fix/p3-correctness`

## Obiettivo

Una nuova esecuzione nella stessa `output_dir` deve sovrascrivere soltanto gli
artefatti rigenerati. Gli altri file devono rimanere invariati. I manifest sono
locali alle sezioni FLOW, FUNZ e PIE; Pathview e `plotTaxonomy()` non producono
manifest. Ogni invocazione scrive `output_dir/<run_id>.log`.

## Checkpoint TDD

1. Baseline GREEN: 28 test rapidi.
2. RED: `c9de983` — `test: define output overwrite lifecycle`.
3. GREEN: `58aaf97` — `fix: overwrite managed outputs and localize manifests`.

I commit sono locali e non è stato eseguito alcun push.

## Evidenza RED

Comando:

```powershell
rtk Rscript tests/test_t4_output_lifecycle.R
```

Il test veniva eseguito e falliva sull'invariante iniziale:

```text
run_artifact_path() still adds a run-specific suffix
```

La causa era il suffisso `__<run_id>` applicato da tutti i writer. I manifest
erano inoltre aggregati per sezione nella root, indicizzati da `manifest_all` e
affiancati da manifest di corsa e fallimento.

## Implementazione GREEN

- TSV, PNG, HTML e output Pathview usano percorsi deterministici senza
  `run_id` e sovrascrivono il target omonimo.
- Non viene eseguita alcuna pulizia ricorsiva: file non selezionati e file
  estranei restano invariati.
- Ogni contesto mantiene un registry in memoria. FLOW, FUNZ e PIE vengono
  serializzati rispettivamente in `manifest_flow.tsv`, `manifest_funz.tsv` e
  `manifest_pie.tsv` nella propria directory fisica.
- I manifest sono sostituiti tramite staging temporaneo e ripristino del file
  precedente se la sostituzione non riesce.
- In caso di errore vengono serializzate le righe già registrate; non viene
  effettuato rollback degli artefatti completati.
- Pathview, tassonomia globale e tassonomia di pathway non scrivono manifest.
- I vecchi artefatti e manifest rimangono intatti e vengono segnalati nel log.
- Successo e fallimento producono un log distinto `output_dir/<run_id>.log` con
  configurazione, versioni, avanzamento, warning, skip, riepilogo e stato.

## Specifica verificata

| Garanzia | Evidenza | Tipo | Risultato |
|---|---|---|---|
| Una rerun usa lo stesso percorso e sovrascrive il contenuto | `tests/test_t4_output_lifecycle.R` | unità | PASS |
| File e manifest di sezioni non selezionate restano byte-identici | `tests/test_t4_output_lifecycle.R` | unità | PASS |
| Solo FLOW, FUNZ e PIE hanno manifest locali | `tests/test_t4_output_lifecycle.R` | unità | PASS |
| Successo e fallimento producono log distinti | test lifecycle e preflight CLI | integrazione CLI | PASS |
| I manifest legacy restano invariati e sono segnalati | `tests/test_t4_output_lifecycle.R` | unità | PASS |
| Due run reali CS8 sovrascrivono FLOW e conservano una sentinella | `tests/test_t4_integration_CS8.R` | integrazione | PASS |
| Run reali Pathview e taxon non producono manifest | `tests/test_t4_integration_CS8.R` | integrazione | PASS |

## Comandi GREEN eseguiti

```powershell
rtk Rscript tests/run_fast_tests.R
rtk Rscript tests/test_t1_sqmtools_oracles.R

$env:SQM_CS8_PROJECT_DIR = "C:\Users\unico\OneDrive - University of Pisa\Documenti\UNIPI\Grani\Progetti\TCE\Shotgun\minion\montescudaio\CS8_All\CS8_All"

rtk Rscript tests/test_t3_integration_CS8.R --oracle=fixture
rtk Rscript tests/test_t3_integration_CS8.R --oracle=live
rtk Rscript tests/test_t4_integration_CS8.R
```

Risultati:

- 29 test rapidi GREEN;
- oracle SQMtools 1.7.2 GREEN;
- integrazione scientifica CS8 fixture e live GREEN;
- integrazione lifecycle CS8 GREEN con quattro run e quattro log;
- nessun artefatto KEGG aggiunto al repository.

Non esiste nel repository un misuratore di coverage R configurato; la verifica
usa la suite rapida completa, i test CLI dei fallimenti e le integrazioni reali
CS8. La correzione dei colori FLOW resta fuori da questa tranche.
