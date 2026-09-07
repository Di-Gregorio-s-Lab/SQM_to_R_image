# Report TDD — parità FLOW / SQMtools-pathview

Data: 2026-09-07  
Oracle: SQMtools 1.7.2 installato localmente

## Obiettivo scientifico

Per ogni pathway, campione e KO rappresentato nel flowplot, il margine
funzionale deve coincidere con `SQM$functions$KEGG$tpm`, la stessa matrice
usata da `SQMtools::exportPathway()` e quindi da pathview. Gli ORF non
definiscono il totale funzionale: servono soltanto a ripartirlo fra i taxa.

## Checkpoint TDD

1. Baseline GREEN — `4191ca4`
   (`test: align enzyme defaults and characterize SQMtools oracles`)
2. Riproduzione RED — `a773346`
   (`test: reproduce CS8 flow pathview divergence`)
3. Correzione GREEN — `998b3b7`
   (`fix(flowplot): align pathway KO totals with SQMtools`)

I commit sono locali; non è stato eseguito alcun push.

## Evidenza RED

Il riproduttore sintetico conteneva due ORF K01563, ma soltanto uno riportava
testualmente la pathway in `KEGGPATH`:

| Campione | Oracle pathview | FLOW precedente | Differenza |
|---|---:|---:|---:|
| S1 | 30 | 10 | -20 |
| S2 | 12 | 5 | -7 |

Nel progetto CS8 erano presenti quattro ORF K01563; il subset testuale
`KEGGPATH` ne conservava uno. I valori ufficiali K01563 erano:

| Campione | TPM ufficiale |
|---|---:|
| CS8T0 | 20.461747379587 |
| CS8T2 | 64.241150618838 |
| CS8T3 | 48.990246839864 |
| CS8T4 | 111.090699101642 |
| CS8T6 | 38.836884026958 |

Il FLOW precedente produceva rispettivamente
`0, 23.082, 17.819, 40.069, 2.426`.

## Causa

Pathview riceveva il `context_sqm` completo e leggeva direttamente
`SQM$functions$KEGG$tpm`. Il FLOW riceveva invece un oggetto già filtrato da
`SQMtools::subsetFun(..., columns = "KEGGPATH")`; perdeva quindi gli ORF
assegnati al KO ma privi della stringa della pathway. Inoltre attribuiva
l'intero TPM dell'ORF a ogni KO degli ORF multi-KO, mentre SQMtools divide
equamente il contributo fra i KO prima dell'aggregazione.

## Correzione GREEN

- La membership dei KO viene estratta dai nodi `ortholog` del KGML letto da
  pathview; i nodi compound sono esclusi.
- Il KGML è scaricato in una directory temporanea e rimosso subito dopo la
  lettura. Nessun XML o PNG KEGG è versionato.
- Il FLOW usa il `context_sqm` completo e seleziona tutti gli ORF associati ai
  KO della pathway, indipendentemente dal testo `KEGGPATH`.
- Per gli ORF multi-KO, il TPM viene diviso equamente fra le annotazioni KO,
  replicando `SQMtools:::aggregate_fun()`.
- Le quote tassonomiche vengono riallineate al valore ufficiale in
  `SQM$functions$KEGG$tpm` e una postcondizione verifica la conservazione del
  margine per campione e KO con tolleranza `1e-8`.
- Gli ORF vengono filtrati per KO prima dell'espansione, evitando di espandere
  l'intero progetto una volta per ogni pathway.

## Verifiche GREEN

```powershell
rtk Rscript tests/run_fast_tests.R
rtk Rscript tests/test_t1_sqmtools_oracles.R

$env:SQM_CS8_PROJECT_DIR = "C:\Users\unico\OneDrive - University of Pisa\Documenti\UNIPI\Grani\Progetti\TCE\Shotgun\minion\montescudaio\CS8_All\CS8_All"

rtk Rscript tests/test_t1_integration_CS8.R --oracle=fixture --check=oracle
rtk Rscript tests/test_t1_integration_CS8.R --oracle=fixture --check=parity
rtk Rscript tests/test_t1_integration_CS8.R --oracle=live --check=oracle
rtk Rscript tests/test_t1_integration_CS8.R --oracle=live --check=parity
```

Risultati:

- 26 file di test rapidi completati con successo;
- oracle SQMtools/pathview e `plotTaxonomy()` verdi;
- oracle fixture CS8 verde;
- parità K01563 CS8 verde;
- KGML live senza drift per i nodi K01563 delle pathway 00361 e 00625;
- controllo live aggiuntivo verde usando tutti i KO ortholog delle due
  pathway.

## Confini della modifica

Questa tranche modifica soltanto il motore dati FLOW e il relativo preflight
(`pathview` è ora una dipendenza della modalità FLOW). Non modifica la
politica di sovrascrittura, la collocazione dei manifest, i log generali, i
nomi dei file o la palette dei flowplot. Questi interventi restano assegnati
alle tranche successive.
