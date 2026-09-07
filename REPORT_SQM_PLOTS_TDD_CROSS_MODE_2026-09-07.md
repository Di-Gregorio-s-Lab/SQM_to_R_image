# Report TDD — armonizzazione KEGG cross-mode su SQMtools

Data: 2026-09-07  
Oracle: SQMtools 1.7.2 installato localmente

## Obiettivo scientifico

FLOW, FUNZ, PIE ed ENZIMI devono usare gli stessi totali funzionali impiegati
da `SQMtools::exportPathway()`/pathview. La membership delle pathway deriva dai
nodi `ortholog` del KGML; gli ORF servono a distribuire tassonomicamente il TPM
ufficiale dei KO, non a ridefinirne il totale. La tassonomia di pathway deve
delegare a `SQMtools::plotTaxonomy()` senza riscalare il subset.

## Checkpoint TDD

1. Baseline GREEN — 26 test rapidi e oracle SQMtools 1.7.2 verdi.
2. Riproduzione RED — `1c38fde`
   (`test: reproduce cross-mode K01563 divergence`).
3. Correzione GREEN — `39f3462`
   (`fix: unify KEGG plot data on SQMtools`).

I commit sono locali; non è stato eseguito alcun push.

## Evidenza RED

Il fixture sintetico conteneva tre ORF K01563, ma soltanto uno riportava
testualmente la pathway in `KEGGPATH` e l'EC in `KEGGFUN`. I vecchi percorsi
quantitativi producevano:

| Modalità | Campione | Oracle | Valore precedente | Differenza |
|---|---|---:|---:|---:|
| FUNZ | S1 | 100 | 10 | -90 |
| FUNZ | S2 | 50 | 5 | -45 |
| PIE | S1 | 100 | 10 | -90 |
| PIE | S2 | 50 | 5 | -45 |
| ENZIMI | S1 | 100 | 10 | -90 |
| ENZIMI | S2 | 50 | 5 | -45 |

Anche la membership tassonomica era errata: un ORF conservato su tre.

## Causa

- Pathview leggeva `SQM$functions$KEGG$tpm`.
- FUNZ e PIE partivano dal subset testuale `KEGGPATH`.
- ENZIMI cercava gli EC in `KEGGFUN` e sommava direttamente i TPM degli ORF.
- La tassonomia percentuale di pathway usava un calcolo personalizzato che
  forzava ogni pathway a 100%.

Questi percorsi non condividevano né membership né denominatore con gli oracle
SQMtools.

## Motore canonico GREEN

Il risultato condiviso espone allocazioni tassonomiche, totali ufficiali,
metadati KO/EC, ORF unici della pathway e audit di conservazione.

- Totali KO: `SQM$functions$KEGG$tpm`.
- Membership: KO deduplicati dei nodi KGML `ortholog`.
- Ripartizione tassonomica: `SQM$orfs$tpm` e `SQM$orfs$tax`.
- Metadati KO/EC: `SQM$misc$KEGG_names`.
- Tassonomia di pathway: `SQMtools::subsetORFs(..., rescale_tpm=FALSE)` e
  `SQMtools::plotTaxonomy(..., rescale=FALSE)` per `abund` e `percent`.

Per un ORF multi-KO, il TPM grezzo è diviso per il numero totale dei suoi KO
prima di filtrare quelli appartenenti alla pathway. Le quote determinano le
proporzioni tassonomiche interne al KO; tali proporzioni vengono poi applicate
al TPM ufficiale del KO. Le postcondizioni verificano, con tolleranza `1e-8`:

1. somma delle quote grezze dell'ORF uguale al TPM originale;
2. somma delle quote tassonomiche del KO uguale al TPM ufficiale.

Un KO ufficiale positivo senza ORF utilizzabili genera un errore esplicito.
Un KO ufficiale a zero non genera contributi positivi.

## Comportamento per modalità

- Pathview resta invariato e continua a usare `SQMtools::exportPathway()`.
- FLOW usa le allocazioni tassonomiche canoniche; Top-N è applicato dopo la
  normalizzazione.
- FUNZ usa direttamente i totali ufficiali dei KO unici della pathway.
- PIE usa le stesse allocazioni di FLOW; la somma delle fette coincide con il
  TPM ufficiale del KO e il totale pathway è la somma dei KO unici.
- ENZIMI costruisce una relazione deduplicata KO→EC. Più KO dello stesso EC si
  sommano; un KO con più EC contribuisce interamente a ciascun EC.
- La tassonomia di pathway include ogni ORF una volta e produce TSV con schema
  `sample, taxon, value, count`. `percent` conserva il denominatore del campione
  e non viene forzato a 100 nella pathway.

`KEGGPATH` resta usato soltanto per individuare e ordinare le pathway Top20.
Non esistono fallback quantitativi a `KEGGPATH` o `KEGGFUN`.

## Risoluzione KEGG

Gli ID della tabella curata hanno precedenza. Per nomi Top20 non curati, il
catalogo `https://rest.kegg.jp/list/pathway/ko` viene scaricato una volta per
processo, normalizzando spazi, maiuscole e il suffisso “Reference pathway”.
Mapping assenti o ambigui falliscono esplicitamente. Il catalogo live verificato
il 2026-09-07 conteneva 503 pathway.

KGML, PNG e cataloghi scaricati restano in directory temporanee e non sono
versionati.

## Verifiche GREEN

```powershell
rtk Rscript tests/run_fast_tests.R
rtk Rscript tests/test_t1_sqmtools_oracles.R

$env:SQM_CS8_PROJECT_DIR = "C:\Users\unico\OneDrive - University of Pisa\Documenti\UNIPI\Grani\Progetti\TCE\Shotgun\minion\montescudaio\CS8_All\CS8_All"

rtk Rscript tests/test_t3_integration_CS8.R --oracle=fixture
rtk Rscript tests/test_t3_integration_CS8.R --oracle=live
```

Risultati:

- 28 file di test rapidi completati con successo;
- oracle pathview e `plotTaxonomy()` SQMtools 1.7.2 verdi;
- fixture CS8 verde per le pathway 00361 e 00625;
- KGML live verde, senza drift dei nodi K01563;
- tutti e quattro gli ORF K01563 recuperati una sola volta;
- `Pathview = FLOW = FUNZ = somma PIE = ENZIMI 3.8.1.5` nei cinque campioni.

| Campione | TPM K01563 / EC 3.8.1.5 |
|---|---:|
| CS8T0 | 20.461747379587 |
| CS8T2 | 64.241150618838 |
| CS8T3 | 48.990246839864 |
| CS8T4 | 111.090699101642 |
| CS8T6 | 38.836884026958 |

## Confini

La tranche non modifica politica di sovrascrittura, collocazione dei manifest,
log generali, nomi dei file o palette. Nessun artefatto KEGG è stato aggiunto
al repository. Questi interventi restano nelle tranche successive.
