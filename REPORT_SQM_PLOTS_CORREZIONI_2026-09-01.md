# Report finale delle correzioni a `sqm_plots.R` — 2026-09-01

## Esito

Il piano approvato è stato implementato senza rifattorizzare il monolite. Le
correzioni scientifiche P0-P3 già chiuse non sono regredite nei casi verificati e i
nuovi contratti su Pathview, filtri tassonomici, pathway vuoti, ordine dei campioni e
provenienza degli artefatti sono coperti da test sintetici e reali.

Il risultato introduce un'identità distinta per ogni invocazione analitica. Le
directory restano inalterate, mentre ogni file prodotto riceve il suffisso
`__<run_id>` prima dell'estensione. I manifest non sono più cumulativi: descrivono
esclusivamente la corsa che li ha generati.

## Perimetro e decisioni applicate

### Pathview

- `SQMtools::exportPathway()` riceve esplicitamente `log_scale = FALSE` e continua a
  usare TPM lineari.
- Il TSV di configurazione registra `log_scale=FALSE`, `pseudocount=NA`,
  `color_source=pathview_native` e `input_scope=complete_all_ko_matrix`.
- Non viene passata una palette esterna: i colori restano quelli nativi di Pathview.
- Il TSV sorgente è denominato
  `pathview_input_all_ko_complete_matrix__<run_id>.tsv`. Il nome dichiara
  esplicitamente che contiene la matrice completa dei KO forniti a Pathview, non
  soltanto i nodi effettivamente disegnati sulla mappa.
- Ogni invocazione Pathview continua a lavorare in una directory temporanea isolata;
  soltanto i file prodotti in quella invocazione vengono trasferiti e manifestati.

Questa è una decisione di contratto: la matrice completa è la sorgente quantitativa
fornita al renderer, ma non pretende di ricostruire la selezione interna dei nodi o la
palette nativa applicata dal servizio Pathview/KEGG.

### Filtraggio tassonomico

- `subsetORFs()` usa `ignore_unclassified_functions = FALSE`.
- Gli ORF senza funzione classificata non vengono eliminati dal filtro tassonomico.
- Restano attive le postcondizioni sull'identità esatta degli `orf_id`.
- I taxa con una sola ORF conservano tabelle ORF, tassonomia e TPM rettangolari a una
  riga.

Il caso reale `Accipitriformes (no class in NCBI)` ora restituisce correttamente una
ORF anziché fallire dentro SQMtools.

### Pathway validi ma vuoti

- `subsetFun()` riceve `allow_empty = TRUE`.
- Il riconoscimento del subset vuoto è centralizzato.
- Una combinazione `contesto × pathway` vuota emette un warning, viene registrata con
  motivazione `empty_subset` e viene saltata senza interrompere gli altri pathway o
  la tassonomia globale.
- Non vengono creati artefatti fittizi né righe artefatto verso file inesistenti.

Nel caso reale `Archaeoglobi`, nove pathway vuoti sono stati saltati e la tassonomia
globale ha comunque prodotto i due artefatti attesi.

### Ordine dei campioni e line plot enzimatici

- Con `--samples`, l'ordine CLI viene conservato esattamente e il manifest registra
  `sample_order_basis=cli`.
- Senza `--samples`, viene mantenuto l'ordine delle colonne SQM e il manifest registra
  `sample_order_basis=sqm_column_order`.
- I campioni duplicati vengono rifiutati perché non possono rappresentare due
  posizioni di visualizzazione distinte.
- I line plot enzimatici restano disponibili. La documentazione chiarisce che la
  linea segue un ordine di visualizzazione scelto dall'utente o ereditato da SQM e
  non dimostra una sequenza temporale.

### Identità della corsa, file e manifest

All'inizio della corsa viene allocato un solo identificatore nel formato:

```text
YYYYMMDDTHHMMSS_UTCpHHMM_<4hex>
YYYYMMDDTHHMMSS_UTCmHHMM_<4hex>  # offset negativo
```

Lo stesso `run_id` viene propagato a tutti gli artefatti e a ogni riga dei manifest.
Le directory analitiche non incorporano il codice della corsa.

I manifest prodotti, quando la relativa sezione esiste, sono:

```text
funz/manifest_funz__<run_id>.tsv
flowplot/manifest_flow__<run_id>.tsv
manifest_taxon__<run_id>.tsv
pathview/manifest_pathview__<run_id>.tsv
pie/manifest_pie__<run_id>.tsv
manifest_all__<run_id>.tsv
manifest_run__<run_id>.tsv
```

`manifest_all` indicizza soltanto i manifest di sezione della corsa corrente. Nessun
manifest storico viene letto o fuso durante la scrittura. `manifest_run` registra
stato, inizio/fine, progetto e output, modalità, argomenti CLI, campioni e loro ordine,
warning, pathway saltati ed eventuale errore.

Se un errore avviene dopo la risoluzione di `output_dir`, gli artefatti parziali
restano al loro posto e vengono inventariati in
`manifest_failed_artifacts__<run_id>.tsv`; `manifest_run__<run_id>.tsv` riceve
`status=failed` e il processo termina comunque con codice non zero. Gli errori CLI che
non permettono di determinare l'output terminano senza creare un manifest.

## Strategia TDD

Il lavoro è stato separato in commit distinti:

- `376fb32 fix: remediate sqm plots audit findings`: correzioni precedenti già
  presenti nel worktree e committate prima del nuovo intervento;
- `68a5a04 test: add sqm run contract regressions`: test RED del nuovo contratto,
  committati prima dell'implementazione;
- implementazione corrente: modifica minima del monolite, aggiornamento della guida,
  delle regole e degli harness di verifica.

La fase RED ha confermato fallimenti per la policy tassonomica, per il valore
Pathview `log_scale` e per l'assenza del contratto `run_id`/manifest. Dopo
l'implementazione gli stessi test sono verdi.

I test aggiunti o estesi verificano:

- formato, unicità, allocazione e propagazione del `run_id`;
- suffisso sui file senza modifica delle directory;
- separazione di due manifest prodotti nello stesso output;
- manifest di fallimento e conservazione degli artefatti parziali;
- errore di parsing successivo a un `output_dir` già determinabile;
- round-trip del writer TSV di emergenza con tab, newline e virgolette;
- `allow_empty=TRUE`, warning e skip selettivo del pathway;
- taxon sintetico e reale con una sola ORF;
- coerenza tra argomenti e configurazione Pathview lineare;
- colori nativi e nome/semantica della matrice KO completa;
- ordine CLI/SQM e relativa provenienza.

## Verifiche eseguite

Ambiente reale verificato: R 4.5.x, SQMtools 1.7.2 e progetto `in/Au_sip` creato con
SqueezeMeta 1.7.3.alpha3. Il warning di compatibilità già noto è rimasto visibile e
viene registrato nel manifest di corsa.

| Verifica | Esito |
|---|---|
| Suite rapida | PASS, 22 file di test R su 22 |
| Gate di copertura P0-P3 | PASS; tutte le funzioni mirate almeno all'80% |
| Copertura globale informativa P3 | 20,66% |
| Integrazione P0 reale | PASS; Bacillota = 88.258 ORF esatte |
| Integrazione P1 reale | PASS; Top-20 globale/Bacillota e massa FLOW conservata |
| Integrazione P2 reale | PASS; 3 KO multi-EC, 705.510 ORF senza KO, 1.317 ORF multi-KO |
| Integrazione P3 reale | PASS; pathway 00633, `S13_1_8`, K10679: KO TPM 18,284 e pathway TPM 25,632 |
| Edge case reale taxon/pathway | PASS; una ORF conservata, nove pathway vuoti saltati, due artefatti taxonomy prodotti |
| Tre CLI nello stesso output | PASS; tre `run_id` e manifest distinti, ordine CLI e SQM registrato |
| Fallimento controllato | PASS; artefatti parziali e manifest `failed` presenti, exit non zero |
| Igiene del diff | PASS; `git diff --check` senza errori |
| Scansione credenziali nei file modificati | PASS; nessuna corrispondenza |

Tre esempi reali di identità allocate nello stesso output sono
`20260901T120723_UTCp0200_a4f9`, `20260901T120757_UTCp0200_95a4` e
`20260901T121338_UTCp0200_93db`. I relativi manifest riportano rispettivamente
ordine CLI, ordine delle colonne SQM e una terza corsa indipendente; nessun file
manifest è stato sovrascritto o fuso.

La copertura globale resta informativa e inferiore all'80% perché il file contiene
ancora l'intera orchestrazione monolitica; i gate di regressione P0-P3 sulle funzioni
scientifiche e operative selezionate superano invece la soglia richiesta.

Una revisione finale indipendente su due assi, regole del repository e specifica
approvata, ha inizialmente rilevato tre problemi: manifest mancante per un errore CLI
successivo a `output_dir`, escaping insufficiente nel writer TSV usato senza `readr` e
un helper legacy che cercava manifest senza `run_id`. I tre punti sono stati corretti,
coperti da test e riesaminati. La seconda revisione non ha trovato finding residui di
severità critica, alta o media.

## Confronto con `AUDIT_SQM_PLOTS_BUGS.md`

Le integrazioni P0-P3 non mostrano regressioni dirette dei bug precedentemente
chiusi:

| Gruppo | Stato dopo questo intervento |
|---|---|
| BUG-P0-01/02 | mapping KEGG e percentuali taxonomy-by-pathway restano corretti |
| BUG-P0-03 | filtro sull'identità ORF resta corretto; ora copre anche il taxon reale a una ORF |
| BUG-P1-01/02 | selezione pathway-only e Top-20 contestuale restano corrette |
| BUG-P1-03 | join FLOW conserva ancora chiavi, righe e massa TPM |
| BUG-P2-01/02 | `Unclassified`/`Other` e associazioni multi-EC restano distinti/completi |
| BUG-P2-03/04/05 | zeri, validazione CLI e audit multi-KO restano coperti |
| BUG-P3-01/02/03/04/05/06 | integrità artefatti, PIE, Pathview ID, preflight, portabilità e gate restano verdi |

La modifica `ignore_unclassified_functions=FALSE` amplia la robustezza del fix
BUG-P0-03 senza cambiare il confine biologico: gli ID restano quelli risolti sulla
tassonomia ORF e i TPM non vengono riscalati.

## Confronto con `AUDIT_SQM_PLOTS_REVIEW_2026-08-31.md`

Il commit precedente `376fb32` aveva già corretto i finding del 31 agosto relativi a
BRITE in `defined`, denominatore della taxonomy globale, semantica PIE,
provenienza/TSV Pathview, palette SQMtools taxonomy, HTML FLOW, enzimi senza segnale,
dimensioni dei PNG e validazione dei TPM.

Questo intervento completa i punti di provenienza rimasti adiacenti a quell'audit:

- la directory temporanea Pathview è ora accompagnata da nomi file per-corsa e da un
  config fedele alla scala lineare;
- la sorgente Pathview è dichiarata come matrice KO completa e i colori nativi sono
  un'eccezione esplicita nelle regole;
- i manifest non possono più mescolare artefatti validi di esecuzioni diverse;
- i casi vuoti non trasformano un input valido in un fallimento globale;
- l'ordine dei line plot è esplicito e tracciato, senza attribuirgli una semantica
  temporale automatica.

Le scelte su colori Pathview, matrice completa e line plot differiscono dalle
raccomandazioni più restrittive degli audit, ma non sono più ambigue: sono decisioni
approvate, implementate e documentate con provenienza verificabile.

## Limiti residui

- `sqm_plots.R` resta un monolite di oltre 4.600 righe. Il refactoring e il riuso di
  una singola tabella ORF-long per tutti i consumer sono esplicitamente rinviati a un
  branch separato.
- Pathview dipende dal servizio KEGG e dalla propria logica interna di mapping e
  colore. Il TSV completo documenta l'input fornito, non garantisce di ricostruire i
  soli nodi disegnati.
- In questo lotto non è stata rigenerata una nuova esecuzione completa `--mode=all`
  né una chiamata Pathview live. Il contratto Pathview è stato verificato con un
  exporter controllato; le quattro integrazioni reali verificano le invarianti sui
  dati `Au_sip` senza dipendere dalla disponibilità del servizio KEGG.
- I line plot collegano campioni nominali secondo l'ordine scelto o ereditato; la loro
  interpretazione scientifica resta responsabilità dell'utente.
- Gli output storici non vengono approvati retroattivamente. Per una nuova analisi si
  devono usare gli artefatti e i manifest identificati dal nuovo `run_id`.

## Conclusione

Il contratto richiesto è implementato e verificato: Pathview dichiara e usa TPM
lineari con colori nativi, i taxa minimi e i pathway vuoti non bloccano più l'intera
analisi, l'ordine dei campioni è esplicito e ogni corsa possiede file e manifest
separati. Non restano attività del piano corrente oltre al refactoring del monolite,
che rimane intenzionalmente fuori perimetro.
