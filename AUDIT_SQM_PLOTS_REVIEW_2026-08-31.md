# Audit indipendente di `sqm_plots.R` — 2026-08-31

## Perimetro e ordine di lavoro

Oggetti esaminati:

- `sqm_plots.R` (4.113 righe);
- `REGOLE_SCRIPT_R_SQUEEZEMETA.md`;
- progetto SqueezeMeta reale `in/Au_sip`;
- suite in `tests/` e output candidati già presenti sotto `out/`;
- solo dopo la revisione indipendente: `AUDIT_SQM_PLOTS_BUGS.md`.

L'ordine richiesto è stato rispettato nella revisione principale: prima sono stati letti
regole, script e input, sono stati eseguiti test e riproduzioni mirate, poi i risultati
sono stati confrontati con l'audit del 2026-08-30.

## Giudizio sintetico

Le correzioni scientifiche principali del vecchio audit non risultano regredite:
mapping KEGG curati, percentuali tassonomiche per pathway, subset ORF per taxon,
`top20` pathway-only e contestuale, conservazione TPM nel join FLOW, multi-EC,
multi-KO e distinzione `Unclassified`/`Other` superano i test correnti.

Lo script non è però ancora approvabile come pienamente conforme. Sono stati
confermati:

- 4 rilievi ad alta severità;
- 5 rilievi medi;
- 3 rilievi bassi/di hardening;
- 1 debito strutturale rilevante.

Il problema numerico più importante è in `taxonomy_global/percent`: sul campione
reale `S13_1_8` la tabella esportata somma a **21,18628**, non a 100, perché
`Unmapped` e `Unclassified` vengono esclusi senza riscalare e senza dichiarare il
denominatore. Il problema concettuale più netto è che la selezione `defined` può
accettare una categoria BRITE come `Transporters` e chiamarla pathway: sul dataset
reale vengono selezionate **28.578 ORF**.

I TPM delle analisi pathway corrette già coperte dall'audit precedente non hanno
mostrato nuovi errori. Alcuni nuovi finding riguardano invece selezione, semantica
delle percentuali, provenienza dei file e metadati: numeri internamente coerenti
possono quindi essere presentati con un significato o una tracciabilità sbagliati.

## Findings confermati

### NEW-P1-01 — `defined` accetta categorie BRITE/non-pathway come pathway

**Severità: alta. Stato rispetto al vecchio audit: nuovo, adiacente a BUG-P1-01.**

Riferimenti:

- `sqm_plots.R:753-823`, `resolve_pathways()`;
- `sqm_plots.R:826-920`, parser gerarchico corretto usato invece da `top20`;
- `sqm_plots.R:1128-1135`, filtro finale con `subsetFun()`.

Meccanismo:

`resolve_pathways()` costruisce i candidati prendendo l'ultima componente di ogni
voce di `sqm$misc$KEGG_paths`, senza richiedere una gerarchia PATHWAY valida a tre
livelli e senza controllare `is_pathway_map`. Il filtro pathway-only esiste, ma viene
applicato soltanto in `select_top_pathways()`.

Riproduzione reale:

```text
requested = Transporters
resolved$canonical_pathway_name = Transporters
resolved$pathway_id = NA
subset ORF = 28578
```

Conseguenza:

FUNZ, FLOW, TAXON e PIE possono produrre directory, titoli e manifest che chiamano
“pathway” una categoria BRITE. I TPM possono essere sommati correttamente, ma
l'identità biologica dell'oggetto analizzato è falsa.

Correzione:

Costruire anche i candidati `defined` da `parse_kegg_pathway_membership()` filtrata
con `is_pathway_map`, verificare l'univocità della gerarchia e rifiutare
esplicitamente le foglie BRITE. Aggiungere una regressione reale/sintetica con
`Transporters` che deve fallire come non-pathway.

### NEW-P1-02 — `taxonomy_global/percent` ha denominatore non dichiarato e non somma a 100

**Severità: alta. Stato: nuovo, adiacente a BUG-P0-02 e BUG-P2-01.**

Riferimenti:

- `sqm_plots.R:2110-2135`, chiamata a `SQMtools::plotTaxonomy()`;
- `sqm_plots.R:2138-2155`, export generico `taxon/value/count`;
- `sqm_plots.R:3967-3988`, scope globale con `ignore_unmapped=TRUE` e
  `ignore_unclassified=TRUE`;
- regole, righe 240-257: denominatore identificabile e controllo delle somme.

Meccanismo:

SQMtools 1.7.2 usa `rescale=FALSE` per default. Lo script rimuove `Unmapped` e
`Unclassified`, ma non riscalda le categorie rimaste. Il TSV chiama la misura
semplicemente `value`, imposta `count=percent` e non esporta denominatore, stato o
quota esclusa.

Riproduzione sul dato reale con `rank=phylum`, `top_n_taxa=3`, `S13_1_8`:

```text
taxa esportati = Bacillota, Nitrososphaerota, Pseudomonadota, Other
somma value = 21.18628
quota non mostrata = 78.81372 punti percentuali
```

Conseguenza:

Il risultato non è una composizione tassonomica dei taxa mostrati, benché il nome
`percent` possa farlo credere. Inoltre `Unclassified`, informazione distinta secondo
le regole, scompare senza una traccia nel TSV o nel manifest.

Correzione:

Scegliere e dichiarare una semantica:

1. composizione dei taxa mostrati: conservare `Unclassified` e/o riscalare a 100;
2. quota sul totale della libreria: mantenere il denominatore originale, ma esportare
   `total_sample_abundance`, `excluded_unmapped`, `excluded_unclassified`, una
   colonna dal nome esplicito e una nota nel manifest/grafico.

In entrambi i casi aggiungere una verifica per sample coerente con il denominatore.

### NEW-P1-03 — PIE ignora `top_n_ko` ma lo dichiara nei manifest

**Severità: alta. Stato: nuovo.**

Riferimenti:

- `sqm_plots.R:3248-3266`, parametro `top_n_ko`;
- `sqm_plots.R:3296-3304`, selezione di tutti i KO e loop completo;
- `sqm_plots.R:3343-3403`, manifest con `top_n_ko` valorizzato.

Evidenza reale in `out/p3_candidate_01_all`:

```text
top_n_ko dichiarato = 5
directory KO prodotte per S13_1_8 = 21
righe manifest PIE per S13_1_8 = 42  # TSV + PNG, un solo rank/dimensione
```

Conseguenza:

Il manifest dichiara una selezione non applicata; l'utente può interpretare i file
come “Top 5” e ottiene molti più output del previsto. I valori all'interno del singolo
pie restano coerenti, ma selezione e provenienza dichiarate sono false.

Correzione:

Decidere esplicitamente tra:

- PIE per tutti i KO: rimuovere `top_n_ko` dal manifest PIE o valorizzarlo `NA/all` e
  documentarlo nell'help;
- PIE Top N: graduatoria deterministica sul TPM aggregato dei campioni mostrati,
  coerente con FUNZ/FLOW. Va definito se i KO fuori Top N sono saltati o rappresentati
  da un oggetto `Other`, perché un pie di `Other KO` avrebbe semantica diversa.

### NEW-P1-04 — Pathview può attribuire alla run corrente file vecchi o campioni sbagliati

**Severità: alta. Stato: nuovo, non una regressione diretta di BUG-P3-01.**

Riferimenti:

- `sqm_plots.R:3198`, snapshot dei file preesistenti;
- `sqm_plots.R:3214-3218`, fallback da delta vuoto a tutti i file presenti;
- `sqm_plots.R:3220-3241`, costruzione del manifest;
- `sqm_plots.R:3229`, ogni file riceve sempre l'intero vettore `selected_samples`.

Meccanismo:

Quando `exportPathway()` sovrascrive file con gli stessi nomi, `setdiff(after,before)`
è vuoto. Il codice interpreta allora **tutti** i file della directory come output
correnti. Un file lasciato da una run precedente entra nel nuovo manifest perché
esiste ed è non vuoto; il pruning di BUG-P3-01 non può riconoscerlo come stale.

Caso minimo riprodotto:

1. directory con `map.png` e `stale_old.png`;
2. exporter fittizio che riscrive solo `map.png`;
3. il manifest corrente contiene entrambi i file.

In modalità `separato`, inoltre, file specifici di un sample ricevono comunque tutti
i `selected_samples`, rendendo imprecisa l'attribuzione per-file.

Correzione:

Eseguire ogni export in una directory temporanea/run-specific e manifestare solo i
file prodotti lì, quindi spostare i file noti nella destinazione; in alternativa
usare una lista deterministica di target più hash/mtime prima-dopo. Per
`split_samples=TRUE`, derivare il sample del singolo file o produrre una directory e
un manifest distinti per sample.

### NEW-P2-01 — I PNG Pathview non hanno il TSV sorgente richiesto

**Severità: media. Stato: nuovo, adiacente a BUG-P3-02.**

Riferimenti:

- `sqm_plots.R:3154-3246`, nessuna chiamata a `write_tsv_safe()`;
- regole, righe 268-276: ogni grafico deve avere il TSV esatto;
- `out/p3_candidate_01_all/pathview/`: due PNG e il manifest, nessun TSV dati.

Conseguenza:

Il colore e il valore assegnato ai KO nel diagramma non sono ricostruibili o
verificabili dai soli output. Il manifest Pathview ha inoltre i campi di audit KO a
`NA`, quindi non collega il grafico alla policy multi-KO dichiarata altrove.

Correzione:

Esportare la stessa matrice KO × sample, con gli stessi valori e trasformazioni,
passata a `exportPathway()/pathview`, più eventuali limiti/color bins. Il TSV deve
essere scritto prima del grafico e il grafico deve essere ricostruibile da esso.

### NEW-P2-02 — La palette obbligatoria viene scartata nei grafici SQMtools

**Severità: media. Stato: nuovo.**

Riferimenti:

- `sqm_plots.R:90-106`, palette da 66 elementi;
- `sqm_plots.R:2123-2132`, passaggio integrale della palette a `plotTaxonomy()`;
- regole, righe 289-308: palette obbligatoria.

SQMtools 1.7.2 accetta `color` in questo percorso solo se `length(color) == N`.
Con `top_n_taxa=3`, la riproduzione corrente emette:

```text
You passed less/more colors than taxa. Using default colors
```

Quindi la palette del progetto non viene usata per i grafici globali e per i grafici
pathway `count=abund` (salvo il caso accidentale `top_n_taxa=66`). Colori assegnati
da un default e da ordini diversi rendono meno affidabile il confronto visivo tra
grafici.

Correzione:

Passare esattamente `N` colori dalla palette, tenendo separati i colori riservati, o
costruire il grafico direttamente dal TSV con una mappa nome→colore deterministica.
Aggiungere un test che ispezioni la scala, non soltanto la classe `ggplot`.

### NEW-P2-03 — Gli HTML FLOW non sono trasferibili/verificabili dal manifest

**Severità: media. Stato: nuovo.**

Riferimenti:

- `sqm_plots.R:521-524`, `saveWidget(..., selfcontained=FALSE)`;
- `sqm_plots.R:2942-2969`, il manifest registra soltanto il file `.html`;
- `sqm_plots.R:3451-3500`, la validazione controlla solo il target dichiarato.

`selfcontained=FALSE` crea una directory `*_files` con dipendenze JavaScript/CSS.
Questi file non sono manifestati. Un HTML può quindi risultare valido nel manifest
ma non funzionare se la sidecar manca; copiare i soli `output_file` non produce un
output trasferibile.

Correzione:

Usare un HTML self-contained quando l'ambiente lo consente, oppure manifestare ogni
dipendenza con relazione al file principale e validare l'intero bundle.

### NEW-P2-04 — Le combinazioni EC senza segnale producono grafici vuoti senza warning/skip

**Severità: media-bassa. Stato: nuovo, adiacente a BUG-P2-03.**

Riferimenti:

- `sqm_plots.R:1687-1740`, griglia completa sample × EC con zeri;
- `sqm_plots.R:2704-2802`, export incondizionato di TSV e PNG;
- regole, righe 321-323: warning e skip della combinazione priva di righe positive.

Nel candidato reale con i 20 EC di default, **8 EC su 20** hanno TPM totale zero,
ma ricevono comunque directory, TSV e grafici separati. Con i default completi ciò
può generare decine di PNG vuoti senza segnalare che l'enzima non è stato osservato.

Correzione:

Conservare, se utile, una riga TSV `status=zero_denominator/no_positive_tpm`, emettere
warning e saltare i PNG della combinazione EC. Il grafico combinato può mantenere o
escludere gli EC zero, ma la scelta deve essere esplicita.

### CONCEPT-P2-05 — Il line plot enzimatico implica un ordine dei campioni non validato

**Severità: media concettuale. Stato: nuovo.**

Riferimenti:

- `sqm_plots.R:123`, `line` è un tipo di grafico predefinito;
- `sqm_plots.R:1770-1786`, i campioni vengono collegati nell'ordine CLI/tabella.

Una linea tra `S13_1_8`, `S13_2_8` e `S13_3_8` suggerisce continuità o trend. Lo
script conosce soltanto nomi categorici; l'ordine di colonna non prova una relazione
temporale o quantitativa.

Correzione:

Rendere `line` opt-in e richiedere una variabile x ordinata/metadata esplicita; il
default sicuro per campioni nominali è il barplot o il point plot non connesso.

### NEW-P3-01 — I manifest PNG enzimatici perdono larghezza e altezza

**Severità: bassa. Stato: nuovo.**

Riferimenti:

- `sqm_plots.R:2679-2700`, helper manifest senza `width`/`height`;
- `sqm_plots.R:2723-2746` e `2774-2800`, loop su `unname(files)` che perde il nome
  della dimensione;
- regole, righe 285-287: dimensioni e DPI quando pertinenti.

Prova nel candidato reale:

```text
output_file = funz/enzimi/insieme/barplot_enzimi_6x4.png
width = NA
height = NA
dpi = 150
```

Correzione: iterare sui nomi di `png_files` come negli altri modi e passare le
dimensioni effettive a `new_manifest_row()`.

### NEW-P3-02 — Validazione input incompleta

**Severità: bassa/hardening. Stato: nuovo.**

Problemi confermati:

- `tax_mode` non viene validato contro `prokfilter/allfilter/nofilter` prima di creare
  `output_dir` e chiamare `loadSQM()` (`sqm_plots.R:3690,3754,3769`);
- `--nome` consuma come valore anche il token successivo `--altra_opzione`
  (`sqm_plots.R:183-190`), producendo un errore indiretto;
- `plot_dpi` e dimensioni non rifiutano `Inf` perché manca `is.finite()`
  (`sqm_plots.R:366-386,3706-3712`);
- i nomi sample sono usati senza sanitizzazione nei filename FLOW
  (`sqm_plots.R:2872,2906,2943`), anche se i sample reali correnti sono sicuri.

Questi casi normalmente terminano con errore e non hanno falsato il dataset attuale,
ma violano il fail-fast e possono lasciare output parziali o messaggi fuorvianti.

### NEW-P3-03 — Due aggregazioni bypassano la validazione dei TPM

**Severità: bassa/hardening; potenziale severità alta su input corrotto. Stato: nuovo.**

Riferimenti:

- `sqm_plots.R:988-1008`, ranking `top20`;
- `sqm_plots.R:1723-1732`, aggregazione enzimi.

Entrambi usano `sum(as.numeric(tpm), na.rm=TRUE)` senza verificare conversione,
finitudine e non negatività. Un valore testuale può diventare `NA` ed essere ignorato;
un `Inf` può dominare il ranking; un negativo può cancellare TPM positivi.

Sul dataset reale corrente non sono state osservate celle negative o non numeriche,
quindi il finding è preventivo e non un falso risultato già prodotto.

Correzione: un'unica validazione numerica subito dopo `loadSQM()`, prima di qualsiasi
ranking/aggregazione, con errore esplicito su non numerico, non finito o negativo.

## Debito strutturale

### DEBT-01 — Monolite da 4.113 righe e ricostruzione ripetuta della tabella centrale

**Severità: media per manutenibilità e rischio di divergenza. Stato: già noto in forma
generale nel vecchio audit, non risolto.**

Lo script contiene nello stesso file CLI, importazione, modello dati, selezione KEGG,
sei famiglie di grafici, Pathview e riconciliazione manifest. `main()` occupa circa
465 righe (`3648-4113`). In `mode=all`, FUNZ, FLOW e PIE richiamano separatamente
`build_orf_long_result()` per lo stesso pathway, invece di costruire una volta la
tabella centrale richiesta dalle regole e passarla ai consumer.

Conseguenze:

- costo ripetuto su subset grandi;
- maggiore rischio che modalità diverse applichino controlli o denominatori diversi;
- test mirati verdi ma copertura globale storicamente intorno al 15%, molto sotto il
  requisito progettuale dell'80%;
- i nuovi finding si concentrano proprio nei percorsi di orchestrazione non coperti.

Direzione: separare almeno `cli`, `sqm_validation`, `pathway_selection`,
`orf_long`, `plots/*`, `manifest`; calcolare/cacheare `orf_long_result` per
context/pathway/sample-set e far dipendere i modi dalla stessa struttura validata.

## Confronto con `AUDIT_SQM_PLOTS_BUGS.md`

Non sono state trovate regressioni dirette dei 17 bug precedenti. Alcune correzioni
restano però incomplete fuori dal perimetro esatto dei vecchi test.

| Finding precedente | Stato corrente | Nota |
|---|---|---|
| BUG-P0-01 mapping KEGG | Risolto | `00710`, `00633`, `00910` sono corretti. |
| BUG-P0-02 percentuali taxonomy pathway | Risolto | La taxonomy globale resta problematica: NEW-P1-02. |
| BUG-P0-03 filtro taxon contig-based | Risolto | `subsetORFs(..., tax_source="orfs")`; Bacillota = 88.258 ORF esatte. |
| BUG-P1-01 `top20` con BRITE | Risolto nel perimetro `top20` | La modalità `defined` ha la lacuna gemella NEW-P1-01. |
| BUG-P1-02 `top20` globale anziché contestuale | Risolto | Ranking eseguito dentro ogni `filter_context`. |
| BUG-P1-03 join FLOW molti-a-molti | Risolto | Metadata KO uno-a-uno e invariante TPM conservata. |
| BUG-P2-01 `Unclassified` in `Other` | Risolto in FLOW/PIE | La globale esclude `Unclassified`: NEW-P1-02, meccanismo diverso. |
| BUG-P2-02 multi-EC | Risolto | Lookup deterministico e join senza aumento di massa. |
| BUG-P2-03 FUNZ con sample zero | Risolto in FUNZ | ENZIMI non adotta la stessa policy: NEW-P2-04. |
| BUG-P2-04 interi CLI frazionari | Risolto | Regex intera prima della conversione. |
| BUG-P2-05 policy multi-KO non tracciata | Risolto per FUNZ/FLOW/PIE | Pathview non ha TSV/audit: NEW-P2-01. |
| BUG-P3-01 manifest con target mancanti | Risolto per esistenza/integrità | Pathview può attribuire un file vecchio ma esistente alla run: NEW-P1-04. |
| BUG-P3-02 TSV PIE incompleto | Risolto | Round-trip TSV→grafico verificato. |
| BUG-P3-03 Pathview con ID `NA` | Risolto | Gate su ID numerico valido e skip controllato. |
| BUG-P3-04 preflight package | Risolto | Preflight mode-aware; test dedicato verde. |
| BUG-P3-05 path personale nel test | Risolto | Fixture Windows sintetico. |
| BUG-P3-06 suite incapace di rilevare falsi risultati | Migliorato, non completo | Mancano regressioni per i nuovi finding. |

## Verifiche eseguite

Ambiente:

- R 4.5.0 tramite `C:\Progra~1\R\R-45~1.0\bin\Rscript.exe`;
- SQMtools 1.7.2;
- warning atteso e lasciato visibile: progetto SqueezeMeta 1.7.3.alpha3 caricato con
  SQMtools 1.7.2.

Esiti:

| Verifica | Esito/evidenza |
|---|---|
| 19 test veloci | PASS |
| Integrazione P0 | PASS; Bacillota 88.258 ORF; 3 sample pathway |
| Integrazione P1 | PASS; global top20=20; Bacillota top20=20; TPM FLOW conservato |
| Integrazione P2 | PASS; 3 KO multi-EC nel 00710; 705.510 ORF senza KO; 1.317 ORF multi-KO |
| Integrazione P3 | PASS; 00633/S13_1_8/K10679: KO TPM 18,284; pathway TPM 25,632 |
| Categoria BRITE in `defined` | FAIL atteso ma accettata; `Transporters`, 28.578 ORF |
| Taxonomy globale percentuale | FAIL semantico; somma 21,18628 e nessun denominatore esportato |
| Palette SQMtools | FAIL; warning e fallback alla palette default |
| PIE `top_n_ko=5` | FAIL; 21 KO reali esportati per il sample esaminato |
| Enzimi di default | 8/20 EC con TPM zero ma grafici prodotti |
| Manifest PNG enzimi | width/height `NA` su file dimensionato |

La run completa `--mode=all` non è stata rigenerata durante questo audit: il
candidato reale già presente e le quattro integrazioni erano sufficienti per
riprodurre i problemi, mentre una nuova run completa produrrebbe migliaia di file
senza aggiungere evidenza ai meccanismi confermati.

## Affidabilità degli output esistenti

Non approvare retroattivamente l'intero `out/p3_candidate_01_all` come pienamente
conforme:

- i TSV pathway percentuali, i subset Bacillota e i FLOW coperti dalle integrazioni
  restano supportati dalle verifiche;
- `taxonomy_global/percent` ha semantica/denominatore non dichiarati;
- i PIE non rispettano il `top_n_ko` scritto nel manifest;
- i Pathview non hanno TSV sorgente e la strategia manifest non garantisce la
  provenienza dalla run;
- i grafici taxonomy basati su SQMtools non usano la palette prescritta;
- i bundle HTML non sono completamente inventariati;
- gli output enzimi includono combinazioni zero e manifest senza dimensioni.

## Priorità di correzione proposta

1. Bloccare categorie BRITE/non-pathway anche in `defined`.
2. Definire e rendere esplicito il denominatore di `taxonomy_global/percent`.
3. Rendere coerenti PIE e `top_n_ko`.
4. Rendere Pathview run-isolated, sample-specific e dotato di TSV sorgente.
5. Correggere palette taxonomy e bundle HTML.
6. Applicare warning/skip agli enzimi senza segnale e riconsiderare il line plot
   predefinito.
7. Chiudere manifest dimensioni e validazioni CLI/TPM.
8. Aggiungere test RED per ogni finding prima del fix e poi spezzare il monolite per
   ridurre la divergenza tra modalità.

## Gate minimo dopo le correzioni

- `Transporters` e altre foglie BRITE rifiutate in `defined`;
- ogni TSV `percent` espone il denominatore e supera la somma attesa per sample;
- numero KO PIE coerente con la semantica dichiarata nel manifest;
- nessun file Pathview preesistente entra nel manifest corrente;
- TSV Pathview ricostruisce i valori del grafico;
- palette effettiva verificata dalla scala ggplot;
- HTML funzionante copiando tutti e soli i target manifestati;
- nessun PNG per combinazioni `status=no_positive_tpm`;
- larghezza/altezza presenti per ogni PNG;
- `tax_mode`, dimensioni, DPI e TPM non finiti rifiutati prima di produrre output;
- 19 test veloci, 4 integrazioni e nuovi test di regressione tutti verdi.
