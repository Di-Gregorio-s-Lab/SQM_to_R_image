# Audit indipendente di `sqm_plots.R` — 2026-09-01

## Perimetro, ordine e snapshot

L'ordine richiesto è stato rispettato:

1. lettura integrale di `sqm_plots.R` e `REGOLE_SCRIPT_R_SQUEEZEMETA.md`;
2. ispezione del progetto reale `in/Au_sip`, della suite e delle API installate di SQMtools 1.7.2;
3. test e riproduzioni mirate senza consultare gli audit precedenti;
4. solo dopo la chiusura della ricerca indipendente, lettura integrale e confronto con
   `AUDIT_SQM_PLOTS_BUGS.md` e `AUDIT_SQM_PLOTS_REVIEW_2026-08-31.md`.

Lo script corrente ha 4.604 righe. Il workspace era già modificato prima di questo
audit; il riferimento riproducibile è quindi l'hash del file analizzato, non `HEAD`:

| File | SHA-256 |
|---|---|
| `sqm_plots.R` | `0ee9487646d4dc1fe7901b508b28e2e18118d949633be3b5bfda84d10cd92473` |
| `REGOLE_SCRIPT_R_SQUEEZEMETA.md` | `99a87109cb55f17363c9751a1d4a2ea32e9438419e759dec4c3b0ddb018f2f4c` |
| `AUDIT_SQM_PLOTS_BUGS.md` | `bea0e1c6a2e2fed2c199288048037b01de747284e6bb6ef0b02a7755dca1c1f6` |
| `AUDIT_SQM_PLOTS_REVIEW_2026-08-31.md` | `bedfc1aba6c089499d8d243ee340084d348d0a24b7ae1e24bdd31a00cdadba5a` |

Ambiente verificato:

- R 4.5.0;
- SQMtools 1.7.2;
- progetto creato con SqueezeMeta 1.7.3.alpha3;
- il warning di compatibilità atteso è rimasto visibile;
- dataset reale: 945.388 ORF e tre campioni (`S13_1_8`, `S13_2_8`, `S13_3_8`).

## Giudizio sintetico

Le correzioni scientifiche centrali dei due audit precedenti non risultano regredite
nei casi coperti: mapping KEGG, selezione pathway-only, Top-20 contestuale, percentuali
tassonomiche per pathway, massa FLOW, multi-KO, multi-EC e distinzione
`Unclassified`/`Other` superano i test sintetici e reali.

Lo script non è ancora pienamente approvabile. Restano aperti:

- 1 finding alto: Pathview esporta un TSV di configurazione falso rispetto alla
  trasformazione realmente applicata;
- 2 finding medio-alti: filtri tassonomici validi e combinazioni pathway valide ma
  vuote possono arrestare un'intera esecuzione;
- 3 finding medi: sorgente/palette Pathview incomplete, manifest cumulativi e line
  plot enzimatico semanticamente non giustificato;
- 2 finding bassi di robustezza CLI/modello dati;
- 1 debito strutturale medio già noto e cresciuto.

Non sono emersi nuovi falsi TPM nei percorsi FUNZ, FLOW, PIE e taxonomy-by-pathway
coperti dalle quattro integrazioni. Il rischio residuo più netto riguarda il significato
e la riproducibilità di Pathview, non le somme TPM già validate negli altri modi.

## Findings aperti

| ID | Severità | Stato rispetto agli audit precedenti | Area |
|---|---|---|---|
| AUD-2026-09-01-P1-01 | alta | nuovo; correzione parziale di NEW-P2-01 | Pathview |
| AUD-2026-09-01-P1-02 | medio-alta | nuovo; adiacente a BUG-P0-03 | filtro taxon |
| AUD-2026-09-01-P1-03 | medio-alta | nuovo | orchestrazione pathway |
| AUD-2026-09-01-P2-01 | media | parzialmente risolto rispetto a NEW-P2-01/02 | Pathview |
| AUD-2026-09-01-P2-02 | media | nuovo; adiacente a BUG-P3-01 | manifest |
| AUD-2026-09-01-P2-03 | media concettuale | già noto, non risolto (CONCEPT-P2-05) | enzimi |
| AUD-2026-09-01-P3-01 | bassa/hardening | nuova estensione di NEW-P3-02 | CLI |
| AUD-2026-09-01-P3-02 | bassa/hardening | nuovo | fallback KEGG |
| DEBT-01 | media strutturale | già noto, non risolto | architettura |

### AUD-2026-09-01-P1-01 — Pathview dichiara una scala logaritmica che non usa

**Severità: alta. Tipo: metadato falso e riproduzione non fedele.**

Riferimenti correnti:

- `sqm_plots.R:3520-3528`, chiamata reale a `SQMtools::exportPathway()`;
- `sqm_plots.R:3653-3665`, scrittura di `pathview_render_config.tsv`.

Il TSV dichiara:

```text
metric = tpm
log_scale = TRUE
pseudocount = 0.001
color_bins = 10
```

La chiamata reale passa soltanto `count="tpm"`, `samples`, `split_samples`, directory e
suffisso. L'introspezione diretta di SQMtools 1.7.2 conferma questi default:

```r
log_scale = FALSE
sample_colors = NULL
max_scale_value = NULL
color_bins = 10
```

Quindi il grafico usa TPM su scala lineare, mentre il file che dovrebbe descriverlo
afferma una trasformazione `log10(TPM + 0.001)`. Una ricostruzione basata sul TSV
produce colori e soglie diversi dal grafico reale. Il difetto è silenzioso.

Correzione richiesta:

- decidere se il rendering desiderato è lineare o logaritmico;
- passare esplicitamente a `exportPathway()` tutti i parametri registrati, oppure
  registrare i valori effettivi (`log_scale=FALSE` e pseudocount non applicabile);
- aggiungere un test con exporter fittizio che confronti argomenti catturati e TSV
  campo per campo. Il test corrente verifica l'esistenza del config, non la verità dei
  suoi valori.

### AUD-2026-09-01-P1-02 — Un taxon valido con funzioni non classificate può fallire

**Severità: medio-alta. Tipo: combinazione valida rifiutata.**

Riferimento: `sqm_plots.R:1333-1363`, in particolare
`ignore_unclassified_functions = TRUE` nella chiamata a `SQMtools::subsetORFs()`.

Riproduzione sul dataset reale:

```text
taxon = Accipitriformes (no class in NCBI)
rank risolto = class
ORF attese = 1
ignore_unclassified_functions=TRUE  -> ERROR: incorrect number of dimensions
ignore_unclassified_functions=FALSE -> OK: 1 ORF
```

Il taxon viene risolto correttamente e non è ambiguo; il fallimento nasce dalla
combinazione scelta dallo script con SQMtools 1.7.2. Il caso Bacillota dell'audit
precedente resta verde, quindi non è una regressione del filtro ORF-vs-contig: è un
edge case non coperto della correzione.

Impatto:

- un filtro tassonomico valido arresta la run prima di produrre manifest;
- il problema interessa anche modalità che non richiedono funzioni classificate,
  per esempio una tassonomia globale del taxon;
- il messaggio interno di SQMtools non identifica il parametro responsabile.

Correzione richiesta:

- usare `ignore_unclassified_functions=FALSE`, coerente con la policy generale di
  non perdere ORF non annotate, oppure applicare una gestione condizionale testata;
- aggiungere una fixture a una ORF senza funzione e una regressione reale o sintetica
  che verifichi l'identità degli `orf_id` restituiti.

### AUD-2026-09-01-P1-03 — Un pathway valido ma vuoto arresta tutta la run

**Severità: medio-alta. Tipo: orchestrazione ed error handling errati.**

Riferimenti:

- `sqm_plots.R:1277-1285`, `subsetFun()` viene chiamata senza `allow_empty`, quindi
  usa il default `FALSE`;
- `sqm_plots.R:4350-4369`, tutti i pathway vengono subsettati anticipatamente con
  `map()`, senza `tryCatch`, controllo di vuoto o skip selettivo;
- regole, righe 335-337: warning e skip della sola combinazione senza righe positive.

Riproduzione reale: il filtro valido `Archaeoglobi (no class in NCBI)` produce uno SQM
non vuoto, ma i pathway definiti predefiniti non vi hanno ORF. La prima chiamata a
`subsetFun(..., allow_empty=FALSE)` solleva un errore e interrompe l'intero contesto.

L'impatto è maggiore del singolo grafico: la costruzione dei subset avviene prima dei
rami di modalità. Anche `--mode=taxon`, che potrebbe produrre la tassonomia globale
del filtro, può essere bloccato da un pathway definito assente.

Correzione richiesta:

- costruire i subset solo per le modalità che li consumano;
- usare `allow_empty=TRUE`, verificare subito le ORF/il segnale e scartare con warning
  soltanto `context × pathway` vuoto;
- serializzare lo skip in un audit/manifest di esecuzione, senza inventare un grafico;
- aggiungere regressioni per taxon valido con zero ORF nel pathway e per
  `mode=taxon`, che deve comunque produrre lo scope globale.

### AUD-2026-09-01-P2-01 — La correzione Pathview resta incompleta

**Severità: media. Stato: parziale rispetto a NEW-P2-01 e NEW-P2-02.**

Aspetti positivi già corretti:

- export eseguito in directory temporanea isolata (`sqm_plots.R:3509-3559`);
- TSV input e TSV configurazione presenti;
- manifest per invocazione `insieme`/`separato` e target reali non vuoti.

Lacune residue:

1. `build_pathview_input_table()` (`3492-3506`) esporta **tutti** i KO presenti in
   `sqm$functions$KEGG$tpm`, mentre il grafico usa solo i KO mappati sui nodi del
   pathway. Il TSV è quindi un superset, non la tabella esatta usata dal grafico come
   richiesto dalle regole.
2. Il renderer riceve lo SQM direttamente, non il TSV scritto. Oggi i valori di base
   provengono dalla stessa matrice, ma non esiste una postcondizione che dimostri che
   nodi, trasformazioni, soglie e colori del grafico siano ricostruibili dal TSV.
3. Non viene passato `sample_colors`. SQMtools 1.7.2 usa per default il rosso per i
   campioni; la palette obbligatoria del progetto non è applicata. L'eccezione aggiunta
   alle regole riguarda solo `SQMtools::plotTaxonomy()` e, per espressa previsione
   delle righe 320-322, non gli altri grafici.
4. I campi di audit multi-KO restano `NA` per Pathview, nonostante il grafico usi la
   matrice KO aggregata.

Correzione richiesta:

- esportare una tabella nodo/pathway × KO × sample con il valore effettivo, il valore
  trasformato, il bin colore e l'indicazione `mapped/plotted`;
- passare esplicitamente `sample_colors` dalla palette e registrare la mappa sample →
  colore;
- verificare con un exporter controllato che ogni valore visualizzato sia presente
  nel TSV e che nessun KO estraneo sia marcato come plottato.

### AUD-2026-09-01-P2-02 — I manifest di sezione mescolano esecuzioni diverse

**Severità: media. Tipo: provenienza e semantica del manifest.**

Riferimenti: `sqm_plots.R:4022-4041` e `4055-4068`.

`write_section_manifest()` legge il manifest preesistente, elimina soltanto i target
mancanti/vuoti e poi unisce tutte le righe ancora valide a quelle correnti. La
deduplicazione è solo per `output_file`.

Conseguenza: un file prodotto da una run precedente, ancora esistente ma non generato
dalla run corrente, resta nel manifest insieme a righe nuove. Lo stesso TSV può quindi
mescolare progetti, campioni, `tax_mode`, dimensioni e selezioni differenti senza
`run_id` o timestamp. Questo non è più il bug storico dei target mancanti
(BUG-P3-01), che risulta corretto; è un problema distinto: il file è valido ma viene
attribuito implicitamente al manifest aggiornato.

La regola richiede un manifest per ogni esecuzione. L'implementazione corrente produce
invece un inventario cumulativo della directory.

Correzione richiesta:

- separare `manifest_current_run.tsv` da un eventuale `artifact_inventory.tsv`;
- aggiungere `run_id`, timestamp e hash/config della run;
- non fondere righe storiche nel manifest della run corrente;
- mantenere il pruning degli artefatti mancanti solo nell'inventario cumulativo.

### AUD-2026-09-01-P2-03 — Il line plot enzimatico resta concettualmente improprio

**Severità: media concettuale. Stato: già noto, non risolto.**

`default_enzyme_plot_types <- c("bar", "line")` è ancora presente a
`sqm_plots.R:138`; `make_enzyme_lineplot()` collega i campioni nell'ordine ricevuto
(`1940-1963`). Lo script non carica metadati temporali né una variabile x quantitativa.

Una linea tra campioni nominali suggerisce continuità o trend non dimostrati. È lo
stesso CONCEPT-P2-05 dell'audit del 31 agosto.

Correzione richiesta: rendere `line` opt-in e richiedere metadati ordinati espliciti;
in assenza di questi usare barplot o punti non connessi.

### AUD-2026-09-01-P3-01 — Liste CLI vuote e campioni duplicati non sono rifiutati

**Severità: bassa/hardening.**

`split_csv_arg()` rimuove le stringhe vuote (`267-270`). Di conseguenza:

- `--taxonomy_counts=` produce `character(0)` e supera
  `all(x %in% allowed)` per verità vacua (`4201-4208`);
- `--flowplot_formats=` ha lo stesso problema (`4210-4217`);
- `--dimensions=` restituisce una lista vuota (`412-435`) e può completare senza
  generare PNG né segnalare il valore mancante.

Inoltre `validate_samples()` controlla solo l'appartenenza. Una lista duplicata come
`S1,S1` passa la validazione e una riproduzione sintetica fallisce più tardi con
`Controllo percentuali fallito per sample`, non con un errore CLI utile.

Correzione richiesta: per ogni lista non opzionale dopo l'esplicitazione CLI,
richiedere `length > 0`, valori univoci e nessun token vuoto; aggiungere test subprocess
che dimostrino il fallimento prima di `loadSQM()` e prima di creare `output_dir`.

### AUD-2026-09-01-P3-02 — Il fallback del nome KO fallisce se `KEGG_names` manca

**Severità: bassa/hardening; nessun impatto sul dataset corrente.**

Riferimenti: `sqm_plots.R:1415-1471` e `1488-1495`.

`validate_sqm_object()` non richiede `sqm$misc$KEGG_names`, coerentemente con il fatto
che è un fallback. Tuttavia, se `KEGGFUN` è vuoto e `KEGG_names` è `NULL`, l'indicizzazione
produce `character(0)` e `dplyr::if_else()` non può riciclarla. Riproduzione sintetica:

```text
ERROR: Can't recycle `true` (size 0) to size 1.
```

Correzione richiesta: costruire sempre un vettore fallback della lunghezza delle righe,
usando prima `KEGG_names` e poi `ko_id`; aggiungere una fixture senza `misc` e con
`KEGGFUN` vuoto.

### DEBT-01 — Monolite e tabella centrale ricostruita tre volte

**Severità: media strutturale. Stato: già noto, peggiorato dimensionalmente.**

Lo script è cresciuto da 4.113 a 4.604 righe. In `mode=all`, lo stesso
`build_orf_long_result()` viene richiamato separatamente da:

- FUNZ, riga 2829;
- FLOW, riga 3138;
- PIE, riga 3744.

Questo contraddice la regola di costruire una sola tabella lunga centrale e aumenta
costo, pressione di memoria e rischio di divergenza tra consumer. Le baseline di
copertura globale osservate restano basse (10,13%-18,50% a seconda dell'esercizio),
pur con coperture mirate sopra l'80%.

Direzione: separare CLI, validazione SQM, selezione pathway, modello ORF-long, renderer
e manifest; costruire/cacheare una sola struttura immutabile per
`context × pathway × sample-set` e passarla ai modi.

## Confronto con gli audit precedenti

### `AUDIT_SQM_PLOTS_BUGS.md` (2026-08-30)

Non sono state osservate regressioni dirette dei 17 bug chiusi.

| Gruppo precedente | Stato corrente | Evidenza sintetica |
|---|---|---|
| BUG-P0-01/02 | risolto | mapping `00710/00633/00910`; percentuali pathway sommano 100 |
| BUG-P0-03 | risolto nel caso principale | Bacillota = 88.258 ORF esatte; nuovo edge case P1-02 |
| BUG-P1-01/02 | risolto | pathway-only e Top-20 ricalcolato nel contesto |
| BUG-P1-03 | risolto | massa TPM FLOW conservata nelle integrazioni |
| BUG-P2-01/02 | risolto | `Unclassified` distinto e multi-EC completo |
| BUG-P2-03/04 | risolto | sample zero tracciato; interi CLI rigorosi |
| BUG-P2-05 | risolto per FUNZ/FLOW/PIE | audit multi-KO presente; Pathview resta parziale |
| BUG-P3-01 | risolto per target mancanti | target relativi, esistenti, non vuoti; nuovo problema cumulativo P2-02 |
| BUG-P3-02/03/04/05 | risolto | PIE round-trip, gate ID, preflight e fixture Windows verdi |
| BUG-P3-06 | migliorato ma non completo | 21 test verdi; i nuovi edge case non hanno regressioni |

### `AUDIT_SQM_PLOTS_REVIEW_2026-08-31.md`

| Finding del 31 agosto | Stato corrente | Nota |
|---|---|---|
| NEW-P1-01 BRITE accettate in `defined` | risolto | `Transporters` viene rifiutato come non-pathway |
| NEW-P1-02 taxonomy globale non dichiarata | risolto | metadati, quota esclusa e invariante 100 verificati su tutti i rank |
| NEW-P1-03 PIE dichiara Top-N non applicato | risolto | help/manifest dichiarano `all_positive_ko`, `top_n_ko=NA` |
| NEW-P1-04 file Pathview vecchi | risolto | export temporaneo isolato per invocazione |
| NEW-P2-01 Pathview senza TSV | parzialmente risolto | TSV presenti, ma non ancora esatti/fedeli: P1-01 e P2-01 |
| NEW-P2-02 palette SQMtools taxonomy | risolto per decisione di specifica | regole aggiornate con eccezione solo per `plotTaxonomy()` |
| NEW-P2-03 HTML FLOW non trasferibile | risolto | HTML self-contained e controllo Pandoc |
| NEW-P2-04 PNG per EC senza segnale | risolto | TSV di stato sì, PNG no |
| CONCEPT-P2-05 line plot | non risolto | resta default |
| NEW-P3-01 dimensioni enzimi | risolto | width/height propagate |
| NEW-P3-02 validazione input | risolto nei casi originari | restano liste vuote/duplicati: P3-01 |
| NEW-P3-03 TPM non validati | risolto | validazione numerica centralizzata prima delle aggregazioni |
| DEBT-01 monolite | non risolto | 4.604 righe e tripla ricostruzione ORF-long |

## Verifiche eseguite

### Suite e integrazioni

| Verifica | Esito |
|---|---|
| `tests/run_fast_tests.R` | PASS, 21 file |
| integrazione P0 su `in/Au_sip` | PASS; Bacillota 88.258 ORF |
| integrazione P1 | PASS; Top-20 globale e Bacillota; massa FLOW conservata |
| integrazione P2 | PASS; 3 KO multi-EC nel 00710; 705.510 ORF senza KO; 1.317 multi-KO |
| integrazione P3 | PASS; 00633/S13_1_8/K10679: KO TPM 18,284; pathway TPM 25,632 |
| taxonomy globale percent, sei rank | PASS; `displayed_percent_sum + excluded_percent = 100` |

### Copertura

Tutti i gate mirati P0-P3 superano l'80% sulle funzioni selezionate. Le baseline
globali informative restano molto inferiori:

| Harness | Copertura globale strumentata |
|---|---:|
| P0 | 10,13% |
| P1 | 10,83% |
| P2 | 18,50% |
| P3 | 16,38% |

Il verde dei gate mirati non copre l'orchestrazione completa né i nuovi casi limite.

### Probe indipendenti aggiuntivi

- confronto dei formali di `SQMtools::exportPathway()` con il TSV config Pathview;
- taxon reale a una ORF con entrambe le policy `ignore_unclassified_functions`;
- taxon reale valido con tutti i pathway definiti assenti;
- oggetto sintetico senza `misc$KEGG_names`;
- campione CLI duplicato;
- audit percentuale globale su phylum, class, order, family, genus e species.

I probe temporanei sono stati rimossi dopo l'esecuzione; nessun file di input è stato
modificato e nessun output storico è stato cancellato.

Non è stata lanciata una nuova run completa `--mode=all` né un export Pathview live:
avrebbero prodotto migliaia di file e richiesto il download KEGG senza aggiungere
evidenza ai meccanismi già dimostrati tramite API installata, exporter controllato e
quattro integrazioni reali. Questo limite va colmato nel gate finale dopo i fix.

## Affidabilità degli output

- FUNZ, FLOW, PIE e percentuali tassonomiche per pathway sono supportati nei casi
  coperti dalle integrazioni correnti.
- La taxonomy globale percentuale ora ha una semantica esplicita di quota della
  libreria totale e supera la sua invariante specifica.
- Non considerare pienamente riproducibili i Pathview correnti finché config,
  trasformazione, palette e tabella dei soli nodi plottati non coincidono.
- I filtri tassonomici comuni come Bacillota sono verificati; non generalizzare tale
  risultato ai taxa rari/non annotati finché P1-02 e P1-03 non sono corretti.
- Non usare un manifest di sezione aggiornato come prova che tutte le righe
  appartengano all'ultima esecuzione; oggi è un inventario cumulativo.
- Come già indicato dagli audit precedenti, non approvare retroattivamente output
  storici: dopo i fix usare una directory nuova per ogni run candidata.

## Priorità di correzione proposta

1. Rendere veritiero e testato il contratto Pathview: parametri espliciti, TSV esatto,
   palette e trasformazioni coerenti.
2. Correggere il filtro taxon con funzioni non classificate e gestire pathway vuoti con
   warning/skip selettivo.
3. Separare manifest della run corrente e inventario storico.
4. Rimuovere `line` dai default senza metadati ordinati.
5. Chiudere hardening CLI e fallback `KEGG_names`.
6. Refactor della tabella ORF-long e dell'orchestrazione, accompagnato da test dei
   nuovi finding e da una run reale pulita mirata prima del gate `mode=all`.

## Gate minimo dopo le correzioni

- il TSV config Pathview coincide esattamente con gli argomenti catturati dal renderer;
- ogni riga Pathview plottata è ricostruibile dal TSV e usa la palette prescritta;
- `Accipitriformes (no class in NCBI)` restituisce l'insieme ORF esatto senza errore;
- un taxon valido senza ORF nel pathway emette warning e non blocca taxonomy globale;
- il manifest della run non contiene righe di run precedenti;
- liste CLI vuote e campioni duplicati falliscono prima di `loadSQM()`;
- fallback KO valido senza `misc$KEGG_names`;
- `line` non è default senza metadati x ordinati;
- 21 test veloci, quattro integrazioni, coverage gate e nuove regressioni tutti verdi.
