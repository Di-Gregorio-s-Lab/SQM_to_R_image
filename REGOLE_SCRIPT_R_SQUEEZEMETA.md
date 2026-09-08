# Regole per gli script R basati su SqueezeMeta

Questa è la reference operativa per i prossimi script R di questo progetto. Le
regole derivano dal confronto tra gli script esistenti e dalle decisioni prese
sulle loro contraddizioni. In caso di conflitto, questo documento ha la
precedenza sul comportamento degli script precedenti.

## 1. Esecuzione e struttura delle directory

Gli script devono essere eseguibili da terminale con l'ambiente Conda `r_env`
già attivo:

```bash
Rscript nome_script.R \
  --project_dir in/Au_sip \
  --output_dir out/nome_script \
  [altre opzioni]
```

Regole:

- non attivare Conda dall'interno dello script;
- non usare percorsi assoluti o dipendere dalla working directory dell'autore;
- rendere obbligatori `--project_dir` e `--output_dir`;
- accettare sia `--nome valore` sia `--nome=valore`;
- implementare `--help` e terminare con stato `0` dopo averlo mostrato;
- rifiutare opzioni sconosciute e valori mancanti;
- validare enumerazioni, interi positivi, campioni, ranghi tassonomici e
  dimensioni prima di caricare o trasformare i dati;
- terminare con stato diverso da zero per input o parametri non validi;
- usare `in/Au_sip` per le prove di integrazione;
- scrivere qualsiasi risultato di prova sotto `out/`, preferibilmente in
  `out/<nome_script>/` o `out/reference_validation/`.

## 2. Importazione dell'oggetto SQM

La sorgente primaria deve essere l'oggetto restituito da `loadSQM()`. Non
leggere direttamente i TSV in `results/tables` se la stessa informazione è già
disponibile nell'oggetto.

```r
sqm <- SQMtools::loadSQM(
  project_path = project_dir,
  tax_mode = tax_mode,
  trusted_functions_only = FALSE,
  load_sequences = FALSE
)
```

Default:

- `tax_mode = "prokfilter"`;
- `trusted_functions_only = FALSE`;
- `load_sequences = FALSE`, per non caricare sequenze inutilizzate e ridurre
  il consumo di memoria.

Usare `load_sequences = TRUE` solo quando lo script elabora realmente le
sequenze. Rendere `tax_mode` configurabile da CLI se sono necessarie anche
`allfilter` o `nofilter`.

Non sopprimere gli avvisi emessi durante il caricamento. Sul dataset di prova è
atteso l'avviso relativo al progetto SqueezeMeta `1.7.3.alpha3` caricato con
SQMtools `1.7.2`: l'avviso deve restare visibile e va trattato come errore solo
se il caricamento o le verifiche strutturali falliscono.

Subito dopo l'importazione verificare almeno:

```r
required_orf_parts <- c("table", "tax", "tpm")
missing_parts <- setdiff(required_orf_parts, names(sqm$orfs))

if (length(missing_parts) > 0L) {
  stop(
    "Parti mancanti in sqm$orfs: ",
    paste(missing_parts, collapse = ", "),
    call. = FALSE
  )
}

if (is.null(rownames(sqm$orfs$table)) ||
    is.null(rownames(sqm$orfs$tax)) ||
    is.null(rownames(sqm$orfs$tpm))) {
  stop("Le tabelle ORF devono avere rownames utilizzabili come orf_id.",
       call. = FALSE)
}
```

## 3. Sorgente corretta per ciascun dato

| Informazione | Sorgente primaria |
|---|---|
| ID dell'ORF | `rownames(sqm$orfs$table)` |
| ID KO | `sqm$orfs$table[["KEGG ID"]]` |
| Nome/descrizione funzionale | `sqm$orfs$table[["KEGGFUN"]]` |
| Codice o codici EC | estrazione da `sqm$orfs$table[["KEGGFUN"]]` |
| Totali KO per campione | `sqm$functions$KEGG$tpm` |
| Appartenenza ai pathway | KO dei nodi `ortholog` del KGML KEGG |
| Tassonomia per ORF | `sqm$orfs$tax` |
| TPM per ORF e campione | `sqm$orfs$tpm` |
| Nomi KO e gerarchie per il ranking Top 20 | `sqm$misc$KEGG_names`, `sqm$misc$KEGG_paths` |

### 3.1 KO

Gli ID KO devono provenire esclusivamente dalla colonna `KEGG ID`. Non cercare
pattern `Kxxxxx` dentro `KEGGFUN`, perché quella colonna contiene la
descrizione funzionale e non è una sorgente affidabile per l'identificativo.

Rimuovere marcatori come `*`, spazi esterni e valori vuoti. Per gestire in modo
robusto uno o più KO si può estrarre il pattern `K[0-9]{5}` dalla sola colonna
`KEGG ID`.

### 3.2 Descrizione funzionale ed EC

`KEGGFUN` è la sorgente primaria sia per il nome funzionale sia per i codici
EC. I codici sono contenuti nel blocco `[EC:...]`, per esempio:

```text
(R,R)-butanediol dehydrogenase [...] [EC:1.1.1.4 1.1.1.- 1.1.1.303]
```

Regole:

- conservare tutti i codici presenti nel blocco, non soltanto il primo;
- conservare anche codici incompleti validi come `1.1.1.-`;
- rappresentare più codici in una colonna testuale con separatore `;`;
- usare `NA` quando il blocco EC non esiste;
- non estrarre codici EC da `KEGG ID`;
- usare `sqm$misc$KEGG_names` solo come fallback se `KEGGFUN` è assente o
  vuoto.

Esempio:

```r
ec_block <- stringr::str_match(
  as.character(orf_table[["KEGGFUN"]]),
  "\\[EC:([^]]+)\\]"
)[, 2]

ec_codes <- ifelse(
  is.na(ec_block),
  NA_character_,
  stringr::str_replace_all(stringr::str_squish(ec_block), "\\s+", ";")
)
```

## 4. Selezione del pathway

L'appartenenza funzionale a un pathway deriva esclusivamente dai KO presenti
nei nodi `ortholog` del relativo KGML KEGG. Estrarre ID `K[0-9]{5}`, rimuovere
i duplicati e selezionare gli ORF la cui colonna `KEGG ID` contiene almeno uno
di questi KO. I nodi compound e gli altri tipi KGML non appartengono al set KO.

`KEGGPATH` e `sqm$misc$KEGG_paths` servono soltanto a classificare, nominare e
ordinare i pathway candidati per il ranking Top 20: non devono filtrare ORF né
stabilire la membership del pathway. Gestire più pathway iterando sugli ID
canonici a cinque cifre e creando una sottodirectory distinta per ciascuno.

## 5. Tabella lunga ORF × campione × KO

Per analisi che collegano funzione, TPM e tassonomia, costruire una sola
tabella lunga centrale e derivare da essa grafici e tabelle aggregate.

Granularità finale:

```text
orf_id × sample × ko_id
```

Colonne minime:

```text
orf_id, sample, tpm, ko_id, kegg_function, ec_codes, KEGGPATH,
superkingdom, phylum, class, order, family, genus, species
```

Procedura:

1. convertire i rownames di `orfs$table`, `orfs$tax` e `orfs$tpm` in
   `orf_id`;
2. verificare che gli ID siano univoci in ogni sorgente;
3. verificare che nessun ORF richiesto sia mancante tra le tre sorgenti;
4. eseguire join espliciti tramite `orf_id`, senza affidarsi alla posizione
   delle righe;
5. trasformare `orfs$tpm` in formato lungo con `sample` e `tpm`;
6. convertire `tpm` in numerico e mantenere soltanto valori non `NA` e
   strettamente positivi;
7. estrarre i KO da `KEGG ID`;
8. dividere il TPM di ogni ORF equamente tra tutti i suoi KO;
9. espandere gli ORF multi-KO e solo dopo filtrare i KO appartenenti al KGML;
10. riscalare le quote tassonomiche di ogni KO al corrispondente margine
    ufficiale in `sqm$functions$KEGG$tpm`;
11. sostituire valori tassonomici mancanti o vuoti con `Unclassified`.

### 5.1 Regola per gli ORF multi-KO

Se un ORF ha più KO, dividere il suo TPM equamente tra tutti i KO annotati
prima di applicare il filtro di membership del pathway. La somma delle quote
prima del filtro deve coincidere con il TPM originale dell'ORF.

Per ogni coppia KO-campione, usare poi queste quote soltanto per determinare le
proporzioni tassonomiche e riscalarle affinché la loro somma coincida con il
margine ufficiale `sqm$functions$KEGG$tpm`. Un margine ufficiale positivo
senza ORF allocabili è un errore, non un valore da stimare o ignorare.

Gli ORF senza KO vengono esclusi dalle analisi KO e dai flowplot, ma possono
restare in analisi esclusivamente tassonomiche se coerente con l'obiettivo
dello script. Il numero di ORF esclusi deve essere riportato.

## 6. Aggregazioni, Top N e percentuali

### 6.1 Top N

Quando più campioni sono mostrati nello stesso grafico:

1. sommare il TPM di ogni categoria su tutti i campioni mostrati;
2. ordinare per TPM decrescente;
3. usare l'ID o il nome della categoria come spareggio deterministico;
4. mantenere lo stesso Top N in ogni campione;
5. aggregare tutte le categorie restanti in `Other`.

Per un grafico relativo a un solo campione, la stessa procedura si applica a
quel campione. `Unclassified` è una categoria informativa distinta e non deve
essere rinominata automaticamente in `Other`, salvo l'eccezione nativa e
tracciata di `taxonomy_global` descritta nella sezione 6.2.

### 6.2 Denominatori

Ogni colonna percentuale deve avere un nome o una descrizione che identifichi
il denominatore.

| Misura | Formula |
|---|---|
| quota KO nel pathway | `KO_TPM_ufficiale / somma_KO_TPM_ufficiali_del_pathway` nello stesso campione |
| composizione tassonomica di un KO | `taxon_KO_TPM_riscalato / KO_TPM_ufficiale` nello stesso campione |
| quota di un arco taxon → KO | `edge_TPM_riscalato / somma_KO_TPM_ufficiali_del_pathway` nello stesso campione |
| composizione tassonomica del pathway | `taxon_TPM_riscalato / somma_KO_TPM_ufficiali_del_pathway` nello stesso campione |

Non dividere mai il valore di un singolo campione per il totale sommato su
tutti i campioni, salvo che l'analisi richieda esplicitamente una percentuale
globale e la etichetti come tale.

Dopo l'aggregazione verificare, con una tolleranza numerica, che le percentuali
attese sommino a 100 per ciascun gruppo pertinente.

Eccezione esplicita: `taxonomy_global/percent` prodotto da
`SQMtools::plotTaxonomy()` mantiene il comportamento nativo non riscalato con
`ignore_unmapped=TRUE` e `ignore_unclassified=TRUE`, senza riscrivere il plot o
la selezione Top N di SQMtools. In SQMtools 1.7.2 una categoria richiesta come
esclusa può restare nei dati mostrati, eventualmente dentro `Other`, se non è
selezionata come riga distinta prima del filtro. Inoltre, nei contesti filtrati
senza riscalatura la matrice percentuale può sommare a meno di 100.

Il TSV deve quindi distinguere categorie richieste, effettivamente escluse e
trattenute; deve riportare `raw_percent_sum`, `displayed_percent_sum`,
`excluded_percent`, `accounted_percent_sum`, `denominator_type` e
`denominator_value`. La verifica pertinente è
`displayed_percent_sum + excluded_percent = raw_percent_sum` per campione.
Quando una categoria positiva richiesta come esclusa resta nei dati mostrati,
lo script emette un warning e prosegue soltanto se la massa è riconciliata in
modo univoco. Differenze inspiegabili o ambigue restano errori bloccanti.

## 7. Output e tracciabilità

### 7.1 Directory
Creare la directory di output con:

```r
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
```

Regole:

- non cancellare ricorsivamente una directory di output già esistente;
- assegnare a ogni corsa un `run_id` nel formato
  `YYYYMMDDTHHMMSS_UTCpHHMM_<4hex>` (`UTCmHHMM` per offset negativi);
- usare nomi deterministici senza `run_id` per gli artefatti: un target
  rigenerato viene sovrascritto, mentre gli altri file restano invariati;
- separare pathway, ranghi e campioni in sottodirectory con nomi sanitizzati;
- usare TSV per le tabelle di supporto;
- ogni grafico deve avere una tabella TSV contenente esattamente i dati
  utilizzati per costruirlo;
- creare manifest locali soltanto per FLOW, FUNZ e PIE. Ogni manifest descrive
  soltanto gli artefatti completati dalla corsa corrente nella propria sezione.

Eccezione di provenienza Pathview: il file
`pathview_input_all_ko_complete_matrix.tsv` documenta la matrice KO
completa passata tramite l'oggetto SQM. È un superset di input e non deve
essere descritto come tabella dei soli nodi effettivamente disegnati da
Pathview.

Il manifest deve contenere almeno:

```text
run_id, script, project_dir, tax_mode, pathway, samples, metric,
top_n_taxa, top_n_ko, output_type, output_file
```

Aggiungere, quando pertinenti, dimensioni, DPI, rango tassonomico e modalità
di raggruppamento. Preferire percorsi relativi alla directory di output per
rendere il manifest trasferibile.

Pathview e gli output prodotti tramite `plotTaxonomy()` non hanno manifest. Non
si creano `manifest_all`, manifest di corsa o manifest di fallimento. Ogni
esecuzione scrive invece `output_dir/<run_id>.log`, con stato, inizio/fine,
durata, argomenti CLI, versioni, campioni, warning, pathway saltati ed eventuale
errore. Una corsa fallita conserva gli artefatti parziali e termina con codice
non zero.

### 7.2 Colori

usa sempre questa palette di colori:

colori_hex <- c(
  "#5d8aa8", "#e32636", "#efdecd", "#ffbf00",
  "#9966cc", "#a4c639", "#cd9575", "#915c83",
  "#008000", "#fbceb1", "#00ffff",
  "#4b5320", "#b2beb5", "#87a96b", "#ff9966", "#a52a2a",
  "#6e7f80", "#ff2052", "#007fff", "#f0ffff", "#89cff0",
  "#f4c2c2", "#21abcd", "#fae7b5", "#ffe135", "#848482",
  "#98777b", "#f5f5dc", "#3d2b1f",
  "#fe6f5e", "#000000", "#ffebcd", "#318ce7", "#ace5ee", "#faf0be",
  "#0000ff", "#a2a2d0", "#6699cc", "#0d98ba", "#8a2be2", "#8a2be2",
  "#de5d83", "#79443b", "#0095b6", "#e3dac9", "#cc0000", "#006a4e",
  "#873260", "#0070ff", "#b5a642", "#cb4154", "#1dacd6", "#66ff00",
  "#bf94e4", "#c32148", "#ff007f", "#08e8de", "#d19fe8", "#f4bbff",
  "#ff55a3", "#fb607f", "#004225", "#cd7f32", "#a52a2a", "#ffc1cc",
  "#e7feff", "#f0dc82"
)

Eccezioni esplicite: i grafici creati da `SQMtools::plotTaxonomy()` conservano
la palette nativa di SQMtools. Non passare a `plotTaxonomy()` l'intero vettore
`colori_hex`, perché la funzione richiede un numero di colori coerente con
`N` e altrimenti lo ignora emettendo un warning. Anche Pathview conserva i
colori nativi: non passare `sample_colors` e registrare
`color_source=pathview_native`. Queste eccezioni non si applicano agli altri
grafici dello script.

### 7.3 Ordine dei campioni nei line plot

Se `--samples` è presente, usare esattamente l'ordine CLI; altrimenti
conservare l'ordine delle colonne SQM. I line plot enzimatici collegano i
campioni secondo questo ordine di visualizzazione e non dichiarano, da soli,
una semantica temporale o una continuità sperimentale.

## 8. Gestione degli errori

Usare `stop(..., call. = FALSE)` per:

- directory del progetto inesistente;
- struttura SQM incompleta;
- campioni o ranghi richiesti assenti;
- pathway inesistente o ambiguo;
- matrici sorgente vuote o chiavi ORF duplicate/mancanti;
- opzioni CLI non valide.

Consentire una membership KGML vuota. Usare `warning(..., call. = FALSE)` e
saltare soltanto la combinazione interessata quando un pathway valido non
contiene KO `ortholog` allocabili o righe positive per uno specifico campione,
rango o grafico. Un pathway vuoto non deve bloccare la tassonomia globale o
altre modalità indipendenti.

I messaggi finali devono indicare chiaramente la directory di output e i
manifest creati.

## 9. Anti-pattern individuati negli script esistenti

| Anti-pattern | Perché evitarlo | Regola sostitutiva |
|---|---|---|
| Usare `orfs$abund` e chiamare i valori TPM | `abund` contiene abbondanze in reads, non TPM | usare `orfs$tpm` |
| Estrarre `Kxxxxx` da `KEGGFUN` | `KEGGFUN` contiene la descrizione | usare `KEGG ID` |
| Cercare gli EC in `KEGG ID` | l'ID KO non contiene l'annotazione EC | estrarre `[EC:...]` da `KEGGFUN` |
| Leggere direttamente `*.orf.tax.*.tsv` e `*.KO.names.tsv` | duplica logica già offerta da `loadSQM()` | usare `orfs$tax`, `KEGGFUN` e `misc` |
| Usare `KEGGPATH` o `subsetFun()` per stabilire la membership | diverge dai nodi realmente disegnati da KEGG | usare soltanto i KO dei nodi `ortholog` KGML |
| Assegnare a ogni KO l'intero TPM di un ORF multi-KO | gonfia le proporzioni grezze | dividere prima tra tutti i KO, poi filtrare il pathway |
| Usare la somma ORF come totale funzionale KO | può divergere dall'oracolo SQMtools | usare `sqm$functions$KEGG$tpm` e riscalare le quote tassonomiche al suo margine |
| Calcolare la contribuzione di un campione sul totale di tutti i campioni | il denominatore non rappresenta il campione | usare il totale del pathway nello stesso campione |
| Selezionare Top N diversi nello stesso grafico multicampione | colori e categorie non sono confrontabili | graduatoria globale sui campioni mostrati |
| Usare insieme `fixed=TRUE` e `ignore_case=TRUE` | R ignora `ignore.case` in questa combinazione | risolvere prima il nome canonico |
| Cancellare tutto il contenuto di `output_dir` | può rimuovere risultati non appartenenti allo script | creare sottodirectory e sovrascrivere solo file noti |
| Usare path assoluti o oggetti preesistenti nella sessione R | lo script non è riproducibile da terminale | CLI completa e `loadSQM()` nello script |

## 10. Checklist prima di considerare uno script completato

- [ ] Lo script mostra un help utilizzabile da terminale.
- [ ] Input e output sono passati tramite CLI.
- [ ] L'importazione usa `loadSQM()` con opzioni esplicite.
- [ ] KO, funzione ed EC provengono dalle colonne corrette.
- [ ] La membership del pathway usa soltanto i KO dei nodi `ortholog` KGML.
- [ ] `KEGGPATH` è usato soltanto per classificare e ordinare i Top 20.
- [ ] I join usano `orf_id` e verificano le chiavi.
- [ ] Il valore quantitativo è realmente TPM.
- [ ] Il TPM multi-KO è diviso prima del filtro pathway.
- [ ] I totali KO coincidono con `sqm$functions$KEGG$tpm`.
- [ ] Top N e denominatori sono coerenti tra campioni.
- [ ] `Unclassified` e `Other` mantengono significati distinti.
- [ ] Ogni grafico ha TSV sorgente e manifest.
- [ ] Nessun output estraneo viene cancellato.
- [ ] La prova usa `in/Au_sip` e scrive soltanto sotto `out/`.

## 11. Baseline del dataset di prova

Con `in/Au_sip` e l'ambiente `r_env`, la validazione iniziale ha osservato:

- `945388` ORF nell'oggetto completo;
- `3` campioni: `S13_1_8`, `S13_2_8`, `S13_3_8`;
- presenza di `orfs$table`, `orfs$tax` e `orfs$tpm`;
- membership del pathway ricavata dai soli nodi `ortholog` KGML, senza usare
  `KEGGPATH` come filtro;
- coincidenza, entro la tolleranza numerica, tra i totali prodotti e i margini
  KO ufficiali `sqm$functions$KEGG$tpm`;
- conservazione del TPM grezzo dividendo gli ORF multi-KO prima del filtro e
  conservazione del margine ufficiale dopo il riscalamento tassonomico;
- presenza di funzioni con più EC, per esempio
  `[EC:1.1.1.4 1.1.1.- 1.1.1.303]`.

Queste osservazioni, validate dal motore canonico il 7 settembre 2026, servono
come smoke test del dataset attuale e non come costanti da inserire negli
script di analisi.
