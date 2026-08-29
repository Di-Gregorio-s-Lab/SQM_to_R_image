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
| Appartenenza ai pathway | `sqm$orfs$table[["KEGGPATH"]]` |
| Tassonomia per ORF | `sqm$orfs$tax` |
| TPM per ORF e campione | `sqm$orfs$tpm` |
| Nomi e gerarchie KEGG di fallback | `sqm$misc$KEGG_names`, `sqm$misc$KEGG_paths` |

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

Usare `subsetFun()` dopo `loadSQM()`, filtrando `KEGGPATH`.

Il confronto richiesto è letterale e case-insensitive. In R,
`grepl(..., fixed = TRUE)` ignora `ignore.case = TRUE`; pertanto non passare
direttamente questa combinazione a `subsetFun()`.

Procedura:

1. estrarre i nomi canonici dei pathway da `sqm$misc$KEGG_paths`;
2. confrontare `tolower(nome_canonico)` con `tolower(nome_richiesto)`;
3. richiedere una sola corrispondenza esatta; con zero o più corrispondenze
   terminare mostrando i candidati utili;
4. passare il nome canonico a:

```r
pathway_sqm <- SQMtools::subsetFun(
  SQM = sqm,
  fun = canonical_pathway_name,
  columns = "KEGGPATH",
  ignore_case = FALSE,
  fixed = TRUE
)
```

Non usare il codice numerico del pathway come filtro se `KEGGPATH` contiene
solo la descrizione testuale. Gestire più pathway iterando su nomi canonici e
creando una sottodirectory distinta per ciascuno.

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
8. espandere gli ORF multi-KO su più righe;
9. sostituire valori tassonomici mancanti o vuoti con `Unclassified`.

### 5.1 Regola per gli ORF multi-KO

Se un ORF ha più KO, ciascun KO riceve l'intero TPM dell'ORF. Non dividere il
TPM per il numero di KO.

Questa scelta può aumentare il totale dopo l'espansione. Di conseguenza:

- il comportamento va dichiarato nei commenti dello script e nei metadati;
- i denominatori delle percentuali KO devono essere calcolati sulla tabella
  espansa, così le categorie mostrate sommano coerentemente a 100%;
- non descrivere il metodo come “TPM diviso equamente tra i KO”.

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
essere rinominata automaticamente in `Other`.

### 6.2 Denominatori

Ogni colonna percentuale deve avere un nome o una descrizione che identifichi
il denominatore.

| Misura | Formula |
|---|---|
| quota KO nel pathway | `KO_TPM / pathway_TPM` nello stesso campione |
| composizione tassonomica di un KO | `taxon_KO_TPM / KO_TPM` nello stesso campione |
| quota di un arco taxon → KO | `edge_TPM / pathway_TPM` nello stesso campione |
| composizione tassonomica del pathway | `taxon_TPM / pathway_TPM` nello stesso campione |

Non dividere mai il valore di un singolo campione per il totale sommato su
tutti i campioni, salvo che l'analisi richieda esplicitamente una percentuale
globale e la etichetti come tale.

Dopo l'aggregazione verificare, con una tolleranza numerica, che le percentuali
attese sommino a 100 per ciascun gruppo pertinente.

## 7. Output e tracciabilità

### 7.1 Directory
Creare la directory di output con:

```r
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
```

Regole:

- non cancellare ricorsivamente una directory di output già esistente;
- sovrascrivere soltanto file omonimi che lo script dichiara di produrre;
- separare pathway, ranghi e campioni in sottodirectory con nomi sanitizzati;
- usare TSV per le tabelle di supporto;
- ogni grafico deve avere una tabella TSV contenente esattamente i dati
  utilizzati per costruirlo;
- ogni esecuzione deve produrre un manifest TSV.

Il manifest deve contenere almeno:

```text
script, project_dir, tax_mode, pathway, samples, metric,
top_n_taxa, top_n_ko, output_type, output_file
```

Aggiungere, quando pertinenti, dimensioni, DPI, rango tassonomico e modalità
di raggruppamento. Preferire percorsi relativi alla directory di output per
rendere il manifest trasferibile.

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

## 8. Gestione degli errori

Usare `stop(..., call. = FALSE)` per:

- directory del progetto inesistente;
- struttura SQM incompleta;
- campioni o ranghi richiesti assenti;
- pathway inesistente o ambiguo;
- matrici vuote o chiavi ORF duplicate/mancanti;
- opzioni CLI non valide.

Usare `warning(..., call. = FALSE)` e saltare soltanto la combinazione
interessata quando un pathway valido non contiene righe positive per uno
specifico campione, rango o grafico.

I messaggi finali devono indicare chiaramente la directory di output e i
manifest creati.

## 9. Anti-pattern individuati negli script esistenti

| Anti-pattern | Perché evitarlo | Regola sostitutiva |
|---|---|---|
| Usare `orfs$abund` e chiamare i valori TPM | `abund` contiene abbondanze in reads, non TPM | usare `orfs$tpm` |
| Estrarre `Kxxxxx` da `KEGGFUN` | `KEGGFUN` contiene la descrizione | usare `KEGG ID` |
| Cercare gli EC in `KEGG ID` | l'ID KO non contiene l'annotazione EC | estrarre `[EC:...]` da `KEGGFUN` |
| Leggere direttamente `*.orf.tax.*.tsv` e `*.KO.names.tsv` | duplica logica già offerta da `loadSQM()` | usare `orfs$tax`, `KEGGFUN` e `misc` |
| Filtrare manualmente `KEGGPATH` in alcuni script e con `subsetFun()` in altri | produce sottoinsiemi non uniformi | risolvere il nome canonico e usare `subsetFun()` |
| Dichiarare che i multi-KO sono divisi, ma duplicare il TPM | documentazione e risultati divergono | assegnare esplicitamente il TPM intero a ciascun KO |
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
- [ ] Il pathway è selezionato tramite nome canonico e `subsetFun()`.
- [ ] I join usano `orf_id` e verificano le chiavi.
- [ ] Il valore quantitativo è realmente TPM.
- [ ] Il comportamento multi-KO è TPM intero per ogni KO.
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
- `473` ORF dopo il filtro canonico
  `Chlorocyclohexane and chlorobenzene degradation`;
- `22` KO nel sottoinsieme;
- coincidenza, entro la tolleranza numerica, tra la somma manuale dei TPM per
  KO e `pathway_sqm$functions$KEGG$tpm` per questo sottoinsieme;
- nessun ORF multi-KO nel pathway di prova, quindi la regola multi-KO deve
  essere controllata anche con un caso minimo costruito in memoria;
- presenza di funzioni con più EC, per esempio
  `[EC:1.1.1.4 1.1.1.- 1.1.1.303]`.

Questi numeri servono come smoke test del dataset attuale, non come costanti da
inserire negli script di analisi.
