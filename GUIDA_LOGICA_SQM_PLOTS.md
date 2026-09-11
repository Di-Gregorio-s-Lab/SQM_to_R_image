# Guida alla logica di `sqm_plots_lean.R`

## Scopo

Questa guida descrive il comportamento effettivo di `sqm_plots_lean.R`: input,
selezioni, trasformazioni dei TPM e file prodotti. Le convenzioni scientifiche
generali restano in `REGOLE_SCRIPT_R_SQUEEZEMETA.md`; in caso di differenza tra
questa guida e l'esecuzione, il codice è la fonte operativa.

Lo script può essere eseguito con `Rscript` oppure caricato con `source()` per
riutilizzare le funzioni. `main()` parte soltanto nell'esecuzione diretta.

## Flusso generale

```text
CLI
 ├─ validazione di opzioni e valori
 ├─ verifica dei package necessari
 ├─ loadSQM() del progetto
 ├─ validazione di ORF, campioni, ranghi e matrici TPM
 ├─ contesto globale oppure un contesto per ogni taxon richiesto
 ├─ calcolo e scrittura del Top N pathway del contesto
 ├─ risoluzione e subset dei pathway defined/top20
 ├─ preparazione della tabella ORF × campione × KO
 ├─ rendering delle modalità richieste
 └─ errors.tsv, manifest di sezione e log della run
```

Il progetto viene caricato una volta. Se `--taxa` è presente, il lavoro viene
ripetuto su subset tassonomici indipendenti senza produrre anche il contesto
globale.

## Avvio e opzioni

Gli argomenti accettano sia `--nome valore` sia `--nome=valore`.

```powershell
Rscript sqm_plots_lean.R `
  --project_dir in/Au_sip `
  --output_dir out/analisi `
  --mode normal
```

Sono obbligatori `--project_dir`, `--output_dir` e `--mode`. Le opzioni
sconosciute, i valori mancanti e le liste con valori non ammessi interrompono
la run.

| Opzione | Default | Significato |
|---|---|---|
| `mode` | obbligatorio | `huge`, `normal` oppure una lista tra `funz,flow,taxon,pie,pathview` |
| `samples` | tutti | Campioni, nell'ordine indicato |
| `tax_mode` | `prokfilter` | Modalità di `loadSQM()`: `prokfilter`, `allfilter` o `nofilter` |
| `taxa` | nessuno | Taxa da analizzare come contesti separati |
| `pathways` | nove ID curati | Pathway della selezione `defined` |
| `pathway_selection_modes` | `defined,top20` | Strategie di selezione dei pathway |
| `pathway_top_n` | `20` | Numero di pathway nel ranking |
| `top_n_ko` | `20` | KO conservati separatamente in FUNZ e FLOW |
| `top_n_taxa` | `15` | Taxa conservati separatamente in FLOW e PIE |
| `taxonomy_ranks` | `phylum,class,order,family,genus,species` | Ranghi tassonomici |
| `taxonomy_counts` | `abund,percent` | Conteggi richiesti a `plotTaxonomy()` |
| `flowplot_formats` | `png,html` | Formati FLOW; i valori ammessi sono `png` e `html` |
| `pathview_sample_modes` | `insieme,separato` | Esportazione Pathview con campioni insieme e/o singoli |
| `enzyme_ecs` | elenco curato | EC inclusi nei grafici enzimatici |
| `enzyme_plot_types` | `bar,line` | Tipi di grafico enzimatico |
| `dimensions` | `12x9,16x9,12x16` | Dimensioni dei PNG in pollici |
| `plot_dpi` | `600` | Risoluzione dei PNG |
| `workers` | min(4, core fisici) | Processi paralleli per i task pathway/modalità |
| `plan_only` | falso | Scrive Top N ed `errors.tsv`, senza grafici |
| `refresh_kegg` | falso | Forza il download del catalogo pathway KEGG |

`pathway_top_n`, `top_n_ko`, `top_n_taxa` e `workers` devono essere interi
positivi; dimensioni e DPI devono essere positivi.

I pathway definiti di default sono `00361`, `00622`, `00623`, `00621`,
`00625`, `00630`, `00633`, `00910` e `00980`. Il mapping contiene anche
`00710`, ma questo ID non fa parte del default.

## Modalità

| Modalità | Lavoro eseguito |
|---|---|
| `huge` | FUNZ, ENZIMI, FLOW, TAXON e PATHVIEW su `defined` e `top20`; PIE sui `defined` per default |
| `normal` | FUNZ, ENZIMI, FLOW, TAXON e PATHVIEW sui soli pathway `defined` |
| `funz` | Barplot KO per pathway e grafici ENZIMI |
| `flow` | Alluvial PNG e/o Sankey HTML per pathway, campione e rango |
| `taxon` | Tassonomia globale e per pathway |
| `pie` | Torta tassonomica per pathway, campione, KO e rango |
| `pathview` | Esportazione KEGG Pathview |

Le modalità elementari possono essere combinate, per esempio
`--mode funz,flow`. `huge` e `normal` devono invece essere usate da sole. Le
vecchie modalità `all` ed `enzimi` non esistono: ENZIMI è incluso in `funz`.

Per default PIE usa soltanto la selezione `defined`, anche in `huge`. PIE
include `top20` soltanto se `--pathway_selection_modes` viene specificato
esplicitamente con `top20`.

## Caricamento e validazione del progetto

Lo script chiama:

```r
SQMtools::loadSQM(
  project_path = <project_dir normalizzata>,
  tax_mode = <tax_mode>,
  trusted_functions_only = FALSE,
  load_sequences = FALSE
)
```

Richiede `sqm$orfs$table`, `sqm$orfs$tax`, `sqm$orfs$tpm` e
`sqm$functions$KEGG$tpm`. Le tre tabelle ORF devono avere row name univoci e
gli stessi insiemi di ID. I TPM devono essere numerici, finiti, non negativi e
senza valori mancanti. Campioni e ranghi richiesti devono esistere.

Se `--samples` manca, vengono usate tutte le colonne di `orfs$tpm`; altrimenti
viene mantenuto l'ordine della CLI.

## Contesti tassonomici

Ogni valore di `--taxa` viene cercato con confronto esatto, ignorando
maiuscole, minuscole e spazi esterni. Deve comparire in un solo rango: un nome
assente o presente in più ranghi è ambiguo e ferma la run.

Gli ORF corrispondenti vengono passati a `SQMtools::subsetORFs()` con TPM e
copy number non riscalati. Gli output del contesto finiscono in:

```text
<output_dir>/taxon_filter/<rango>/<taxon>/
```

## Selezione e ranking dei pathway

### `defined`

Un pathway può essere richiesto tramite ID numerico di cinque cifre o nome
KEGG completo. La risoluzione usa il catalogo `pathway/ko` di KEGG e richiede
una corrispondenza unica. Il nome viene poi passato letteralmente a:

```r
SQMtools::subsetFun(
  SQM = sqm,
  fun = pathway_name,
  columns = "KEGGPATH",
  ignore_case = FALSE,
  fixed = TRUE,
  allow_empty = FALSE
)
```

Un subset non preparabile viene registrato in `errors.tsv`; gli altri pathway
continuano.

### `top20`

Il ranking è ricalcolato in ogni contesto e usa i campioni selezionati. Parte
da `sqm$misc$KEGG_paths`, conserva soltanto gerarchie KEGG valide di tre
livelli appartenenti alle sei radici PATHWAY e scarta le categorie “not
included”. Per ogni appartenenza pathway-KO somma il TPM ufficiale del KO in
`sqm$functions$KEGG$tpm`.

Lo stesso KO può contribuire per intero a più pathway dei quali è membro. Gli
spareggi sono risolti per nome, radice e categoria. Se uno stesso nome pathway
appartiene a gerarchie diverse, la run si ferma.

Il ranking viene sempre scritto, anche con `--plan_only` e in modalità
`normal`:

```text
<context_output_dir>/top20.tsv
```

Le colonne sono `rank`, `pathway_id`, `pathway_name`, `pathway_root`,
`pathway_category`, `total_tpm`, `samples`, `filtered_taxon` e
`filtered_taxon_rank`. Le righe senza ID KEGG risolvibile restano nel TSV ma
non vengono trasformate in task `top20`.

## Allocazione ORF → KO

FUNZ, FLOW e PIE condividono `prepare_analysis()`. Per ogni ORF lo script
estrae tutti gli ID `Kxxxxx` da `KEGG ID`; per ogni campione:

1. divide il TPM dell'ORF per il numero totale dei KO annotati su quell'ORF;
2. conserva le quote dei KO presenti nel pathway;
3. per ogni coppia campione-KO riscalza le quote affinché la loro somma sia
   uguale al TPM KO ufficiale del subset in `sqm$functions$KEGG$tpm`.

In formula, per un KO con quote osservate positive:

```text
TPM allocato ORF-KO = quota ORF grezza × TPM KO ufficiale / somma quote grezze del KO
```

Questa procedura conserva quindi la massa ufficiale di ogni KO, senza
assegnare l'intero TPM ORF a ogni annotazione. Se un KO ha TPM ufficiale
positivo ma nessun ORF permette di allocarlo, la preparazione fallisce.

Le tassonomie mancanti diventano `Unclassified`. Nome funzionale ed EC vengono
da `sqm$misc$KEGG_names`; gli EC sono estratti esclusivamente dal blocco
`[EC:...]`.

## Output analitici

### FUNZ

FUNZ ordina i KO in base alla somma dei TPM ufficiali su tutti i campioni
selezionati, mantiene i primi `top_n_ko` e collassa gli altri in `Other`. La
percentuale è calcolata sul totale dei KO del pathway nello stesso campione.

```text
funz/pathway/<definiti|top20>/<pathway>/barplot_ko_data.tsv
funz/pathway/<definiti|top20>/<pathway>/barplot_ko_<dimensione>.png
```

Il TSV contiene `sample`, `ko_id`, eventuali `kegg_function` ed `ec_codes`,
`tpm`, `pathway_tpm` e `percent`.

### ENZIMI

ENZIMI usa l'intera matrice ufficiale KO del contesto, non i soli ORF di un
pathway. Collega KO ed EC tramite `misc$KEGG_names`, filtra gli EC richiesti e
somma i TPM dei KO per `sample × ec_code`. Se più KO condividono un EC, i loro
TPM vengono sommati.

```text
funz/enzimi/insieme/
funz/enzimi/separato/<EC>/
```

Ogni directory contiene `enzimi_data.tsv` e i PNG `bar` e/o `line` nelle
dimensioni richieste.

### FLOW

FLOW seleziona Top KO e Top taxa sulla somma di tutti i campioni scelti, poi
crea una tabella distinta per campione. `Unclassified` resta separato;
`Other` è riservato al collasso delle categorie fuori Top N.

Ogni arco contiene TPM e tre percentuali, tutte riferite al totale pathway del
campione:

```text
flow_percent  = 100 × TPM arco taxon-KO / TPM pathway
taxon_percent = 100 × TPM taxon / TPM pathway
ko_percent    = 100 × TPM KO / TPM pathway
```

```text
flowplot/<definiti|top20>/<pathway>/<rango>/flow_<campione>_data.tsv
flowplot/<definiti|top20>/<pathway>/<rango>/flow_<campione>_<dimensione>.png
flowplot/<definiti|top20>/<pathway>/<rango>/flow_<campione>.html
```

Il PNG è un alluvial taxon → KO. L'HTML è un Sankey autosufficiente. Entrambi
mostrano quote tassonomiche e funzionali; gli hover HTML aggiungono EC, nome
funzionale, TPM e percentuale dell'arco.

### TAXON

TAXON delega sia la vista globale sia quella per pathway a
`SQMtools::plotTaxonomy()` con `N=15`, `rescale=FALSE`,
`ignore_unmapped=TRUE`, `ignore_unclassified=TRUE` e
`no_partial_classifications=FALSE`.

```text
taxonomy_global/<abund|percent>/<rango>/taxonomy_<dimensione>.png
taxonomy_by_pathway/<definiti|top20>/<pathway>/<abund|percent>/<rango>/taxonomy_<dimensione>.png
```

Questa modalità produce soltanto PNG: non scrive TSV tassonomici e rimuove un
eventuale `taxonomy_data.tsv` preesistente nella directory di destinazione.

### PIE

Per ogni campione, KO del pathway e rango, PIE riscalza la distribuzione
tassonomica allocata al TPM KO ufficiale, conserva i primi `top_n_taxa` e
collassa gli altri in `Other`.

```text
pie/<definiti|top20>/<pathway>/<campione>/<KO>/<rango>/pie_data.tsv
pie/<definiti|top20>/<pathway>/<campione>/<KO>/<rango>/pie_<dimensione>.png
```

`pct` è una proporzione tra 0 e 1 sul totale del KO. Il TSV registra anche
`ko_sample_tpm`, `pathway_sample_tpm`, `ko_pathway_percent`, nome funzionale ed
EC. Le combinazioni senza TPM positivo non producono artefatti.

### PATHVIEW

Pathview riceve il contesto SQM, l'ID pathway, `count="tpm"`,
`split_samples=FALSE` e `log_scale=FALSE`.

```text
pathview/<definiti|top20>/<insieme|separato>/<pathway>/
```

`insieme` passa tutti i campioni in una chiamata; `separato` esegue una
chiamata per campione e aggiunge la directory del campione. Prima
dell'esportazione vengono rimossi soltanto due vecchi sidecar noti,
`pathview_input_all_ko_complete_matrix.tsv` e `pathview_render_config.tsv`.
Gli altri file della directory non vengono ripuliti.

## Cache KEGG

Il catalogo è salvato in:

```text
_cache/kegg/pathway_catalog.tsv
```

Una cache valida viene riusata; una cache corrotta viene scaricata di nuovo.
`--refresh_kegg` forza il refresh di questo catalogo.

Lo script contiene anche gli helper per cache KGML
`_cache/kegg/ko<pathway_id>.xml`, ma la pipeline corrente non li usa per
selezionare i KO: i KO del pathway arrivano dal subset SQM. Di conseguenza
`--refresh_kegg` non forza oggi un download KGML durante una run normale.

## Manifest, errori e log

I renderer interni accumulano manifest per i soli file prodotti:

```text
funz/manifest_funz.tsv
flowplot/manifest_flow.tsv
pie/manifest_pie.tsv
```

`output_file` viene convertito in path relativo alla radice della sezione. Le
righe duplicate per lo stesso file vengono eliminate. TAXON e PATHVIEW non
hanno manifest.

Gli errori dei singoli task vengono raccolti in `<output_dir>/errors.tsv` con
`task_id`, `mode`, `pathway_id` ed `error`. La pipeline continua con gli altri
task; il processo termina con codice `1` se la tabella finale contiene almeno
un errore, altrimenti con `0`.

Quando `output_dir` è ricavabile dalla CLI, l'esecuzione crea inoltre
`<output_dir>/<timestamp>_<id>.log`, con argomenti, warning, errori bloccanti e
stato finale. Gli output già presenti non vengono cancellati in blocco; i file
con lo stesso nome possono essere sovrascritti.

## Dipendenze effettive

`SQMtools` e `ggplot2` sono sempre verificati. `ggalluvial` serve quando FLOW
può essere eseguito; `plotly` e `htmlwidgets` servono per FLOW HTML.

Nell'implementazione corrente `pathview` viene richiesto per ogni run non
`plan_only`, anche quando la modalità scelta non esporta mappe. `plan_only`
richiede soltanto `SQMtools` e `ggplot2`.

## Proprietà e limiti da ricordare

- Il progetto SqueezeMeta in ingresso non viene modificato.
- Il Top N pathway è sempre calcolato e scritto.
- I task sono paralleli, ma risultati ed errori vengono riordinati per ID per
  mantenere un esito deterministico.
- La selezione pathway usa `subsetFun()`; il KGML non decide l'appartenenza dei
  KO nella pipeline corrente.
- TAXON conserva il comportamento nativo di `SQMtools::plotTaxonomy()` e non
  espone il TSV sottostante.
- La directory di output non viene pulita: una run successiva può lasciare
  artefatti storici non riscritti.
