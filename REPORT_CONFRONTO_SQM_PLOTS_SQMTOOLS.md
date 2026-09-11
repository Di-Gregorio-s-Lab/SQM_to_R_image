# Confronto tra `sqm_plots_lean.R`, `plotTaxonomy()` ed `exportPathway()`

## Esito

`sqm_plots_lean.R` usa due logiche distinte:

- **TAXON** delega realmente il rendering a `SQMtools::plotTaxonomy()`; la logica nativa è rispettata *dopo* che lo script ha scelto l'oggetto SQM da passargli.
- **PATHVIEW** delega realmente a `SQMtools::exportPathway()`; i valori e il mapping visualizzati seguono il KGML KEGG nativo.
- **FUNZ, FLOW e PIE** non delegano né a `plotTaxonomy()` né a `exportPathway()`: costruiscono tabelle proprie da ORF, TPM e KO. I loro totali KO sono riallineati alla matrice KO del sottoinsieme SQM, ma la definizione di “KO del pathway” non è quella KGML di `exportPathway()`.
- **`barplot_ko_pathway.R` legacy** segue il subset testuale di `KEGGPATH`, ma duplica il TPM di un ORF per ogni suo KO: non riproduce né la matrice KO SQM né la logica di Pathview.

La differenza che può cambiare materialmente i risultati è questa:

```text
FUNZ / FLOW / PIE / TAXON per pathway
  nome pathway in KEGGPATH -> ORF selezionati -> tutti i KO di quegli ORF

PATHVIEW
  ID pathway -> KGML KEGG -> KO dei nodi ortholog del KGML
```

Un ORF selezionato perché il suo campo `KEGGPATH` contiene un pathway può avere altri KO annotati. `subsetORFs()` li riaggrega tutti nella matrice KO del sottoinsieme; quindi tali KO possono apparire in FUNZ/FLOW/PIE e contribuire al TAXON per pathway, pur non appartenendo ai nodi KGML usati da PATHVIEW. L'uguaglianza fra i due output non è quindi garantita, neppure con gli stessi campioni e lo stesso ID pathway.

## Base comune e punto in cui i flussi divergono

Lo script carica SQMtools 1.7.2 e costruisce un oggetto SQM globale oppure, con `--taxa`, un contesto ottenuto con `subsetORFs(..., rescale_tpm = FALSE)`. In entrambi i casi il subset ricalcola le matrici tassonomiche e funzionali dal sottoinsieme di ORF, senza rinormalizzare i TPM.

Per ogni pathway `defined` o `top20`, `sqm_plots_lean.R` chiama:

```r
subsetFun(
  SQM = sqm, fun = pathway_name, columns = "KEGGPATH",
  ignore_case = FALSE, fixed = TRUE, allow_empty = FALSE
)
```

In SQMtools 1.7.2 `subsetFun()` applica `grepl()` al campo `KEGGPATH` degli ORF e passa gli ORF trovati a `subsetORFs()`. Quest'ultima ricostruisce `$taxa` e `$functions$KEGG` per gli ORF trattenuti. È questa matrice KO ricostruita che lo script chiama “ufficiale” nelle tabelle FUNZ/FLOW/PIE: è ufficiale **per il sottoinsieme di ORF testuale**, non per il KGML.

Al contrario, la chiamata PATHVIEW è:

```r
exportPathway(
  SQM = sqm, pathway_id = pathway_id, count = "tpm", samples = selected_samples,
  split_samples = FALSE, log_scale = FALSE, output_dir = final_dir
)
```

Qui `sqm` è il contesto globale o tassonomico, non `analysis$pathway_sqm`. `exportPathway()` parte dall'intera matrice `SQM$functions$KEGG$tpm` del contesto, scarica il KGML di `ko<pathway_id>` e lascia a `pathview::node.map()` il filtro sui nodi `ortholog`.

`subsetFun()` viene comunque eseguito prima anche per un task PATHVIEW: costruisce `prepared`, determina l'ID e prepara l'analisi comune. Non alimenta i valori esportati da `exportPathway()`, ma può impedirne l'esecuzione se il subset o l'allocazione ORF→KO falliscono; soltanto le entry presenti in `prepared` arrivano al ciclo PATHVIEW.

## Script legacy `barplot_ko_pathway.R`

Lo script legacy applica `loadSQM()` e poi lo stesso `subsetFun(...,
columns = "KEGGPATH", fixed = TRUE)` testuale usato dal lean. Sul sottoinsieme:

```text
ORF selezionato × campione × ogni KO separato da virgola
  -> conserva il TPM intero dell'ORF su ogni riga KO
  -> somma per campione × KO
  -> Top-N nel gruppo di campioni e Other
```

Non divide il TPM di un ORF multi-KO, né lo riscalda a
`functions$KEGG$tpm`. Per esempio, un ORF con TPM 10 e due KO genera 10 per
ciascun KO e un totale legacy di 20. Di conseguenza `sample_pathway_total_tpm`
nel manifesto e le percentuali derivate possono superare la massa TPM degli
ORF del pathway.

Il legacy separa `KEGG ID` soltanto sulle virgole e rimuove gli asterischi; il
lean estrae gli ID che corrispondono a `K[0-9]{5}`. Entrambi possono includere
KO non presenti nel KGML, perché selezionano gli ORF da `KEGGPATH`; soltanto
il lean conserva i totali della matrice KO del proprio subset.

Se `--samples` è assente, il legacy crea un gruppo e un PNG per ogni campione,
con Top-N calcolato separatamente. Se presente, crea un unico grafico combinato
nel loro ordine. Il lean usa invece tutti i campioni selezionati nel medesimo
barplot FUNZ. Il legacy non invoca `plotTaxonomy()` o `exportPathway()` e non
supporta contesti tassonomici, FLOW, PIE, TAXON o PATHVIEW.

Il suo unico TSV è `manifest_<pathway>.tsv`: contiene tutte le righe KO grezze
per ogni PNG/dimensione, il flag `included_in_top_n` e percentuali calcolate
sul totale legacy duplicato. Non è una tabella analitica equivalente a
`barplot_ko_data.tsv` del lean.

## Confronto per modalità

| Modalità / output | Dati effettivamente prodotti | Relazione con `plotTaxonomy()` | Relazione con `exportPathway()` | Differenze da considerare |
|---|---|---|---|---|
| `taxon`, vista globale | PNG da `SQM$taxa[[rank]][[abund|percent]]` dell'oggetto di contesto | **Diretta.** Chiama `plotTaxonomy(N=15, rescale=FALSE, ignore_unmapped=TRUE, ignore_unclassified=TRUE, no_partial_classifications=FALSE)`. | Nessuna. | Il Top-15 nativo è scelto sulla matrice completa prima della selezione di `samples`; i campioni CLI vengono applicati dopo. Con `--taxa`, la matrice è quella ricalcolata sugli ORF filtrati. |
| `taxon`, per pathway | Stesso PNG nativo, ma sul risultato di `subsetFun(... KEGGPATH ...)` | **Diretta dopo il subset.** Il comportamento di Top-N, `Other` e categorie speciali è nativo. | Non è una vista del KGML. | Il campo `KEGGPATH`, non il KGML, decide gli ORF; il grafico può quindi descrivere ORF/KO ulteriori rispetto alla mappa PATHVIEW. |
| `pathview` | File Pathview e legende, per tutti i campioni insieme o una chiamata per campione | Nessuna. | **Diretta per i valori esportati.** `count="tpm"`, `log_scale=FALSE`; i KO dello stesso nodo KGML sono sommati. | `subsetFun()` viene eseguito in preparazione, ma `exportPathway()` riceve il contesto `sqm`, non `pathway_sqm`. Un KO presente in più nodi viene copiato in ciascun nodo; non esiste il collasso `Other`. Un errore di subset/allocazione può quindi bloccare PATHVIEW prima dell'export. `separato` esegue chiamate monoscampione con `split_samples=FALSE`. |
| `funz` | TSV e barplot KO per pathway; Top-KO globale sui campioni selezionati, resto in `Other` | Nessuna. | Solo la sorgente TPM è affine: usa la matrice KO ricalcolata dal subset, non il mapping KGML né l'aggregazione KO→nodo. | I totali per KO/campione sono quelli della matrice KO del subset; l'assegnazione ORF→KO è usata per i metadati e deve poter spiegare ogni KO positivo. Top-KO e `Other` rendono il risultato volutamente diverso dalla mappa Pathview. |
| Legacy `barplot_ko_pathway.R` | PNG KO e manifest per pathway; Top-N/`Other` per gruppo di campioni | Nessuna. | Nessuna: non usa KGML o nodi Pathview. | Seleziona gli ORF con `KEGGPATH` come il lean, ma assegna il TPM intero a ogni KO multi-annotato. Diverge quindi da PATHVIEW sia nel set di KO sia nei valori; il manifesto usa totali e percentuali già potenzialmente duplicati. |
| `flow` | TSV, alluvial PNG e Sankey HTML taxon→KO; valori TPM allocati e percentuali sul TPM del pathway | Non usa `$taxa`: legge la tassonomia degli ORF e usa TPM, non `abund`/ `percent` nativi. | Non usa KGML né nodi Pathview. | L'allocazione divide inizialmente il TPM ORF tra i KO annotati e poi lo riscalda al TPM KO del subset. Top taxa/KO sono calcolati **solo sui campioni selezionati** e `Unclassified`/`Unmapped` rimangono categorie separate; è diverso dal Top-N e dai filtri di `plotTaxonomy()`. |
| `pie` | TSV e torta per campione×KO×rango; TPM tassonomico di un KO | Non usa `$taxa` né `plotTaxonomy()`. | Non usa KGML. | La distribuzione ORF→taxon viene riscalata esattamente al TPM del KO nel sottoinsieme; il Top taxa è locale al singolo campione e il resto è `Other`. Non equivale né alla tassonomia `abund/percent` nativa né a un nodo Pathview multi-KO. |
| `funz` — ENZIMI | TSV e grafici EC dall'intera matrice KO del contesto | Nessuna. | Nessuna. | Somma TPM per relazione KO→EC estratta da `misc$KEGG_names`; non è limitato al pathway e un KO con più EC contribuisce a ciascun EC. |
| `normal` | Raggruppa FUNZ, ENZIMI, FLOW, TAXON e PATHVIEW sui soli pathway `defined` | Non aggiunge logica. | Non aggiunge logica. | Eredita le differenze di ogni renderer. |
| `huge` | Come `normal`, su `defined` e `top20`; aggiunge PIE, che per default resta su `defined` | Non aggiunge logica. | Non aggiunge logica. | La selezione `top20` non è il KGML: usa `misc$KEGG_paths` e TPM KO del contesto. |
| Combinazioni elementari | Eseguono soltanto i renderer elencati | Vale la riga del renderer. | Vale la riga del renderer. | `pathview` deve essere incluso esplicitamente per avere mappe native. |
| `--plan_only` | Solo `top20.tsv` ed `errors.tsv` | Nessuna visualizzazione. | Nessuna esportazione. | È un ranking dei pathway, non una verifica di equivalenza fra i renderer. |

## Confronto quantitativo dei dati

### Selezione del pathway

`compute_top20()` somma i TPM KO per le associazioni in `sqm$misc$KEGG_paths`. Questo è un ranking basato sulla gerarchia annotata nell'oggetto SQM. Quando un pathway viene poi renderizzato, lo script cambia criterio e usa la ricerca letterale del suo nome in `ORF$KEGGPATH`. PATHVIEW usa invece i KO del KGML.

Esistono dunque tre livelli, non uno solo:

1. **ranking Top20:** KO→pathway da `misc$KEGG_paths`;
2. **FUNZ/FLOW/PIE/TAXON per pathway:** pathway→ORF da `ORF$KEGGPATH`, poi riaggregazione di tutti i KO di quegli ORF;
3. **legacy:** pathway→ORF da `ORF$KEGGPATH`, poi TPM ORF replicato per ogni KO;
4. **PATHVIEW:** pathway→nodi→KO dal KGML KEGG corrente.

Un confronto dei totali è significativo solo dopo avere verificato i quattro criteri di membership per il pathway in esame. Lo script contiene un helper per scaricare e parsare il KGML (`get_pathway_kos()`), ma il loader creato in `main_impl()` non viene consumato da `run_pipeline()`; il KGML non governa oggi FUNZ/FLOW/PIE/TAXON.

### Valori KO

Per FUNZ/FLOW/PIE, `allocate_ko_tpm()` parte dal TPM degli ORF, lo divide per il numero di KO annotati nell'ORF e, per ogni coppia KO×campione, applica un fattore che riporta la somma al TPM di `pathway_sqm$functions$KEGG$tpm`. Quando la preparazione riesce, le somme dei flussi/torte per un KO coincidono quindi con quel TPM del sottoinsieme.

Questa conservazione non implica equivalenza con PATHVIEW, perché:

- il sottoinsieme KO può essere diverso dal set KGML;
- Pathview somma più KO che condividono un nodo/reazione;
- Pathview può mostrare lo stesso KO in più nodi, senza ripartirlo;
- FUNZ/FLOW comprimono categorie in `Other`, PATHVIEW no.

Il legacy è ulteriormente distinto: non conserva la massa né al livello ORF
né al livello KO. Il suo `Other` e le sue percentuali sono coerenti con il
proprio totale duplicato, non con il TPM funzionale SQM o con la matrice dei
nodi Pathview.

### Tassonomia

Il solo renderer semanticamente equivalente a `plotTaxonomy()` è TAXON. FLOW e PIE sono analisi tassonomiche proprie: attribuiscono TPM KO agli ORF e ai loro taxa, mentre `plotTaxonomy()` disegna una matrice tassonomica già aggregata di abbondanze o percentuali. Inoltre `plotTaxonomy()` seleziona il Top-N su tutti i campioni presenti nella matrice; FLOW seleziona Top KO e Top taxa sui soli campioni richiesti e PIE lo fa per un solo campione.

Le opzioni native impostate da TAXON meritano attenzione: `rescale=FALSE` preserva i valori di SQM, e i flag `ignore_unmapped` e `ignore_unclassified` agiscono dopo una prima selezione Top-N. Se una categoria speciale non rientra nel Top-N iniziale, può essere già confluita in `Other`; i flag non sono una sottrazione globale preventiva.

## Risposta operativa

- Per una **mappa KEGG ufficiale**, l'output da usare è PATHVIEW: aderisce a `exportPathway()` e al KGML, con i parametri espliciti dello script.
- Per una **composizione tassonomica nativa**, usare TAXON: aderisce a `plotTaxonomy()` sull'oggetto SQM che lo script gli passa.
- FUNZ, FLOW e PIE sono utili per attribuire i TPM del subset agli ORF e ai taxa, ma vanno etichettati come output della logica custom di `sqm_plots_lean.R`, non come una riproduzione di PATHVIEW o `plotTaxonomy()`.
- Il barplot legacy serve soltanto come output storico: non è quantitativamente confrontabile con FUNZ o PATHVIEW finché sono presenti ORF multi-KO.
- Prima di interpretare differenze quantitative fra PATHVIEW e gli altri renderer, confrontare per ciascun pathway gli insiemi di KO da `KEGGPATH`, dalla matrice del `subsetFun()` e dal KGML. Senza questa verifica, una differenza può essere dovuta alla membership e non ai TPM.

## Evidenze consultate

- `GUIDA_LOGICA_SQM_PLOTS.md`, in particolare selezione dei pathway, allocazione ORF→KO, TAXON e PATHVIEW.
- `REPORT_PLOTTAXONOMY_SQMTOOLS.md`: lettura di `$taxa`, Top-N prima della selezione di `samples`, aggregazione `Other` e flag delle categorie speciali.
- `REPORT_EXPORTPATHWAY_SQMTOOLS.md`: matrice KEGG, download KGML, mapping ai nodi `ortholog` e somma KO→nodo.
- `sqm_plots_lean.R`: `subset_sqm_ids()` (righe 802–818), `prepare_analysis()` (820–846), `render_taxonomy()` (616–635), `export_pathview_isolated()` (637–651), `compute_top20()` (767–796) e `run_pipeline()` (1007–1147).
- `archive/legacy-scripts/barplot_ko_pathway.R`: `load_and_filter_sqm()` (181–197), `build_orf_level_tpm_table()` (221–251), `build_ko_summary_table()` (253–265), `get_top_ko_ids()` (286–294), `build_plot_table()` (296–375) e `build_manifest_rows()` (454–482).
- Namespace R locale: SQMtools 1.7.2; `subsetFun_()` e `subsetORFs_()` ispezionate per confermare la selezione `KEGGPATH` e la riaggregazione.
