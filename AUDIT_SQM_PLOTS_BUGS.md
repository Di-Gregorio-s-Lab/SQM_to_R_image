# Audit bug `sqm_plots.R`

> Stato del documento: handoff operativo, audit del 2026-08-30. Tutti i bug elencati sono ancora aperti nello snapshot corrente, salvo diversa annotazione. I numeri di linea si riferiscono allo script di 2866 righe presente al momento dell'audit e potranno spostarsi durante le correzioni: usare anche il nome della funzione indicato.

Questo documento guida la prossima sessione di correzione di `sqm_plots.R`. Non e' un semplice elenco: per ogni problema riporta il riferimento al codice, la prova osservata, l'effetto sui risultati, la direzione di correzione e il test che deve impedirne la ricomparsa.

## Leggere prima di usare gli output esistenti

Gli output gia' prodotti in `out/` che dipendono dai punti sotto **non devono essere approvati retroattivamente**. In particolare, fino alla correzione e rigenerazione non sono affidabili:

- gli output selezionati per ID numerico per i tre pathway con mapping errato;
- tutti i grafici e TSV percentuali in `taxonomy_by_pathway/`;
- tutti gli output sotto `taxon_filter/`, indipendentemente dalla modalita';
- le selezioni e gli output `top20`, soprattutto quando combinati con `--taxa`;
- i TSV e diagrammi FLOW quando un KO ha piu' descrizioni;
- i raggruppamenti `Other`/`Unclassified`, i metadati EC e i manifest storici.

Dopo i fix usare una **nuova directory di output**, non fondere una run corretta con `out/smoke_all` o con altri risultati storici. Conservare questi ultimi solo come prova del difetto. Non cancellarli e non considerarli validati solo perche' i test correnti passano.

## Stato verificato

- Script auditato: `sqm_plots.R`, 2866 righe.
- Input reale disponibile: `in/Au_sip`.
- Baseline osservata: 945388 ORF, sample `S13_1_8`, `S13_2_8`, `S13_3_8`.
- Pathway di prova osservato: 473 ORF, 22 KO.
- Ambiente R usato nelle verifiche precedenti: `C:\Progra~1\R\R-4.5.0\bin\Rscript.exe`.
- SQMtools osservato: 1.7.2; progetto SqueezeMeta: 1.7.3.alpha3.
- Il warning di compatibilita' tra SQMtools 1.7.2 e il progetto generato con SqueezeMeta 1.7.3.alpha3 e' atteso nello snapshot corrente: non sopprimerlo senza verificarne la causa.
- I test esistenti passano, ma non coprono i bug principali; almeno un test codifica una mappatura KEGG errata.

Fonti di verita' da usare durante la correzione:

- regole analitiche locali: `REGOLE_SCRIPT_R_SQUEEZEMETA.md`;
- mapping KEGG ufficiali: <https://www.kegg.jp/pathway/map00710>, <https://www.kegg.jp/pathway/map00633> e <https://www.kegg.jp/pathway/map00910>;
- implementazioni installate di SQMtools 1.7.2 per verificare la semantica di `plotTaxonomy()`, `subsetTax()` e `subsetORFs()`; non dedurla soltanto dal nome delle funzioni.

## Indice rapido dei problemi aperti

| ID | Priorita' | Severita' | Risultato a rischio |
|---|---:|---|---|
| BUG-P0-01 | P0 | critica | identita' biologica di tre pathway e Pathview |
| BUG-P0-02 | P0 | critica | percentuali tassonomiche per pathway |
| BUG-P0-03 | P0 | critica | qualsiasi output filtrato per taxon |
| BUG-P1-01 | P1 | alta | significato della selezione `top20` |
| BUG-P1-02 | P1 | alta | `top20` dichiarato specifico per taxon |
| BUG-P1-03 | P1 | alta | TPM e archi nei FLOW |
| BUG-P2-01 | P2 | media-alta | distinzione biologica `Unclassified`/`Other` |
| BUG-P2-02 | P2 | media | completezza dei metadati multi-EC |
| BUG-P2-03 | P2 | media | esecuzione FUNZ con sample vuoti |
| BUG-P2-04 | P2 | media | validazione degli interi CLI |
| BUG-P2-05 | P2 | media | tracciabilita' multi-KO e ORF esclusi |
| BUG-P3-01 | P3 | media | attendibilita' dei manifest |
| BUG-P3-02 | P3 | media | riproducibilita' dei PIE dai TSV |
| BUG-P3-03 | P3 | media | robustezza Pathview con ID assente |
| BUG-P3-04 | P3 | media-bassa | diagnostica delle dipendenze R |
| BUG-P3-05 | P3 | bassa | igiene e portabilita' del test Windows |
| BUG-P3-06 | P3 | media | capacita' della suite di rilevare falsi risultati |

## P0 - Correggere prima

### BUG-P0-01: mappature KEGG errate nei pathway curati

Severita': critica. Stato: corretto e verificato su `fix/p0-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:59-69`: `known_pathways`.
- `sqm_plots.R:163-170`: help CLI stampa gli stessi codici.
- `sqm_plots.R:481-528`: `resolve_pathways()` usa `known_pathways` per accettare codici numerici e assegnare `pathway_id`.
- `sqm_plots.R:2095-2169`: `run_pathview_mode()` passa `pathway_id` a `SQMtools::exportPathway()`.
- `tests/test_pathview_sample_modes.R:65-66`: test storico codifica `Nitrogen metabolism = 00643`.

Prova:

- KEGG REST ufficiale associa:
  - `map00710` a `Carbon fixation in photosynthetic organisms`, non `00622`.
  - `map00633` a `Nitrotoluene degradation`, non `00642`.
  - `map00910` a `Nitrogen metabolism`, non `00643`.
- I codici presenti nello script puntano invece a pathway diversi:
  - `00622` = Xylene degradation.
  - `00642` = Ethylbenzene degradation.
  - `00643` = Styrene degradation.
- Dimensioni reali interessate gia' osservate nel dataset: 1776 ORF, 135 ORF, 1786 ORF.

Impatto:

- Selezione CLI numerica, Pathview e manifest dichiarano pathway biologicamente diversi da quelli elaborati.
- Il problema e' silenzioso: i nomi possono sembrare corretti mentre l'ID esportato e' sbagliato.

Direzione di correzione:

- Riallineare `known_pathways` agli ID KEGG ufficiali.
- Aggiornare i test che assumono `00643` per `Nitrogen metabolism`.
- Verificare separatamente la risoluzione degli ID per pathway curati e per nomi completi non presenti nella tabella curata; un valore `NA` non deve mai raggiungere Pathview.

Test di regressione:

- Test unitario su `known_pathways` per i tre mapping sopra.
- Test su `resolve_pathways()` con codice numerico e nome canonico.
- Test Pathview che confermi `Nitrogen metabolism -> 00910`.

### BUG-P0-02: percentuali tassonomiche per pathway con denominatore sbagliato

Severita': critica. Stato: corretto e verificato su `fix/p0-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:1353-1367`: `build_taxonomy_plot()` crea `sqm_subset` con `subsetSamples()` e chiama `SQMtools::plotTaxonomy()`.
- `sqm_plots.R:1373-1389`: `extract_taxonomy_plot_data()` prende `plot_object$data` senza ricalcolare il denominatore.
- `sqm_plots.R:1926-2073`: `run_taxonomy_mode()` usa la stessa funzione sia per tassonomia globale sia per `taxonomy_by_pathway`.
- `sqm_plots.R:2738-2766`: nel ramo pathway viene passato `pathway_info$pathway_sqm`, ma il dato percentuale viene comunque generato da `plotTaxonomy()`.
- `REGOLE_SCRIPT_R_SQUEEZEMETA.md:240-250`: le percentuali devono dichiarare e rispettare il denominatore; la composizione tassonomica del pathway deve essere sul pathway nello stesso campione.

Prova:

- File storico: `out/smoke_all/taxonomy_by_pathway/Chlorocyclohexane_and_chlorobenzene_degradation/percent/phylum/taxonomy_Chlorocyclohexane_and_chlorobenzene_degradation_percent_phylum_data.tsv`.
- Somme percentuali osservate per sample: `0.01916393`, `0.00667423`, `0.03347328`, non 100.

Impatto:

- I grafici `taxonomy_by_pathway/.../percent/...` non mostrano la composizione del pathway, ma una quota rispetto a un totale piu' ampio.
- I TSV sorgente e i PNG derivati possono sembrare validi ma raccontano un denominatore diverso da quello atteso.

Direzione di correzione:

- Per lo scope `taxonomy_by_pathway`, costruire il TSV a partire dagli ORF del pathway e calcolare `taxon_TPM / pathway_TPM` per sample.
- Mantenere `plotTaxonomy()` solo dove il suo denominatore coincide con il significato richiesto, oppure annotare esplicitamente il denominatore se resta globale.

Test di regressione:

- Per ogni sample con denominatore pathway positivo in `taxonomy_by_pathway` e `count=percent`, la somma deve essere 100 entro tolleranza.
- Un sample con denominatore zero deve essere marcato senza dati o saltato con warning tracciabile; non deve essere forzato a 100 e non deve interrompere tutta la run.
- Il TSV deve esporre o documentare il denominatore usato.

### BUG-P0-03: filtro taxon risolto sulle ORF ma applicato tramite contig

Severita': critica. Stato: corretto e verificato su `fix/p0-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:688-719`: `resolve_taxa_filters()` cerca il taxon in `sqm$orfs$tax`.
- `sqm_plots.R:722-732`: `subset_sqm_by_taxon()` chiama `SQMtools::subsetTax()`.
- `sqm_plots.R:2581-2588`: `run_pipeline()` risolve i taxa e poi applica `subset_sqm_by_taxon()`.

Percorso reale:

1. `main()` legge `--taxa` in `requested_taxa` (`sqm_plots.R:2511-2512`).
2. `run_pipeline()` chiama `resolve_taxa_filters(sqm, requested_taxa)` (`sqm_plots.R:2581-2582`).
3. `resolve_taxa_filters()` cerca il rank in `sqm$orfs$tax` (`sqm_plots.R:688-696`).
4. `run_pipeline()` chiama `subset_sqm_by_taxon()` (`sqm_plots.R:2583-2585`).
5. `subset_sqm_by_taxon()` delega a `SQMtools::subsetTax()` (`sqm_plots.R:722-732`), che filtra secondo la tassonomia dei contig.

Prova Bacillota gia' osservata:

- Match ORF attesi: 88258.
- ORF restituite: 92735.
- ORF non-Bacillota incluse: 5012.
- ORF Bacillota escluse: 535.

Impatto:

- Tutti gli output sotto `taxon_filter/...` possono includere ORF di taxa non richiesti ed escludere ORF richieste.
- Il bug si propaga a FUNZ, FLOW, TAXON, PATHVIEW e PIE per contesti filtrati.

Direzione di correzione:

- Se il filtro e' definito su `sqm$orfs$tax`, il subset deve selezionare esattamente gli `orf_id` corrispondenti.
- La direzione gia' identificata e' usare `SQMtools::subsetORFs(..., tax_source = "orfs", rescale_tpm = FALSE, rescale_copy_number = FALSE)`.

Test di regressione:

- Su `Bacillota`, gli ID ORF nel subset devono coincidere esattamente con gli ID ORF che hanno `phylum == "Bacillota"` in `sqm$orfs$tax`.
- Nessuna ORF fuori taxon inclusa; nessuna ORF del taxon esclusa.

## Evidenza di chiusura P0 - 2026-08-30

Checkpoint GREEN:

- `827922a`: mapping KEGG ufficiali e gate Pathview con ID valido;
- `8c71f94`: filtro tassonomico applicato agli ID ORF esatti;
- `9e96d1a`: percentuali tassonomiche pathway/sample calcolate dai TPM ORF.

Verifiche osservate:

- `tests/run_fast_tests.R`: 7 file di test completati, exit `0`;
- `tests/check_p0_coverage.R`: gate `covr` 3.6.5 superato; copertura funzioni P0 tra 84.21% e 100%, baseline informativa dell'intero script 10.71%;
- `tests/test_p0_integration_Au_sip.R`: exit `0`, 88258 ORF Bacillota restituite senza inclusioni o esclusioni, tre sample del pathway `00361` con somma percentuale 100;
- run CLI mirata in `out/p0_candidate_20260830_01`: exit `0`; il TSV percentuale pathway contiene denominatore, stato e flag `plotted`, con somma 100 per `S13_1_8`, `S13_2_8` e `S13_3_8`;
- review indipendente: nessun rilievo CRITICAL/HIGH.

Limitazione non-P0 riconfermata dalla run candidata: sul lungo percorso Windows del workspace il nome del PNG pathway e' stato troncato e il target previsto nel manifest non esiste. Il difetto resta sotto BUG-P3-01/portabilita' degli output; non altera il TSV e le invarianti scientifiche usate per chiudere P0, ma impedisce di considerare l'intera directory candidata un output finale approvato.

## P1 - Corretto su `fix/p1-correctness`

### BUG-P1-01: `top20` include BRITE e categorie non-pathway

Severita': alta. Stato: corretto e verificato su `fix/p1-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:74-81`: radici KEGG PATHWAY ammesse.
- `sqm_plots.R:547-648`: parsing gerarchico e membership ORF/pathway deduplicata.
- `sqm_plots.R:654-742`: ranking pathway-only con metadati gerarchici e spareggio deterministico.

Prova:

- Nel top 20 reale gia' osservato, 11 elementi non sono map pathway ma BRITE/categorie funzionali: `Transporters`, `DNA repair and recombination proteins`, `Enzymes with EC numbers`, `Function unknown`, `Peptidases`, `Transcription factors`, ecc.

Impatto:

- La modalita' `top20` non significa "top 20 pathway KEGG".
- Pathview viene disabilitato per `top20` quando manca `pathway_id`, ma FUNZ/FLOW/TAXON/PIE possono comunque lavorare su insiemi non-pathway.

Direzione di correzione:

- Conservare e usare la radice gerarchica del campo `KEGGPATH`, distinguendo vere pathway map da BRITE/categorie non-pathway.
- Il ranking `top20` deve contenere solo pathway biologici ammessi dalla semantica del modo.

Test di regressione:

- Dataset sintetico con una vera pathway e una voce BRITE ad alto TPM: la BRITE non deve entrare nel `top20`.
- Test sul dataset reale: nessuno dei nomi non-pathway elencati deve comparire tra i pathway selezionati.

### BUG-P1-02: `top20` specifico per taxon e' in realta' globale

Severita': alta. Stato: corretto e verificato su `fix/p1-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:744-771`: `select_context_pathway_groups()` mantiene `defined` e calcola `top20` sullo SQM ricevuto.
- `sqm_plots.R:3057-3064`: `defined` viene risolto una sola volta sul progetto completo.
- `sqm_plots.R:3112-3148`: gruppi, entry e subset pathway vengono costruiti dentro ogni `filter_context`.

Impatto:

- Quando l'utente usa `--taxa` insieme a `--pathway_selection_modes=top20`, i pathway scelti sono i top del progetto intero, non del taxon filtrato.
- Gli output sotto `taxon_filter/.../top20` hanno un'etichetta contestuale ma una selezione globale.

Direzione di correzione:

- Il ranking `top20` deve essere calcolato sullo stesso `sqm` che verra' usato per generare gli output del contesto.
- Separare chiaramente `defined` globale da `top20` contestuale.

Test di regressione:

- Dataset sintetico con due taxa e pathway dominanti diversi: il `top20` del taxon deve cambiare rispetto al globale.
- Manifest e directory devono riflettere la selezione effettivamente usata.

### BUG-P1-03: join KO molti-a-molti raddoppia il TPM nei flow

Severita': alta. Stato: corretto e verificato su `fix/p1-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:1341-1392`: `build_flow_ko_metadata()` produce una riga per KO e concatena descrizioni distinte ordinate.
- `sqm_plots.R:1394-1455`: `join_flow_ko_metadata()` impone relazione molti-a-uno e postcondizioni su chiavi, righe e TPM per sample.
- `sqm_plots.R:1457-1496`: `build_flow_table_for_rank()` usa il lookup univoco senza duplicare gli archi.

Prova:

- Riproduzione minima osservata:
  - TPM reale: 30.
  - TPM prodotto: 60.
  - Archi reali: 2.
  - Righe prodotte: 4.
- Nel dataset reale sono stati osservati 123 KO con piu' descrizioni funzionali.

Impatto:

- `flow_percent` puo' sommare a 100 dopo la duplicazione perche' viene ricalcolato sul totale duplicato, ma la massa TPM degli archi e' falsa.
- Il Sankey e il TSV possono avere archi duplicati o gonfiati.

Direzione di correzione:

- La tabella metadata KO usata nel join deve essere una riga per KO.
- Conservare tutte le descrizioni distinte con una policy deterministica, per esempio valori unici ordinati e concatenati in un solo campo; non esplodere nuovamente le righe di abbondanza.
- Dopo il join, il TPM totale per sample deve restare invariato rispetto a prima del join.

Test di regressione:

- Test minimo con un KO associato a due `kegg_function`: il TPM post-join deve restare uguale al TPM pre-join.
- Assert esplicito in `build_flow_table_for_rank()` o nel test: nessun aumento di righe per chiave `sample/taxon/KO` causato dal metadata join.

## Evidenza di chiusura P1 - 2026-08-30

Branch: `fix/p1-correctness`, creato da `72115a9`.

Checkpoint TDD e implementazione:

- `dcfdc46`, `f8913f7`, `94d0821`: RED sintetici rispettivamente per gerarchia `top20`, selezione contestuale e join KO FLOW;
- `3d43c41`: RED di integrazione reale P1;
- `d4b19d2`: filtro delle sole gerarchie KEGG PATHWAY;
- `c681ab2`: calcolo di `top20` dentro ogni contesto SQM;
- `5e3a7b8`: metadata KO uno-a-uno e conservazione della massa FLOW;
- `13500e7`: gate di copertura P1;
- `f34eadf`: parsing gerarchico in batch, introdotto dopo che la prima integrazione ha evidenziato una regressione prestazionale sul progetto completo.

Risultati finali osservati, sempre tramite `rtk`:

- `tests/run_fast_tests.R`: exit `0`, 10 file di test veloci completati;
- `tests/check_p0_coverage.R`: exit `0`, funzioni P0 tra 84.21% e 100%; baseline globale informativa 9.74%;
- `tests/check_p1_coverage.R`: exit `0`, funzioni P1 tra 84.62% e 100%; baseline globale informativa P1 13.59%; `covr` 3.6.5;
- `tests/test_p0_integration_Au_sip.R`: exit `0`, 88258 ORF Bacillota e invarianti P0 ancora verdi;
- `tests/test_p1_integration_Au_sip.R`: exit `0`, `global_top20=20`, `Bacillota_top20=20`, confronto con riferimento indipendente superato e primi tre Bacillota `Quorum sensing`, `ABC transporters`, `Two-component system`;
- sul pathway `00361`, chiavi e TPM FLOW per sample restano invariati dopo il join dei metadata KO;
- il warning di compatibilita' progetto SqueezeMeta 1.7.3.alpha3/SQMtools 1.7.2 resta visibile;
- run CLI mirata: exit `0` in `out/p1_candidate_20260830_01`, nove TSV FLOW e nove PNG per i tre pathway Bacillota; tutti i target del manifest esistono, le chiavi `sample/taxon/KO` sono univoche e `flow_percent` somma a 100 entro `1e-6` in ogni TSV;
- review indipendente: nessun rilievo CRITICAL, HIGH o MEDIUM.

Non e' stata eseguita una run completa `--mode=all`: P2 resta aperto. La run candidata P1 non approva gli output storici. Il difetto dei nomi PNG troncati su percorsi Windows lunghi resta P3; il fatto che i nomi brevi del candidato P1 esistano non costituisce una correzione di quel difetto.

## P2 - Corretto su `fix/p2-correctness`

### BUG-P2-01: `Unclassified` confluisce in `Other`

Severita': media-alta. Stato: corretto e verificato su `fix/p2-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:504-549`: helper condivisi per Top N sui soli taxa classificati, categoria `Unclassified` preservata e valore sorgente `Other` rifiutato.
- `sqm_plots.R:1755-1774`: FLOW applica il collasso tassonomico condiviso.
- `sqm_plots.R:2210-2224`: PIE applica la stessa semantica.
- `tests/test_p2_reserved_taxa.R`: regressione sintetica comune a FLOW e PIE.
- `REGOLE_SCRIPT_R_SQUEEZEMETA.md:226-238`: `Unclassified` e `Other` devono restare categorie distinte.

Impatto:

- Si perde informazione biologica: non classificato non equivale a classificato ma fuori Top N.

Test di regressione:

- Dataset dove `Unclassified` e' fuori Top N: deve rimanere categoria distinta e non sommarsi a `Other`.

Direzione di correzione:

- Escludere `Unclassified` dal collasso Top N e riservargli sempre un livello distinto, sia nei FLOW sia nei PIE.
- Usare `Other` soltanto per taxa classificati ma fuori Top N.

### BUG-P2-02: associazioni multi-EC conservate solo parzialmente

Severita': media. Stato: corretto e verificato su `fix/p2-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:485-497`: estrazione esclusiva dal blocco `[EC:...]` di `KEGGFUN`.
- `sqm_plots.R:1194-1235`: `extract_ko_ec_lookup()` produce una riga per KO con codici distinti, ordinati e concatenati.
- `sqm_plots.R:1259-1455`: FUNZ usa il lookup molti-a-uno con postcondizioni su chiavi, righe e TPM.
- `sqm_plots.R:2185-2240`: PIE usa e propaga lo stesso elenco EC completo.
- `tests/test_p2_multi_ec.R`: regressioni per codici duplicati, mancanti, incompleti e join conservativi.

Prova:

- Il valore storico `78` non aveva una granularita' documentata. Il controllo token-level indipendente corrente osserva 493 KO con piu' EC distinti nel progetto completo e 3 nel pathway `00710`; l'integrazione ricalcola il valore senza hardcodarlo nello script.

Impatto:

- Legende, subtitle, file name e TSV possono mostrare un solo EC anche quando il KO e' multi-EC.

Direzione di correzione:

- Costruire un lookup a una riga per KO che conservi tutti gli EC distinti in ordine deterministico.
- Evitare un join uno-a-molti: la completezza dei metadati non deve duplicare il TPM.

Test di regressione:

- KO con due EC: il TSV e la legenda devono preservare l'informazione multi-EC o dichiarare una policy deterministica.

### BUG-P2-03: sample senza TPM positivo interrompe FUNZ

Severita': media. Stato: corretto e verificato su `fix/p2-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:1259-1455`: `build_ko_plot_table()` aggiunge denominatore, stato, flag `plotted` e sentinelle zero-denominator.
- `sqm_plots.R:2309-2383`: FUNZ scrive sempre il TSV e salta soltanto il PNG quando tutti i sample sono vuoti.
- `tests/test_p2_funz_zero.R`: copre sample positivo+zero, asse completo e caso all-zero.

Impatto:

- Un singolo sample senza TPM positivo nel pathway puo' fermare l'intera sezione FUNZ invece di produrre warning/skip controllato.

Test di regressione:

- Dataset con un sample vuoto e uno non vuoto: lo script deve completare e registrare chiaramente il sample saltato o vuoto.

### BUG-P2-04: CLI accetta interi frazionari troncati

Severita': media. Stato: corretto e verificato su `fix/p2-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:577-598`: `parse_positive_integer_arg()` valida la stringa grezza prima della coercizione.
- `sqm_plots.R:3312-3326`: i tre argomenti interi sono risolti prima di creare `output_dir` e prima di `loadSQM()`.
- `tests/test_p2_cli_integer.R`: copre valori validi, frazionari, notazione scientifica, segni, whitespace, zero e overflow, inclusa la precedenza sugli effetti collaterali.

Prova:

- `--top_n_ko=1.9` diventa `1` e supera `validate_positive_integer()`.

Impatto:

- Input CLI non valido produce risultati validi solo per troncamento implicito.

Test di regressione:

- `--top_n_ko=1.9`, `--top_n_taxa=2.5`, `--pathway_top_n=3.1` devono fallire con messaggio chiaro.

### BUG-P2-05: policy multi-KO ed esclusioni ORF non compaiono nei manifest

Severita': media. Stato: corretto e verificato su `fix/p2-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:374-413`: validazione e normalizzazione dell'audit destinato ai manifest.
- `sqm_plots.R:415-482`: `new_manifest_row()` aggiunge i sei campi P2 in modo append-only.
- `sqm_plots.R:1064-1177`: audit pre-espansione strutturato e tabella ORF x sample x KO con TPM intero per associazione.
- `sqm_plots.R:2340-2379`, `2603-2677` e `3051-3152`: propagazione nei manifest FUNZ, FLOW e PIE.
- `tests/test_p2_ko_provenance.R`: conteggi indipendenti, replica TPM e round-trip dei manifest dei tre modi.
- `REGOLE_SCRIPT_R_SQUEEZEMETA.md:208-222`: il comportamento multi-KO deve essere esplicito e tracciabile.

Impatto:

- I manifest non documentano quante ORF sono state escluse ne' che il TPM viene replicato per ogni KO associato.

Nota concettuale da rendere esplicita:

- `REGOLE_SCRIPT_R_SQUEEZEMETA.md:208-222` stabilisce gia' la policy del progetto: ogni KO associato riceve l'intero TPM dell'ORF, senza frazionamento. Questa scelta non conserva la massa quando si somma tra KO, quindi i denominatori KO devono essere calcolati sulla tabella espansa e la policy deve essere dichiarata. Non cambiarla silenziosamente durante il fix: il difetto confermato e' che oggi policy e numero di ORF escluse non sono tracciati.

Direzione di correzione:

- Portare `excluded_orfs_without_ko` fuori dall'attributo transitorio e inserirlo nel manifest o in un audit TSV.
- Aggiungere al manifest un campo stabile per la policy multi-KO e, se rilevante, il numero di associazioni prodotte.

Test di regressione:

- Manifest o audit TSV devono includere policy multi-KO e conteggio ORF senza KO per sezione/pathway.

## Evidenza di chiusura P2 - 2026-08-30

Branch: `fix/p2-correctness`, creato da `812171f`.

Checkpoint TDD e implementazione:

- `ab0d9ea`: preparazione harness P2;
- `07e773c` / `57d2290`: RED/GREEN per categorie tassonomiche riservate;
- `facc8f7` / `840e456`: RED/GREEN per lookup multi-EC deterministico;
- `ef63d06` / `9ee2ef2`: RED/GREEN per sample FUNZ con denominatore zero;
- `30a97cf` / `f61454f`: RED/GREEN per parsing CLI degli interi;
- `6089c40`, `0aa0955` / `116d4c2`: RED e GREEN per audit multi-KO e propagazione nei manifest;
- `8a194aa`: compatibilita' del lookup EC con fixture tassonomiche prive della colonna `ec_codes`;
- `39f3848`: audit KO estratto senza espandere l'intero SQM;
- `1a00ee3`: integrazione reale P2;
- `38b9daa`: gate di copertura P2;
- `b0db380`: verificatore parametrico delle tre run candidate.

Risultati finali osservati, sempre tramite `rtk`:

- `tests/run_fast_tests.R`: exit `0`, 15 file di test veloci completati;
- `tests/check_p0_coverage.R`: exit `0`, funzioni P0 tra 84.21% e 100%, baseline globale informativa 8.94%;
- `tests/check_p1_coverage.R`: exit `0`, funzioni P1 tra 84.62% e 100%, baseline globale informativa 13.34%;
- `tests/check_p2_coverage.R`: exit `0`, funzioni P2 tra 85.71% e 100%, baseline globale informativa 16.09%; `covr` 3.6.5;
- `tests/test_p0_integration_Au_sip.R`: exit `0`, 88258 ORF Bacillota e invarianti P0 ancora verdi;
- `tests/test_p1_integration_Au_sip.R`: exit `0`, `global_top20=20`, `Bacillota_top20=20` e invarianti P1 ancora verdi;
- `tests/test_p2_integration_Au_sip.R`: exit `0`, `multi_ec_00710=3`, `excluded_orfs_without_ko=705510`, `multi_ko_orfs=1317`; distinzione FLOW `Unclassified`/`Other`, lookup EC e audit confrontati con riferimenti indipendenti;
- il warning di compatibilita' progetto SqueezeMeta 1.7.3.alpha3/SQMtools 1.7.2 e' rimasto visibile in tutte le integrazioni reali;
- run FLOW `00361` in `out/p2_candidate_20260830_01_flow`: exit `0`, tre TSV e tre PNG, 114 righe FLOW complessive, chiavi univoche e somme `flow_percent=100` per sample;
- run FUNZ `00710` in `out/p2_candidate_20260830_01_funz`: exit `0`, 99 righe pathway, denominatori coerenti e 3 KO multi-EC completi;
- run PIE `00633` / `S13_1_8` in `out/p2_candidate_20260830_01_pie`: exit `0`, 11 TSV e 11 PNG; `K10679` conserva `1.-.-.-;1.5.1.34` nel TSV e nel manifest;
- `tests/verify_p2_candidate_outputs.R`: exit `0` sulle tre directory candidate;
- review spec indipendente: nessun rilievo CRITICAL, HIGH, MEDIUM o LOW; review standard: nessun rilievo CRITICAL/HIGH.

Note di review non bloccanti: `sqm_plots.R` resta un monolite preesistente oltre il limite generale indicato da AGENTS.md e merita un refactor separato; nei contesti FLOW/PIE senza alcuna ORF con KO l'audit viene calcolato ma, non essendoci alcun artefatto, non nasce una riga manifest che lo serializzi. Quest'ultimo caso limite non riguarda i flussi P2 richiesti e puo' essere reso esplicito in un futuro audit TSV senza cambiare il significato dei manifest come inventario di file.

Non e' stata eseguita una run `--mode=all`. Gli output storici non sono stati modificati o approvati. Le directory candidate P2 restano evidenza locale, non output storici approvati. Il problema Windows dei nomi PNG troncati resta P3 e non e' stato corretto ne' usato per ampliare P2.

## P3 - Igiene output, dipendenze e portabilita'

### BUG-P3-01: manifest storici mantengono riferimenti a file mancanti

Severita': media. Stato: verificato.

Riferimenti:

- `sqm_plots.R:2379-2389`: `read_section_manifest()` legge manifest esistente.
- `sqm_plots.R:2392-2394`: `merge_section_manifest()` fa `bind_rows()` e `distinct(output_file)` senza verificare l'esistenza del target.
- `sqm_plots.R:2408-2414`: `write_section_manifest()` riscrive il manifest unito.

Prova:

- L'audit precedente ha registrato 15 riferimenti a file mancanti nei manifest. Il difetto strutturale che li conserva e' stato riconfermato nel codice; ripetere il conteggio sulla directory scelta come baseline prima del fix, perche' il contenuto di `out/` puo' cambiare tra sessioni.

Impatto:

- Il manifest non e' un inventario affidabile degli output presenti.

Test di regressione:

- Dopo ogni run, ogni `output_file` nei manifest deve esistere rispetto alla base del manifest.

### BUG-P3-02: TSV PIE incompleto rispetto ai dati usati nel grafico

Severita': media. Stato: verificato.

Riferimenti:

- `sqm_plots.R:1392-1432`: `build_pie_chart_table()` produce taxa, TPM e percentuali.
- `sqm_plots.R:2243-2313`: `run_pie_mode()` usa anche pathway denominator, contributo KO, nome KO ed EC per subtitle/caption/file stem.

Impatto:

- Il TSV non basta a ricostruire completamente il grafico e le sue annotazioni.

Test di regressione:

- TSV PIE deve contenere almeno denominatore pathway del sample, contributo KO al pathway, nome KO ed EC usati nel grafico.

### BUG-P3-03: `pathview_is_exportable("defined", NA)` restituisce true

Severita': media. Stato: corretto insieme a BUG-P0-01 su `fix/p0-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:275-276`: `pathview_is_exportable()`.
- `sqm_plots.R:2111-2116`: `run_pathview_mode()` esegue `stop()` se `pathway_id` e' `NA` o vuoto.
- `tests/test_pathway_selection.R:85-86`: test corrente accetta `defined, NA` come exportable.

Impatto:

- Un pathway valido per nome completo ma fuori dai nove curati puo' superare il gate in `run_pipeline()` e arrivare a `run_pathview_mode()` senza ID esportabile, causando un errore invece di uno skip controllato.

Direzione di correzione:

- Definire l'esportabilita' in base alla presenza di un `pathway_id` KEGG valido, non in base al solo tipo `defined`/`top20`.
- Se il nome e' risolvibile in modo affidabile, assegnare l'ID prima del gate; altrimenti emettere warning e saltare Pathview senza fermare le altre modalita'.

Test di regressione:

- Casi minimi espliciti: `top20 + NA -> FALSE`, `defined + NA -> FALSE`, `defined + "00710" -> TRUE`.
- Un pathway senza ID deve produrre uno skip controllato della sola sezione Pathview, non uno `stop()` dell'intera run.

### BUG-P3-04: preflight package inefficace e incompleto

Severita': media-bassa. Stato: verificato.

Riferimenti:

- `sqm_plots.R:1-40`: i package sono caricati prima di `check_required_packages()`.
- `sqm_plots.R:2421-2423`: `main()` chiama `check_required_packages()` solo dopo il caricamento globale.
- Package usati ma da verificare nel preflight: `forcats`, `rlang`, dipendenze Pathview.

Impatto:

- Se manca un package caricato in testa, lo script fallisce prima del messaggio di preflight.

Test di regressione:

- Preflight deve fallire con messaggio controllato quando manca una dipendenza dichiarata.

### BUG-P3-05: test con percorso personale assoluto

Severita': bassa. Stato: verificato.

Riferimento:

- `tests/test_windows_project_path.R:12`.

Impatto:

- Il fixture incorpora un percorso e un nome utente personali. Il path oggi e' trattato solo come stringa e quindi il test puo' comunque passare altrove, ma l'intento di normalizzazione Windows puo' essere verificato senza dipendere ne' divulgare la struttura del computer dell'autore.

Test di regressione:

- Usare un path Windows sintetico e non personale; nessun test deve richiedere che un percorso assoluto specifico dell'autore esista.

### BUG-P3-06: la suite corrente non esercita i percorsi che hanno prodotto risultati falsi

Severita': media. Stato: verificato.

Riferimenti:

- `tests/test_pathway_selection.R`.
- `tests/test_pathview_sample_modes.R`.
- `tests/test_enzyme_mode.R`.
- `tests/test_windows_project_path.R`.

Copertura mancante osservata:

- conservazione del TPM nei FLOW;
- denominatore delle percentuali tassonomiche per pathway;
- filtro tassonomico ORF rispetto al contig;
- distinzione `Unclassified`/`Other`;
- sample con TPM pathway nullo;
- rifiuto di interi CLI frazionari;
- integrita' dei target nei manifest;
- una verifica integrata sui dati reali in `in/Au_sip`.

Impatto:

- I quattro test passano pur in presenza dei bug P0 e P1; il verde corrente non costituisce una validazione scientifica degli output.

Direzione di correzione:

- Scrivere prima una riproduzione minima fallente per ogni fix e mantenerla nella suite.
- Separare test sintetici veloci dalla run di integrazione reale, ma rendere entrambe ripetibili e documentate.

## Ordine di implementazione per la prossima sessione

### Fase 0 - Bloccare i falsi risultati critici

1. Aggiungere test fallenti per mapping KEGG, percentuali pathway e filtro ORF.
2. Correggere `known_pathways` e il test che codifica `00643`.
3. Calcolare esplicitamente le percentuali tassonomiche sul totale pathway/sample, con gestione del denominatore zero.
4. Sostituire il subset contig-based con un subset costruito dall'insieme esatto di ORF risolto.
5. Eseguire test sintetici e controlli reali Bacillota prima di procedere.

### Fase 1 - Ripristinare il significato di `top20` e la massa FLOW

1. Conservare la radice gerarchica di `KEGGPATH` e filtrare BRITE/non-pathway prima del ranking.
2. Calcolare `top20` dentro ogni `filter_context`, non una volta sul progetto globale.
3. Rendere il lookup descrittivo KO uno-a-uno e aggiungere un'invariante di conservazione TPM pre/post join.

### Fase 2 - Sistemare categorie, metadati e input limite

1. Separare sempre `Unclassified` da `Other`.
2. Definire una rappresentazione deterministica multi-EC e rendere esplicita nei metadati la policy multi-KO gia' fissata dalle regole del progetto.
3. Gestire sample senza TPM positivo con warning/skip tracciato.
4. Validare la stringa numerica CLI prima di convertirla in intero.

### Fase 3 - Rendere output e suite verificabili

1. Scrivere solo riferimenti esistenti nei manifest e verificare anche `manifest_all.tsv`.
2. Aggiungere al TSV PIE tutti i dati usati per titolo, sottotitolo e caption.
3. Correggere il gate Pathview per ID assenti.
4. Spostare il preflight prima dei `library()` o usare accessi namespaced con una lista completa e mode-aware delle dipendenze.
5. Rendere portabile il test Windows e aggiungere i test mancanti elencati in BUG-P3-06.

Regola di avanzamento: non iniziare una fase dichiarando chiusa la precedente finche' i suoi test sintetici non sono verdi e l'invariante reale associata non e' stata controllata. Non modificare i test per adattarli a un risultato noto errato.

## Comandi di verifica consigliati

Usare sempre `rtk` nella shell della sessione.

```bat
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\test_pathway_selection.R
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\test_pathview_sample_modes.R
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\test_enzyme_mode.R
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\test_windows_project_path.R
```

Eseguire poi i nuovi test aggiunti per ciascun bug. Per l'integrazione reale usare una directory nuova:

```bat
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe sqm_plots.R --project_dir=in\Au_sip --output_dir=out\regression_audit_fixed --mode=all --pathway_selection_modes=defined,top20
```

Non riutilizzare `out\regression_audit_fixed` per confrontare implementazioni diverse: scegliere un nuovo suffisso per ogni run candidata. La run completa puo' essere costosa; durante lo sviluppo usare test sintetici e una run mirata sul pathway `00361`, poi eseguire `--mode=all` come gate finale.

Controlli manuali o scriptabili dopo la run:

- Sommare `value` per `sample` nei TSV `taxonomy_by_pathway/.../percent/...`: ogni somma deve essere 100 entro tolleranza.
- Verificare che ogni `output_file` nei manifest esista.
- Verificare che i TSV flow conservino il TPM totale prima e dopo il join dei metadati KO.
- Verificare che `top20` non contenga voci BRITE o categorie non-pathway.
- Verificare che gli output sotto `taxon_filter/...` derivino esattamente dall'insieme ORF del taxon richiesto, senza inclusioni **e senza esclusioni**.
- Registrare versioni di R, SQMtools e progetto SqueezeMeta insieme alla run candidata.

## Checklist prossima sessione

- [x] Aggiornare test di regressione per BUG-P0-01, BUG-P0-02 e BUG-P0-03 prima o insieme alle correzioni.
- [x] Correggere le mappature KEGG e aggiornare i test che codificano ID errati.
- [x] Ricalcolare le percentuali tassonomiche pathway con denominatore pathway/sample.
- [x] Rendere il filtro taxon coerente con la sorgente ORF usata per risolvere il taxon.
- [x] Limitare `top20` a vere pathway map.
- [x] Calcolare `top20` nel contesto filtrato quando `--taxa` e' attivo.
- [x] Rendere uno-a-uno il metadata join KO nei flow e verificare conservazione TPM.
- [ ] Tenere `Unclassified` distinto da `Other`.
- [ ] Definire e documentare policy multi-EC e multi-KO nei TSV/manifest.
- [ ] Gestire sample pathway vuoti senza interrompere l'intera run.
- [ ] Validare input CLI numerici prima di `as.integer()`.
- [ ] Pulire i manifest da target inesistenti o fallire esplicitamente.
- [ ] Ricontare e registrare i riferimenti manifest mancanti sulla baseline usata nella nuova sessione.
- [ ] Rendere completo il TSV PIE rispetto alle annotazioni del grafico.
- [ ] Rendere il preflight package effettivo.
- [ ] Rimuovere path personali dai test.
- [ ] Eseguire tutti i test esistenti e nuovi, quindi una run reale pulita su `in/Au_sip`.
- [ ] Rigenerare gli output affetti in una directory nuova; non promuovere quelli storici.

Done only when:

- I mapping KEGG sono conformi agli ID ufficiali.
- Le percentuali tassonomiche per pathway sommano a 100 per ogni sample con denominatore positivo; i sample vuoti sono gestiti esplicitamente.
- Il subset Bacillota coincide esattamente con gli ID ORF attesi.
- `top20` contiene solo vere foglie pathway e viene calcolato nel contesto corretto.
- Il TPM degli archi flow resta invariato dopo i join.
- `Unclassified` resta distinto da `Other`.
- Ogni target nei manifest esiste.
- I TSV contengono i metadati e i denominatori necessari a ricostruire i grafici.
- La policy multi-KO e' esplicita e i multi-EC non vengono troncati.
- I test vecchi e nuovi passano sui dati sintetici e la run reale pulita `in/Au_sip` completa senza falsi riferimenti nei manifest.
