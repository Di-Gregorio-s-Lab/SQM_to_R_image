# Audit bug `sqm_plots.R`

> Stato del documento: audit e registro di chiusura del 2026-08-30. P0, P1, P2 e P3 sono stati corretti e verificati sui rispettivi branch, come annotato nelle sezioni di evidenza. Gli output storici restano prove read-only e non sono stati approvati retroattivamente. I numeri di linea originari possono essersi spostati durante le correzioni: usare anche il nome della funzione indicato.

Questo documento ha guidato la correzione di `sqm_plots.R` e ora ne conserva le evidenze. Per ogni problema riporta il riferimento al codice, la prova osservata, l'effetto sui risultati, la correzione e il test che deve impedirne la ricomparsa.

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

## Indice rapido dei problemi auditati

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

Limitazione non-P0 riconfermata al checkpoint P0: sul lungo percorso Windows del workspace il nome del PNG pathway era stato troncato e il target previsto nel manifest non esisteva. Il difetto restava allora sotto BUG-P3-01/portabilita' degli output; non alterava il TSV e le invarianti scientifiche usate per chiudere P0, ma impediva di considerare l'intera directory candidata un output finale approvato. E' stato successivamente corretto e verificato in P3.

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

Al checkpoint P1 non era stata eseguita una run completa `--mode=all` e P2 restava aperto. La run candidata P1 non approvava gli output storici. Il difetto dei nomi PNG troncati su percorsi Windows lunghi restava P3; il fatto che i nomi brevi del candidato P1 esistessero non costituiva una correzione. Il difetto e' stato successivamente corretto e verificato in P3.

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

Al checkpoint P2 non era stata eseguita una run `--mode=all`. Gli output storici non erano stati modificati o approvati. Le directory candidate P2 restano evidenza locale, non output storici approvati. Il problema Windows dei nomi PNG troncati restava P3 e non era stato corretto ne' usato per ampliare P2; e' stato successivamente corretto e verificato in P3.

## P3 - Corretto su `fix/p3-correctness`

### BUG-P3-01: manifest storici mantengono riferimenti a file mancanti

Severita': media. Stato: corretto e verificato su `fix/p3-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:3451-3550`: `manifest_target_status()`, `validate_current_manifest_targets()` e `prune_stale_manifest_targets()` distinguono target nuovi non validi da righe storiche stale.
- `sqm_plots.R:3553-3638`: merge, scrittura delle sezioni, riconciliazione completa e rigenerazione di `manifest_all.tsv` operano soltanto su target relativi, regolari, esistenti e non vuoti.
- `sqm_plots.R:412-515`: token deterministico, budgeting del path Windows e postcondizione successiva a `ggsave()` impediscono che un PNG troncato entri nel manifest.
- `tests/test_p3_manifest_integrity.R` e `tests/test_p3_portable_png.R`: regressioni per target stale/nuovi, traversal, atomicita' del manifest, collisioni, determinismo e soglia di 240 caratteri.

Prova:

- La baseline storica read-only ricontata prima del fix contiene 18 target mancanti nei manifest di sezione e un riferimento a manifest mancante in `manifest_all.tsv`. Nessuno di questi file o manifest storici e' stato modificato.
- Il nuovo verificatore parametrico ha osservato zero target mancanti nei tre candidati P3.

Impatto:

- Il manifest non e' un inventario affidabile degli output presenti.

Correzione verificata:

- I target della run corrente vengono validati prima della riscrittura; un target nuovo mancante, assoluto o fuori root interrompe la run lasciando invariato il manifest esistente.
- Le righe storiche stale vengono rimosse con warning senza cancellare o ricreare gli output interessati.
- Tutte le sezioni presenti vengono riconciliate e `manifest_all.tsv` viene rigenerato da zero includendo soltanto manifest validi e non vuoti.
- I PNG mantengono il nome logico quando il path previsto non supera 240 caratteri; oltre soglia viene accorciato soltanto il filename con token deterministico a 12 cifre esadecimali e viene verificato il file fisico esatto e non vuoto.

### BUG-P3-02: TSV PIE incompleto rispetto ai dati usati nel grafico

Severita': media. Stato: corretto e verificato su `fix/p3-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:2338-2428`: `build_pie_chart_table()` esporta, in modo append-only, pathway, ID, selezione, rank, nome KO, denominatori, contributo KO e ordine tassonomico.
- `sqm_plots.R:2432-2490`: `make_pie_plot(plot_data)` riceve soltanto la tabella esportata e deriva da essa titolo, sottotitolo, caption, legenda e ordine.
- `sqm_plots.R:3320-3370`: `run_pie_mode()` scrive il TSV e costruisce il grafico dalla stessa tabella esportata; i test la rileggono per verificare il round-trip.
- `tests/test_p3_pie_export.R` e `tests/test_p3_integration_Au_sip.R`: regressione sintetica e confronto reale dei denominatori calcolati indipendentemente.

Impatto:

- Il TSV non basta a ricostruire completamente il grafico e le sue annotazioni.

Correzione verificata:

- Il TSV PIE conserva le colonne precedenti e aggiunge `pathway`, `pathway_id`, `pathway_selection`, `rank`, `ko_name`, `ko_sample_tpm`, `pathway_sample_tpm`, `ko_pathway_percent` e `taxon_order`.
- Un round-trip TSV produce le stesse annotazioni, lo stesso ordine e la stessa massa TPM della tabella in memoria.

### BUG-P3-03: `pathview_is_exportable("defined", NA)` restituisce true

Severita': media. Stato: corretto insieme a BUG-P0-01 su `fix/p0-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:347-356`: `pathview_is_exportable()` richiede un ID conforme a `^[0-9]{5}$`.
- `sqm_plots.R:4024-4035`: la pipeline emette warning e salta soltanto Pathview quando l'ID non e' esportabile.
- `tests/test_p0_pathway_mapping.R` e `tests/test_pathway_selection.R`: `top20 + NA` e `defined + NA` sono falsi, mentre `defined + "00910"` e' vero.

Impatto:

- Un pathway valido per nome completo ma fuori dai nove curati puo' superare il gate in `run_pipeline()` e arrivare a `run_pathview_mode()` senza ID esportabile, causando un errore invece di uno skip controllato.

Correzione verificata:

- Definire l'esportabilita' in base alla presenza di un `pathway_id` KEGG valido, non in base al solo tipo `defined`/`top20`.
- Se il nome e' risolvibile in modo affidabile, assegnare l'ID prima del gate; altrimenti emettere warning e saltare Pathview senza fermare le altre modalita'.

Test di regressione mantenuto:

- Casi minimi espliciti: `top20 + NA -> FALSE`, `defined + NA -> FALSE`, `defined + "00710" -> TRUE`.
- Un pathway senza ID deve produrre uno skip controllato della sola sezione Pathview, non uno `stop()` dell'intera run.

### BUG-P3-04: preflight package inefficace e incompleto

Severita': media-bassa. Stato: corretto e verificato su `fix/p3-correctness` (2026-08-30).

Riferimenti:

- `sqm_plots.R:7-29`: elenco comune e `required_packages_for_mode()` selezionano le dipendenze contestuali di FLOW, HTML, PIE e Pathview.
- `sqm_plots.R:29-47`: `check_required_packages()` aggrega le dipendenze mancanti in un unico errore controllato.
- `sqm_plots.R:49-86`: parsing bootstrap base-R, help/no-arg e preflight vengono eseguiti prima di qualsiasi `library()`.
- `tests/test_p3_preflight.R`: subprocess con libreria utente vuota verifica help/no-arg e fallimento analitico prima di `output_dir` e `loadSQM()`.

Impatto:

- Se manca un package caricato in testa, lo script fallisce prima del messaggio di preflight.

Correzione verificata:

- `--help` e l'invocazione senza argomenti funzionano senza caricare dipendenze analitiche.
- Le dipendenze comuni e contestuali vengono verificate prima della creazione di `output_dir` e prima di `loadSQM()`; una mancanza produce un solo messaggio con modalita' e package coinvolti.

### BUG-P3-05: test con percorso personale assoluto

Severita': bassa. Stato: corretto e verificato su `fix/p3-correctness` (2026-08-30).

Riferimento:

- `tests/test_windows_project_path.R`: usa una root temporanea sintetica con spazi e non contiene percorsi personali.

Impatto:

- Il fixture incorpora un percorso e un nome utente personali. Il path oggi e' trattato solo come stringa e quindi il test puo' comunque passare altrove, ma l'intento di normalizzazione Windows puo' essere verificato senza dipendere ne' divulgare la struttura del computer dell'autore.

Correzione verificata:

- Nessun test contiene `C:\Users\unico`; il fixture Windows non richiede che un percorso assoluto specifico dell'autore esista.

### BUG-P3-06: la suite corrente non esercita i percorsi che hanno prodotto risultati falsi

Severita': media. Stato: corretto e verificato progressivamente sui branch `fix/p0-correctness` - `fix/p3-correctness` (2026-08-30).

Riferimenti:

- `tests/run_fast_tests.R`: 19 file sintetici, con le integrazioni reali P0-P3 escluse esplicitamente.
- `tests/check_p0_coverage.R` - `tests/check_p3_coverage.R`: gate mirati `covr` sulle funzioni corrette.
- `tests/test_p0_integration_Au_sip.R` - `tests/test_p3_integration_Au_sip.R`: invarianti reali caricate da `in/Au_sip`.
- `tests/verify_p2_candidate_outputs.R` e `tests/verify_p3_candidate_outputs.R`: controllo parametrico degli artefatti candidati.

Copertura originariamente mancante e ora mantenuta:

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

Correzione verificata:

- Ogni slice P0-P3 conserva la riproduzione sintetica introdotta prima del fix.
- Test veloci, coperture mirate, integrazioni reali e verificatori candidate sono separati, ripetibili e documentati nelle evidenze di chiusura.

## Evidenza di chiusura P3 - 2026-08-30

Branch: `fix/p3-correctness`, creato da `01a31fd`.

Checkpoint TDD e implementazione:

- `ab2b49f`: preparazione harness P3;
- `54e1c9a` / `f2fa1e0`: RED/GREEN per integrita' dei target manifest;
- `a9884c1` / `18cf930`: RED/GREEN per nomi PNG portabili su Windows;
- `b086d7d` / `654e1ea`: RED/GREEN per export PIE completo e round-trip tabella-grafico;
- `827f770` / `130e5d5`: RED/GREEN per preflight package completo e mode-aware;
- `47b7af5`: fixture Windows sintetico senza percorso personale;
- `47880e9`: riconciliazione sicura delle righe manifest storiche;
- `a1b3993`: integrazione reale, copertura P3 e verificatore parametrico dei candidati;
- `6f2e2cc`: fix emerso in review per manifest legacy malformati senza colonna `output_file`.

Comandi finali eseguiti tramite `rtk`, tutti con exit `0`:

```powershell
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\run_fast_tests.R
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\check_p0_coverage.R
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\check_p1_coverage.R
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\check_p2_coverage.R
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\check_p3_coverage.R
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\test_p0_integration_Au_sip.R
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\test_p1_integration_Au_sip.R
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\test_p2_integration_Au_sip.R
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\test_p3_integration_Au_sip.R
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\verify_p3_candidate_outputs.R out\p3_candidate_01_windows out\p3_candidate_01_pie out\p3_candidate_01_all
rtk proxy C:\Progra~1\R\R-4.5.0\bin\Rscript.exe tests\test_p3_manifest_integrity.R
```

Risultati osservati:

- test veloci: 19 file completati;
- `covr` 3.6.5: tutte le funzioni mirate P0-P3 superano l'80%; baseline globali informative rispettivamente 8.19%, 12.21%, 15.62% e 15.09%;
- integrazioni P0, P1 e P2 ancora verdi senza regressioni;
- integrazione P3 reale su `00633` / `S13_1_8`: `K10679`, KO TPM `18.284`, pathway TPM `25.632`, confrontati con denominatori calcolati indipendentemente;
- warning di compatibilita' SqueezeMeta 1.7.3.alpha3 / SQMtools 1.7.2 ancora visibile;
- candidato Windows `out/p3_candidate_01_windows`: 8 target validi, inclusi 6 PNG; nessun file troncato senza estensione;
- candidato PIE `out/p3_candidate_01_pie`: 22 target validi, inclusi 11 TSV PIE completi;
- candidato finale `out/p3_candidate_01_all`: 56 target validi e cinque sezioni (`flow`, `funz`, `taxon`, `pathview`, `pie`) in `manifest_all.tsv`;
- Pathview live per `00361` completato senza mock e tutti i target Pathview registrati esistono e sono non vuoti;
- la run ridotta ma reale `--mode=all` e' completata con successo; non era mai stata eseguita nelle chiusure P0-P2.
- review finale owner/spec: nessun rilievo CRITICAL o HIGH residuo. Un rilievo medio emerso durante la chiusura riguardava i manifest legacy privi di `output_file`; e' stato chiuso con `6f2e2cc` e test dedicato in `tests/test_p3_manifest_integrity.R`.

Il primo tentativo richiesto in `out/p3_candidate_20260830_01_windows` ha esercitato correttamente il fail-fast: la root assoluta era gia' lunga 222 caratteri e non poteva contenere prefisso minimo, token, dimensione ed estensione entro il limite 240. La directory parziale e' rimasta come evidenza locale e non e' stata classificata come candidata verde. Sono quindi state usate le tre root piu' corte sopra, senza cambiare la soglia o le directory interne.

Baseline storica AC-008, conservata read-only: 18 target mancanti nei manifest di sezione e un riferimento a manifest mancante in `manifest_all.tsv`. Gli output storici non sono stati modificati, rigenerati o approvati. Restano fuori scope il refactor del monolite `sqm_plots.R` e la serializzazione dell'audit nei contesti FLOW/PIE completamente privi di KO e quindi privi di artefatti.

Review finale read-only del 2026-08-31: l'asse di specifica ha confermato la risoluzione dei rilievi documentali dopo la correzione dell'evidenza PIE; l'asse standard non ha rilevato violazioni. Non restano rilievi CRITICAL/HIGH. Il monolite `sqm_plots.R` e' stato riconfermato come debito preesistente fuori scope, non come blocco P3.

## Traccia di implementazione completata

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

La sequenza e' stata eseguita per slice RED-GREEN, senza modificare i test per adattarli a risultati noti errati. Ogni fase e' stata chiusa soltanto dopo i test sintetici e l'invariante reale associata.

## Comandi di verifica originari

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

## Checklist completata

- [x] Aggiornare test di regressione per BUG-P0-01, BUG-P0-02 e BUG-P0-03 prima o insieme alle correzioni.
- [x] Correggere le mappature KEGG e aggiornare i test che codificano ID errati.
- [x] Ricalcolare le percentuali tassonomiche pathway con denominatore pathway/sample.
- [x] Rendere il filtro taxon coerente con la sorgente ORF usata per risolvere il taxon.
- [x] Limitare `top20` a vere pathway map.
- [x] Calcolare `top20` nel contesto filtrato quando `--taxa` e' attivo.
- [x] Rendere uno-a-uno il metadata join KO nei flow e verificare conservazione TPM.
- [x] Tenere `Unclassified` distinto da `Other`.
- [x] Definire e documentare policy multi-EC e multi-KO nei TSV/manifest.
- [x] Gestire sample pathway vuoti senza interrompere l'intera run.
- [x] Validare input CLI numerici prima di `as.integer()`.
- [x] Pulire i manifest da target inesistenti o fallire esplicitamente.
- [x] Ricontare e registrare i riferimenti manifest mancanti sulla baseline usata nella nuova sessione.
- [x] Rendere completo il TSV PIE rispetto alle annotazioni del grafico.
- [x] Rendere il preflight package effettivo.
- [x] Rimuovere path personali dai test.
- [x] Eseguire tutti i test esistenti e nuovi, quindi una run reale pulita su `in/Au_sip`.
- [x] Generare i candidati corretti in directory nuove; non promuovere gli output storici.

Done verificato:

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
