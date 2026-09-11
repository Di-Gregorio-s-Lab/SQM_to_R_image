# `loadSQM()` — duplicati nei bin e `sum(nobin)`

## Esito

I due messaggi sono collegati. `loadSQM()` legge
`intermediate/18.<progetto>.contigsinbins` e avverte quando una stessa coppia
**contig + Method** ricorre più volte. Il loader costruisce poi le abbondanze
dei bin riallineandole alle righe di `results/18.<progetto>.bintable`.
Un bin della tabella che non è più rappresentabile dall'assegnazione dei
contig produce una riga `NA`; `colSums(x)` e quindi `sum(nobin)` diventano
`NA`, causando `valore mancante dove è richiesto TRUE/FALSE`.

Questo è un difetto/limite del loader: l'avvertimento segnala dati di binning
incoerenti, ma non interrompe il caricamento prima del calcolo con `NA`.

## Evidenza primaria

- In [SQMtools v1.7.2 — `loadSQM.R`](https://github.com/jtamames/SqueezeMeta/blob/v1.7.2/lib/SQMtools/R/loadSQM.R#L363-L395), il controllo è
  `max(table(paste(inBins$X..Contig, inBins$Method))) > 1`; subito dopo il
  loader crea `SQM$contigs$bins` e chiama `get.bin.abunds()`.
- In [SQMtools v1.7.2 — `bin_methods.R`](https://github.com/jtamames/SqueezeMeta/blob/v1.7.2/lib/SQMtools/R/bin_methods.R#L270-L278), il codice riallinea
  `x` ai nomi della bintable e poi esegue esattamente:

  ```r
  nobin = colSums(SQM$contigs$abund) - colSums(x)
  if(sum(nobin)>0) { x['No_bin',] = nobin }
  ```

  Se il riallineamento inserisce `NA`, la condizione `if` non è valutabile.
- Lo [script SqueezeMeta 18.getbins.pl](https://github.com/jtamames/SqueezeMeta/blob/v1.7.2/scripts/18.getbins.pl#L66-L78) ricrea entrambi i file;
  scrive una riga per ogni contig trovato in ciascun FASTA di bin
  ([righe 120–140](https://github.com/jtamames/SqueezeMeta/blob/v1.7.2/scripts/18.getbins.pl#L120-L140)). Bin FASTA sovrapposti o artefatti di
  restart possono quindi creare il duplicato.

## Verifica non distruttiva

Adattare `project` e verificare prima quale prefisso (`18` o `19`) esiste nella
propria versione:

```r
project <- "C:/.../COL_SA_mg_mt"
prefix <- "18"
in_bins <- read.delim(file.path(project, "intermediate",
  sprintf("%s.COL_SA_mg_mt.contigsinbins", prefix)), skip = 1,
  comment.char = "", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
bin_table <- read.delim(file.path(project, "results",
  sprintf("%s.COL_SA_mg_mt.bintable", prefix)), skip = 1,
  comment.char = "", quote = "", row.names = 1, check.names = FALSE)

contig_col <- grep("Contig", names(in_bins), value = TRUE)[1]
method_col <- "Method"
bin_col <- grep("Bin", names(in_bins), value = TRUE)[1]
key <- in_bins[c(contig_col, method_col)]
duplicates <- in_bins[duplicated(key) | duplicated(key, fromLast = TRUE), ]
assigned_missing_from_bintable <- setdiff(unique(in_bins[[bin_col]]), rownames(bin_table))
duplicates
assigned_missing_from_bintable
```

`contig_col` e `bin_col` tollerano le varianti di nome generate da R; verificare
comunque `names(in_bins)` prima di interpretare l'output.

## Rimedi supportati

1. Correggere la causa nel risultato di binning: per ciascun `Method`, un
   contig deve appartenere a un solo bin. Rigenerare insieme
   `contigsinbins` e `bintable` dalla stessa esecuzione/insieme di FASTA,
   preferibilmente ripetendo il passaggio di raccolta bin in uno stato pulito;
   non modificare soltanto uno dei due TSV.
2. Se il progetto è stato riavviato con binner o opzioni diverse, escludere
   output intermedi/Fasta obsoleti prima di rigenerare. Un caso ufficiale di
   errore in fase `binning info` dopo restart multipli da step 16 è discusso
   in [issue #128](https://github.com/jtamames/SqueezeMeta/issues/128).
3. Aggiornare SQMtools da solo non è una cura: lo stesso `sum(nobin)` è
   presente anche nel [sorgente master](https://github.com/jtamames/SqueezeMeta/blob/master/lib/SQMtools/R/bin_methods.R#L278-L286).
   Inoltre, la firma pubblica di [`loadSQM()`](https://github.com/jtamames/SqueezeMeta/blob/v1.7.2/lib/SQMtools/R/loadSQM.R#L183-L193) non include
   `load_bins`; `load_sequences = FALSE` riduce solo le sequenze. Non esiste
   quindi il workaround ufficiale `load_bins = FALSE`.

Se servono solo profili tassonomici/funzionali già prodotti, [`loadSQMlite()`](https://github.com/jtamames/SqueezeMeta/blob/v1.7.2/lib/SQMtools/R/loadSQMlite.R#L1-L6)
legge `results/tables` e non contiene dettagli di ORF, contig o bin: è un
workaround analitico, non una riparazione dei bin.
