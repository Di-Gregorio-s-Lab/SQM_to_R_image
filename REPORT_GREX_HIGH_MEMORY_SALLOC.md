# Grex: nodi Slurm con almeno 200 GB RAM

Fonte unica: [documentazione ufficiale Grex — Slurm partitions](https://um-grex.github.io/docs/running-jobs/slurm-partitions/) e [interactive jobs](https://um-grex.github.io/docs/running-jobs/interactive-jobs/), consultate l'11 settembre 2026.

## Partizioni/nodi idonei

La tabella riporta la RAM **per nodo** esattamente come dichiarata da Grex. Per una richiesta di almeno 200 GB, chiedere un singolo nodo della partizione indicata e non superarne il limite.

| Classe Grex | Partizione | Nodi | CPU/GPU per nodo | RAM/nodo |
|---|---:|---:|---|---:|
| CPU generale | `largemem` | 12 | 40 CPU, Cascade Lake | 380 Gb |
| CPU generale | `genoa` | 27 | 192 CPU, AMD EPYC 9654 | 750 Gb |
| CPU generale | `genlm` | 3 | 192 CPU, AMD EPYC 9654 | 1500 Gb |
| CPU generale (speciale) | `test` | 1 | 18 CPU, Cascade Lake | 512 Gb |
| GPU generale | `lgpu` | 2 | 2× L40s/48 GB; 64 CPU | 380 Gb |
| CPU contribuita | `mcordcpu` | 5 | 168 CPU, AMD EPYC 9634 | 1500 Gb |
| CPU contribuita | `chrim` | 4 | 192 CPU, AMD EPYC 9654 | 750 Gb |
| CPU contribuita | `chrimlm` | 1 | 192 CPU, AMD EPYC 9654 | 1500 Gb |
| CPU contribuita | `hsc` | 1 | 192 CPU, AMD EPYC 9654 | 1500 Gb |
| CPU contribuita | `pgs` | 1 | 192 CPU, AMD EPYC 9655 | 750 Gb |
| GPU contribuita | `livi` | 1 | 16× V100/32 GB; 48 CPU | 1500 Gb |
| GPU contribuita | `agro` | 2 | 2× A30/24 GB; 24 CPU | 250 Gb |
| GPU contribuita | `mcordgpu` | 2 | 4× A30/24 GB; 32 CPU | 512 Gb |

Le risorse provengono dalle tabelle ufficiali CPU/GPU generali e contribuite. `hsc` è elencata nella tabella “Contributed CPU” con 1500 Gb; la relativa nota a piè pagina la chiama però “GPU node”, quindi verificare con `sinfo` prima di richiedere GPU.

## Vincoli di accesso e scheduler

- Grex richiede in molti casi `--partition=`: non seleziona automaticamente una partizione compatibile, e un job non può eseguire su più partizioni contemporaneamente. Si possono elencare più partizioni, ad esempio `--partition=skylake,largemem`.
- La documentazione non prescrive un flag Slurm `--constraint` per questi nodi: l'instradamento documentato è per `--partition=`. I tipi CPU/GPU nella tabella sono le caratteristiche hardware, non nomi di constraint da inventare.
- `test` consente oversubscription ed è pensata anche per turnaround di job interattivi/OOD: non è la scelta per una riserva esclusiva affidabile da 200 GB.
- Una partizione GPU richiede anche una richiesta GPU (`--gpus=N`, `--gpus-per-node=N` o `--gpus-per-task=N`); un job solo CPU viene rifiutato.
- Le partizioni contribuite sono riservate ai rispettivi gruppi proprietari. Per utenti esterni, Grex espone le alternative preemptible `livi-b`, `agro-b`, `mcordgpu-b` e `genoacpu-b` (quest'ultima copre i nodi CPU AMD Genoa contribuiti); il minimo garantito dichiarato è un'ora e il job può poi essere terminato per preemption.

## `salloc` ufficiale e richiesta da 200 GB

Grex prescrive `salloc` per un job interattivo e richiede le risorse come opzioni a riga di comando. La sua sintassi-esempio ufficiale è:

```bash
salloc --nodes=1 --ntasks=1 --cpus-per-task=6 --mem=12000M --partition=skylake --time=0-2:00:00
```

Per una sessione CPU interattiva da 200 GB su una partizione generale idonea, l'applicazione diretta della stessa sintassi è:

```bash
salloc --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=200G --partition=largemem --time=0-2:00:00
```

Sostituire `largemem` solo con una partizione a cui l'account è autorizzato; scegliere `--cpus-per-task` e `--time` in base al programma. Grex consiglia per l'interattivo richieste minime (in genere meno di 3 ore e meno di 4 GB per core) per ridurre l'attesa: 200 GB è quindi una richiesta grande e può restare in coda più a lungo.
