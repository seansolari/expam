## Batching metagenomes and removing contaminated reads

### Batching samples

As we have a large collection of metagenomes to run, create batches of (paired) sequence files with each batch having ~70Gb of sequence data.

``` bash
scripts/preprocessing/batcher.py --chunk 70 -d sample_batches/
```

### Filtering contaminated reads

Run [bowtie2](https://github.com/BenLangmead/bowtie2) to align metagenomic reads against human chromosomes and remove any unambiguously mapping (contaminated) reads.

``` bash
scripts/preprocessing/run_bowtie2.py -o clean_reads/ --batch batch[BATCH ID].txt
```

### Running *expam*

These clean reads can then be used as input to `expam`.

```bash
expam_limit -x 0.7 -t 5.0 -o batchID.txt expam classify -db pibd/ --out ./expam_results/batchID -d clean_reads/ --paired --alpha 0.1 > screen_log_batchID.txt
```

**Next step...**

[Read processing and cleaning metagenomics output](./2_preprocessing.md)