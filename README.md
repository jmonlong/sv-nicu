Annotation pipeline for structural variants found in long reads (e.g. ONT). 
Uses public SV catalogs, repeat, gene, regulatory tracks to annotate and explore the SVs in an HTML report (and paired TSV file).

## Running the pipeline

Within the associated Docker container (see [Dockerfile](Dockerfile)):

```sh
snakemake --snakefile /scripts/Snakefile --config vcf=sniffles.vcf.gz bam=reads.bam bai=reads.bam.bai html=sv-report.html tsv=sv-annotated.tsv --cores 2
```

## Inputs

- Aligned reads, sorted and indexed (`bam=` and `bai=`)
- SV calls from Sniffles in VCF (`vcf=`)
- Reference file in indexed FASTA (`ref=`)

## Outputs

The pipeline creates two outputs that contain almost exactly the same information: a HTML report (`html=`) and a TSV file (`tsv=`).
Both uses the following column names

Examples are available in the [*examples* folder](examples), e.g. the [HTML report](examples/sv-report.html).

### Column names

- `pLI` prob of  loss-of-function intolerance described above (`-1` if no information available).
- `type` SV type
- `ac` allele count: `1` for het or `2` for hom.
- `qual` the call quality. Corresponds to the read support (RE field in the Sniffles' VCF)
- `freq` allele frequency of similar SV in gnomAD-SV (>10K genomes sequenced with Illumina whole-genome sequencing).
- `dgv` does the variant overlap any variant in DGV in any way?
- `clinsv` any similar clinical SV variants (nstd102 in dbVar)? "Similar" defined as reciprocal overlap > 50%.
- `ctcf` does the variant overlap a CTCF binding site? From ENCODE track for kidney.
- `cres` does the variant overlap a regulatory region? From ENCODE track for kidney.
- `simp.rep` is the variant in or close to a simple repeat (see simple repeat track in the UCSC Genome Browser).
- `cons` does the variant overlap a conserved element (as defined by 100 vertebrate phastCons track)
- `impact` potential impact based on gene annotation: *coding*, *UTR*, *promoter*, *intronic*.
- `genes`/`gene`/`gene_type` information about the gene(s) overlapped by the variant and their impact.
- `gene.dist` distance to the nearest gene which is specified by `gene.near`.
- `sel.gene.dist` distance to the nearest gene of interest which is specified by `sel.gene.near`.
- `cds.dist` distance to the nearest coding regions, useful for intronic variants for example.
- `nb.pc.genes` number of protein-coding genes overlapped by the variant.
- `cov` median scaled coverage for large variants (>100 kbp), if *indexcov* results are available. Between parenthesis is the number of bin overlapping the variant (the higher the more confident).
- `large` a boolean value flagging large SVs (>100 kbp). Because they often overlap hundreds of genes and are potentially false-positives, it is convenient to remove them from the table using this column. There is a tab specifically about large SVs that is more appropriate to investigate those.

### TSV file

In addition to the HTML report, a gene-centric TSV file is written with all the annotation described above.
Gene-centric means that there is one row for each gene-variant pair.
This helps filter on gene features (e.g. gene of interest, pLI).


## Methods

### SV calling 

The structural variants were called from the nanopore reads using [Sniffles](https://github.com/fritzsedlazeck/Sniffles).
[indexcov](https://github.com/brentp/goleft/tree/master/indexcov) provides quick estimates of read coverage to confirm large CNVs or identify chromosomal aberrations. 

To select for higher confidence SVs, increase the minimum quality (corresponding to the read support, RE in Sniffles' VCF). 
In the tables, the `qual` column has been winsorized at 100 to help selecting the practical ranges.

In the report, we don't consider variation in alternate chromosomes or the mitochondrial genome.

### SV database and frequency annotation

The variants were compared to catalogs of known SVs.
The frequency estimates are based on the [gnoma-SV catalog](https://macarthurlab.org/2019/03/20/structural-variants-in-gnomad/), as the maximum frequency of variants with reciprocal overlap >10%.
We filter out calls that match the GIAB ([Zook et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32541955/)) public catalog and 11 in-house control genomes (nanopore + Sniffles) to remove variants that pass the frequency filter based on gnomAD-SV simply because they can't be detected by short reads.
Finally, variants are flagged if overlapping [DGV](http://dgv.tcag.ca/dgv/app/home) (any overlap) or [Clinival SVs (nstd102 dbVar)](https://www.ncbi.nlm.nih.gov/dbvar/studies/nstd102/) (reciprocal overlap>50%).

### Gene annotation

The [GENCODE](https://www.gencodegenes.org/) gene annotation was used to flag variants as *coding*, *UTR*, *promoter*, *intronic* (prioritized in this order) and compute the distance to the nearest gene.
While we consider lncRNA, miRNA in the annotation (`genes` column), most tables focus on protein-coding genes.

The *pLI* score was computed by the [gnomAD project](https://gnomad.broadinstitute.org/).
It represents the probability that the gene is intolerant to loss-of-function variants.
A variant affecting a gene with a high pLI score (e.g. >0.9) is more likely to have a biological impact.
The variants in each section are ordered to show those affecting genes with the highest pLI first.

### Genes of interest

The report uses a few gene lists by default.
To provide an additional list, format gene names into a file (one gene per line) and add `gene_list=FILE` to the `--config` options.

### Filters

In most tables, we removed common variants, i.e. either:

- frequency higher than 1% in gnomAD-SV
- seen in the SV catalog from long-read studies

### Tiers

- Tier 1: Deletions + Insertions in Exons and aneuploidy. 
- Tier 2: INS + DEL+ tandemDUP in genes (Introns + Exons).
- Tier 3: Genome wide range. The genome wide scope of SV including rearrangements (Inversions, Translocations). These events will be hard to interpret but maybe useful to track over time. 

### Database object

The known SV and other annotation tracks are prepared into one file using [prepare-sv-report-data.R](prepare-sv-report-data.R). 
Note that this script requires SV calls on our 11 control genomes, not yet available in this repo.
