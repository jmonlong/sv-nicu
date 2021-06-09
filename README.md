Annotation pipeline for structural variants found in long reads (e.g. ONT). 
Uses public SV catalogs, repeat, gene, regulatory tracks to annotate and explore the SVs in an HTML report (and paired TSV file).

Jump to:
- [Docker container](#docker-container)
- [Running the pipeline](#running-the-pipeline)
- [Inputs](#inputs)
- [Outputs](#outputs)
    - [Column names](#column-names)
    - [TSV file](#tsv-file)
- [Methods](#methods)
    - [Tiers definition](#tiers)
- [Recommendations for interpretation](#recommendations-for-interpretation)

## Docker container

`quay.io/jmonlong/nicu`, see [quay.io page](https://quay.io/repository/jmonlong/svnicu)

For now, I build and push the container [manually](buil_docker.sh). 
Once I've put the necessary annotation files somewhere (or made them inputs), it will be automatically built.

## Running the pipeline

Two modes depending on the input calls: single VCF or multiple VCFs (to be concatenated, e.g. one per chromosome).
Within the associated Docker container (see [Dockerfile](Dockerfile)):

```sh
## single-VCF input
snakemake --snakefile /scripts/Snakefile --config vcf=sniffles.vcf.gz bam=reads.bam ref=hg37.fa html=sv-report.html tsv=sv-annotated.tsv karyo_tsv=chr-arm-karyotype.tsv --cores 2

## VCF list input
snakemake --snakefile /scripts/Snakefile --config vcf_list=vcf.list.txt bam=reads.bam ref=hg37.fa html=sv-report.html tsv=sv-annotated.tsv karyo_tsv=chr-arm-karyotype.tsv --cores 2
```

It's also possible to force the sample name in the merged VCF using `sample=`.

## Inputs

- Aligned reads, sorted and indexed. Either:
   - One BAM file (`bam=`)
   - A TSV file with 2 columns (chromosome + path to BAM for this chromosome) (`bam_list=`)
- SV calls from Sniffles in VCF. Either
   - One VCF (`vcf=`)
   - Multiple VCFs that will be merged (`vcf_list=`)
- Reference file in indexed FASTA (`ref=`)

## Outputs

The pipeline creates three outputs that contain almost exactly the same information: a HTML report (`html=`) and two TSV files containing the annotated SVs (`tsv=`) and synthetic karyotype at chr-arm level (`karyo_tsv=`).
Examples are available in the [*examples* folder](examples), e.g. the [HTML report](examples/sv-report.html).

### Column names

In the HTML report and [TSV file with annotated SVs](examples/sv-annotated.tsv):

- `pLI` prob of  loss-of-function intolerance described above (`-1` if no information available).
- `type` SV type
- `ac` allele count: `1` for het or `2` for hom.
- `qual` the call quality. Corresponds to the read support (RE field in the Sniffles' VCF)
- `freq` allele frequency of similar SV in gnomAD-SV (>10K genomes sequenced with Illumina whole-genome sequencing).
- `hpgp` frequency in our 11 in-house nanopore controls.
- `dgv` does the variant overlap any variant in DGV in any way?
- `clinsv` any similar clinical SV variants (nstd102 in dbVar)? "Similar" defined as reciprocal overlap > 50%.
- `ctcf` does the variant overlap a CTCF binding site? From ENCODE track for kidney.
- `cres` does the variant overlap a regulatory region? From ENCODE track for kidney.
- `simp.rep` is the variant in or close to a simple repeat (see simple repeat track in the UCSC Genome Browser).
- `cons` does the variant overlap a conserved element (as defined by 100 vertebrate phastCons track)
- `impact` potential impact based on gene annotation: *coding*, *UTR*, *promoter*, *intronic*.
- `genes`/`gene` information about the gene(s) overlapped by the variant and their impact.
- `gene.dist` distance to the nearest gene which is specified by `gene.near`.
- `sel.gene.dist` distance to the nearest gene of interest which is specified by `sel.gene.near`.
- `cds.dist` distance to the nearest coding regions, useful for intronic variants for example.
- `nb.pc.genes` number of protein-coding genes overlapped by the variant.
- `cov` median scaled coverage for large variants (>100 kbp), if *mosdepth* results are available. Between parenthesis is the number of bin overlapping the variant (the higher the more confident).
- `large` a boolean value flagging large SVs (>100 kbp). Because they often overlap hundreds of genes and are potentially false-positives, it is convenient to remove them from the table using this column. There is a tab specifically about large SVs that is more appropriate to investigate those.

In the [synthetic karyotype TSV file](examples/chr-arm-karyotype.tsv):

- `chr` and `arm` chromosome name and either *p*/*q* for the arm.
- `median.cov` the median read coverage, normalized i.e. *1* means diploid.
- `aberrant` a flag for chr-arm median coverage is closer to an aneuploid state (<0.75 or >1.25).


### Annotated SV TSV file

The TSV file with annotated SV is actually gene-centric and contains all the annotations described above.
Gene-centric means that there is one row for each gene-variant pair.
This helps filter on gene features (e.g. gene of interest, pLI) although large SVs spanning many genes are "duplicated" in the file.


## Methods

### SV calling 

The structural variants were called from the nanopore reads using [Sniffles](https://github.com/fritzsedlazeck/Sniffles).
[mosdepth](https://github.com/brentp/mosdepth) provides quick estimates of read coverage to confirm large CNVs or identify chromosomal aberrations. 

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
To provide an additional list, format gene names into a file in either of the following formats:

1. one gene per line. This list will be called "*custom*" in the output tsv/report
2. a tab-separated file with one column for gene name and one column for list name (no headers).

To use this gene list, add `gene_list=FILE` to the `--config` options.

### Filters

In most tables, we removed common variants, i.e. either:

- frequency higher than 1% in gnomAD-SV
- seen in the SV catalog from long-read studies

### Tiers

- Tier 1: Deletions + Insertions in Exons; SVs >100kbp overlapping genes and with AF<10%; chr or chr-arm aberration (e.g. aneuploidy)
- Tier 2: INS + DEL + tandemDUP in genes or promoter (introns, UTR, promoter).
- Tier 3: Genome wide range. The genome wide scope of SV including rearrangements (Inversions, Translocations). These events will be hard to interpret but maybe useful to track over time. 

### Database object

The known SV and other annotation tracks are prepared into one file using [prepare-sv-report-data.R](prepare-sv-report-data.R). 
Note that this script requires SV calls on our 11 control genomes, not yet available in this repo.

### Local re-assembly of tier 1 SVs

Supporting reads identified by Sniffles are extracted and assembled using [shasta](https://github.com/chanzuckerberg/shasta). 
The assembled sequenced is then aligned to the reference genome with [minimap2](https://github.com/lh3/minimap2) and the (potentially) fine-tuned SV identified using [SVIM-asm](https://github.com/eldariont/svim-asm).
If this re-assembled SV is similar to the original Sniffles call, it is used to update its breakpoint definition.

This module is enabled by default.
To disable it, add `amb=false` to the `--config` options.

## Recommendations for interpretation

- Large CNVs: use the `cov` column as an orthogonal support. If coverage around 1, it's likely a false positive.
- `_asm` in the SV id (`svid` column) means the allele could be re-assembled from raw reads. It's a good sign.
- `qual` column can be used as confidence, the higher the better.
- Take the genotype information (`ac` column for allele count) with a grain of salt.
- Using the HTML report
    - Click on the links to jump to UCSC Genome Browser (`coord` column) or gnomAD (`pLI` column) and other gene information (`gene` column).
    - Remove (FP?) large CNVs from the tables using the `large` column. 
