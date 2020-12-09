## config
## inputs:
## - 'vcf' VCF from Sniffles OR `vcf_list` text file with list of VCFs to concatenate
## - 'bam' BAM header
## - 'bai' BAM index
## - 'ref' reference FASTA
## output
## - 'html' HTML report
## - 'tsv' TSV file
if 'html' not in config:
    config['html'] = 'sv-report.html'
if 'tsv' not in config:
    config['tsv'] = 'sv-annotated.tsv'
if 'karyo_tsv' not in config:
    config['karyo_tsv'] = 'chr-arm-karyotype.tsv'

## sample name
SAMPLE = 'FAST0X'
if 'sample' in config:
    SAMPLE = config['sample']

## input VCF or VCF list
VCF = ''
VCFS = []
if 'vcf' in config:
    VCF = config['vcf']
if VCF == '' and 'vcf_list' in config:
    inf = open(config['vcf_list'], 'r')
    for line in inf:
        VCFS.append(line.rstrip().rstrip('.vcf'))
    inf.close()
    VCF = SAMPLE + '.sniffles.merged.vcf'

## additional gene list?
GENE_LIST = 'NA'
if 'gene_list' in config:
    GENE_LIST = config['gene_list']

## optional: scripts and database for annotation
if 'rmd' not in config:
    config['rmd'] = '/scripts/sv-report.Rmd'
if 'db' not in config:
    config['db'] = '/scripts/sv_annotation_database.RData'
if 'extract_tier1_r' not in config:
    config['extract_tier1_r'] = '/scripts/extract_tier1_svs.R'
if 'extract_reads_py' not in config:
    config['extract_reads_py'] = '/scripts/extract_reads.py'
if 'merge_contigs_py' not in config:
    config['merge_contigs_py'] = '/scripts/merge_contigs.py'

EXTRACT_TIER1_R = config['extract_tier1_r']
EXTRACT_READS_PY = config['extract_reads_py']
MERGE_CONTIGS_PY = config['merge_contigs_py']

rule make_report:
    input:
        rmd=config['rmd'],
        vcf_t1='svs_tier1.tsv',
        vcf_t1_ft='svs_tier1_finetuned.vcf',
        vcf_other='svs_other.tsv',
        idxcov='indexcov_res/indexcov_res-indexcov.bed.gz',
        sv_db=config['db']
    output:
        html=config['html'],
        tsv=config['tsv'],
        karyo=config['karyo_tsv']
    shell:
        """
        cp {input.rmd} report.Rmd
        Rscript -e 'rmarkdown::render("report.Rmd", output_file="{output.html}")' {input.vcf_t1}  {input.vcf_t1_ft}  {input.vcf_other}  {input.idxcov}  {input.sv_db} {GENE_LIST} {output.tsv} {output.karyo}
        """

rule rehead_vcf:
    input: '{vcf}.vcf'
    output: '{vcf}.rn.vcf'
    params: temp='{vcf}.temp.txt'
    shell:
        """
        echo {SAMPLE} > {params.temp}
        bcftools reheader -s {params.temp} {input} > {output}
        """

rule merge_vcfs:
    input: expand('{vcf}.rn.vcf', vcf=VCFS)
    output: VCF
    shell: "bcftools concat {input} > {output}"

## extract tier 1 SVs to be finetuned or visualized
rule extract_tiers:
    input:
        vcf=VCF,
        sv_db=config['db']
    output: 
        vcf_t1='svs_tier1.tsv',
        vcf_other='svs_other.tsv'
    script: EXTRACT_TIER1_R

if 'amb' in config and not config['amb']:
    rule skip_amb:
        input: 'svs_tier1.tsv'
        output: 'svs_tier1_finetuned.vcf'
        shell: 'touch {output}'
else:
    ## extract reads supporting tier1 svs, to be used for assembly 
    checkpoint extract_reads:
        input:
            svs='svs_tier1.tsv',
            bam=config['bam'],
            bai=config['bai']
        output: "svs_tier1_reads.tsv"
        script: EXTRACT_READS_PY

    ## assemble reads supporting a SV
    rule asm_shasta:
        input: '{svid}.reads.fa'
        output: '{svid}.contigs.fa'
        threads: 2
        params:
            odir= '{svid}_shasta_out'
        shell:
            """
            rm -rf {params.odir}
            shasta --input {input} --config /scripts/Nanopore-Sep2020.conf --assemblyDirectory {params.odir} --threads {threads}
            cp {params.odir}/Assembly.fasta {output}
            """

    ## merge assembled contigs for all svs
    def aggregate_contigs(wildcards):
        reads_tsv = checkpoints.extract_reads.get(**wildcards).output[0]
        res = []
        with open(reads_tsv, 'r') as inf:
            heads = next(inf)
            for line in inf:
                line = line.rstrip().split('\t')
                contig_fa = line[3].replace('.reads.fa', '.contigs.fa')
                res.append(contig_fa)
        return res
    rule merge_contigs:
        input: aggregate_contigs
        output: "all.fa"
        script: MERGE_CONTIGS_PY

    ## align contigs to reference
    rule map_contigs:
        input:
            asm="all.fa",
            ref=config['ref']
        output: "all.sam"
        shell: "minimap2 --paf-no-hit -a -x asm5 --cs -r2k -t 2 {input.ref} {input.asm} > {output}"

    ## sort aligned contigs
    rule sort_bam:
        input: "{reads}.sam"
        output:
            bam="{reads}.sorted.bam",
            bai="{reads}.sorted.bam.bai"
        shell:
            """
            samtools sort -o {output.bam} {input}
            samtools index {output.bam}
            """

    ## call SVs from aligned contigs
    rule call_contigs:
        input:
            bam="all.sorted.bam",
            bai="all.sorted.bam.bai",
            ref=config['ref']
        output: 'svs_tier1_finetuned.vcf'
        params: odir='svimasm_out'
        shell:
            """
            svim-asm haploid {params.odir} {input.bam} {input.ref}
            cp {params.odir}/variants.vcf {output}
            """

## quick read coverage track from the BAM index
rule indexcov:
    input:
        bam=config['bam'],
        bai=config['bai']
    output: 'indexcov_res/indexcov_res-indexcov.bed.gz'
    params:
        dir='indexcov_res/',
        bam='temp.bam',
        bai='temp.bai'
    shell:
        """
        rm -f {params.bam} {params.bai}
        ln -s {input.bam} {params.bam}
        ln -s {input.bai} {params.bai}
        goleft indexcov --directory {params.dir} --sex X,Y {params.bam}
        """

