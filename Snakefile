## config
## inputs:
## - 'vcf' VCF from Sniffles OR `vcf_list` text file with list of VCFs to concatenate
## - 'bam_list' BAM list. TSV file with 1st column is chr name, 2nd column is path to BAM file
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

## input BAM list: chr -> path
BAMS = {}
if 'bam_list' in config:
    inf = open(config['bam_list'], 'r')
    for line in inf:
        line = line.rstrip().split('\t')
        BAMS[line[0]] = line[1]
    inf.close()
def input_chr_bam(wildcards):
    return BAMS[wildcards.chr]

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

if 'ref' not in config:
    config['ref'] = 'ref.fa'

## additional gene list?
GENE_LIST = 'NA'
if 'gene_list' in config:
    GENE_LIST = config['gene_list']

# should the tier1 SVs be re-assembled locally? (default is true)
if 'amb' in config:
    if type(config['amb']) == str:
        config['amb'] = config['amb'] in ['T', 'true', 'True']
else:
    config['amb'] = True

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

BIN=10000

rule make_report:
    input:
        rmd=config['rmd'],
        vcf_t1='svs_tier1.tsv',
        vcf_t1_ft='svs_tier1_finetuned.vcf',
        vcf_other='svs_other.tsv',
        cov=expand('mosdepth.bin{bin}bp.regions.bed.gz', bin=BIN),
        sv_db=config['db']
    output:
        html=config['html'],
        tsv=config['tsv'],
        karyo=config['karyo_tsv']
    benchmark: 'benchmark_sv_annotation/make_report.benchmark.tsv'
    shell:
        """
        cp {input.rmd} report.Rmd
        Rscript -e 'rmarkdown::render("report.Rmd", output_file="{output.html}")' {input.vcf_t1}  {input.vcf_t1_ft}  {input.vcf_other}  {input.cov}  {input.sv_db} {GENE_LIST} {output.tsv} {output.karyo}
        """

rule rehead_vcf:
    input: '{vcf}.vcf'
    output: '{vcf}.rn.vcf'
    params: temp='{vcf}.temp.txt'
    shell:
        """
        echo {SAMPLE} > {params.temp}
        bcftools reheader -s {params.temp} {input} > {output}
        rm {params.temp}
        """

rule merge_vcfs:
    input: expand('{vcf}.rn.vcf', vcf=VCFS)
    output: VCF
    shell:
        """
        grep "#" {input[0]} > {output}
        cat {input} | grep -v '#' >> {output}
        """

## extract tier 1 SVs to be finetuned or visualized
rule extract_tiers:
    input:
        vcf=VCF,
        sv_db=config['db']
    output: 
        vcf_t1='svs_tier1.tsv',
        vcf_other='svs_other.tsv'
    script: EXTRACT_TIER1_R

if not config['amb']:
    rule skip_amb:
        input: 'svs_tier1.tsv'
        output: 'svs_tier1_finetuned.vcf'
        shell: 'touch {output}'
else:
    ## extract reads supporting tier1 svs, to be used for assembly 
    rule extract_reads_chr:
        input:
            svs='svs_tier1.tsv',
            bam=input_chr_bam
        output: "svs_tier1_reads_{chr}.tsv"
        benchmark: 'benchmark_sv_annotation/extract_reads.{chr}.benchmark.tsv'
        shell: "python3 {EXTRACT_READS_PY} -b {input.bam} -c {wildcards.chr} -v {input.svs} -o {output}"

    checkpoint extract_reads:
        input: expand('svs_tier1_reads_{chr}.tsv', chr=BAMS.keys())
        output: "svs_tier1_reads.tsv"
        benchmark: 'benchmark_sv_annotation/extract_reads.benchmark.tsv'
        run:
            outf = open(output[0], 'w')
            inf = open(input[0], 'r')
            outf.write(next(inf))
            inf.close()
            for ff in input:
                inf = open(ff, 'r')
                header = next(inf)
                for line in inf:
                    outf.write(line)
                inf.close()
            outf.close()

    ## assemble reads supporting a SV
    rule asm_shasta:
        input: '{svid}.reads.fa'
        output: '{svid}.contigs.fa'
        threads: 2
        params:
            odir= '{svid}_shasta_out'
        benchmark: 'benchmark_sv_annotation/asm_shasta.{svid}.benchmark.tsv'
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
        benchmark: 'benchmark_sv_annotation/map_contigs.benchmark.tsv'
        shell: "minimap2 --paf-no-hit -a -x asm5 --cs -r2k -t 2 {input.ref} {input.asm} > {output}"

    ## sort aligned contigs
    rule sort_bam:
        input: "{reads}.sam"
        output:
            bam="{reads}.sorted.bam",
            bai="{reads}.sorted.bam.bai"
        benchmark: 'benchmark_sv_annotation/sort_bam.{reads}.benchmark.tsv'
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
        benchmark: 'benchmark_sv_annotation/call_contigs.benchmark.tsv'
        shell:
            """
            svim-asm haploid {params.odir} {input.bam} {input.ref}
            cp {params.odir}/variants.vcf {output}
            """

rule mosdepth_chr:
    input: input_chr_bam
    output: 'mosdepth.{chr}.bin{bin}bp.regions.bed.gz'
    params:
        prefix='mosdepth.{chr}/mosdepth.{chr}.bin{bin}bp',
        out_dir='mosdepth.{chr}'
    benchmark: 'benchmark_sv_annotation/mosdepth.{chr}.bin{bin}bp.benchmark.tsv'
    shell:
        """
        mkdir -p {params.out_dir}
        mosdepth -b {wildcards.bin} -c {wildcards.chr} -n -x -m {params.prefix} {input}
        cp {params.prefix}.regions.bed.gz {output}
        """

rule merge_mosdepth:
    input: expand('mosdepth.{chr}.bin{{bin}}bp.regions.bed.gz', chr=BAMS.keys())
    output: 'mosdepth.bin{bin}bp.regions.bed.gz'
    shell:
        """
        zcat {input} | gzip > {output}
        """

rule clean:
    shell:
        """
        rm -rf mosdepth.*.bin*bp.regions.bed.gz all.fa all.sam all.sorted.bam all.sorted.bam.bai *.reads.fa *.contigs.fa *_shasta_out svs_tier1_reads_*.tsv *.rn.vcf
        """
    
