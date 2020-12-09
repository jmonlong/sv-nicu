## Load packages
library(sveval)
library(GenomicRanges)
library(dplyr)

message('gnomAD-SV...')
if(!file.exists('gnomad_v2.1_sv.controls_only.sites.vcf.gz')){
  download.file('https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.controls_only.sites.vcf.gz', 'gnomad_v2.1_sv.controls_only.sites.vcf.gz')
}
gnomad = readSVvcf('gnomad_v2.1_sv.controls_only.sites.vcf.gz', other.field='AF')
## for CNVs, for now let's just remove CN2 and mark them as common (AF=1)
gnomad = subset(gnomad, alt != '<CN=2>' & type!='BND' & type!='CPX')
gnomad$AF = ifelse(gnomad$alt %in% c('<CN=0>', '<CN=1>'), 1, gnomad$AF)
gnomad$type = ifelse(gnomad$alt %in% c('<CN=0>', '<CN=1>'), 'DEL', gnomad$type)
gnomad$AF = ifelse(gnomad$type == 'MCNV', 1, gnomad$AF)
gnomad$type = ifelse(gnomad$type == 'MCNV', 'DUP', gnomad$type)
gnomad$AF = as.numeric(gnomad$AF)
gnomad$ac = 1
gnomad$alt = NULL
gnomad = unique(gnomad)

message('GIAB long read SV catalog')
if(!file.exists('HG002_SVs_Tier1_v0.6.vcf.gz')){
  download.file('ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz', 'HG002_SVs_Tier1_v0.6.vcf.gz')
}
giab = readSVvcf('HG002_SVs_Tier1_v0.6.vcf.gz', sample=NULL)
giab$ac = 1

message('In-house nanopore samples')
vcfs = list.files('.', '.sniffles.vcf.gz')
## vcfs = setdiff(vcfs, paste0(c('HG002', 'HG003', 'HG004'), '.sniffles.vcf.gz')) ## Temporary to test
hpgp = lapply(vcfs, function(vcff){
  svs = readSVvcf(vcff, qual.field='RE')
  svs$sample = gsub('.sniffles.vcf.gz', '', vcff)
  svs
})
hpgp = do.call(c, hpgp)
hpgp = subset(hpgp, ac>0 & seqnames(hpgp) %in% c(1:22,'X','Y','MT'))
hpgp$ac = 1
hpgp$qual = hpgp$ref = hpgp$alt = NULL

message('Gene annotation...')
if(!file.exists('gencode.v35lift37.annotation.gtf.gz')){
  download.file('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh37_mapping/gencode.v35lift37.annotation.gtf.gz', 'gencode.v35lift37.annotation.gtf.gz')
}
types.ranked = c('CDS', 'UTR', 'promoter', 'gene')
types.labels = c('coding', 'UTR', 'promoter', 'intronic')
genc = rtracklayer::import('gencode.v35lift37.annotation.gtf.gz')
genc = subset(genc, type %in% types.ranked)
prom = promoters(subset(genc, type=='gene'))
prom$type = 'promoter'
genc = c(genc, prom)
mcols(genc) = mcols(genc)[,c('type', 'gene_name', 'gene_type')]
seqlevels(genc) = gsub('chr', '', seqlevels(genc))

message('Clinical SVs')
if(!file.exists('nstd102.GRCh37.variant_call.tsv.gz')){
  download.file('https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/tsv/nstd102.GRCh37.variant_call.tsv.gz', 'nstd102.GRCh37.variant_call.tsv.gz')
}
clinsv = read.table('nstd102.GRCh37.variant_call.tsv.gz', header=TRUE, as.is=TRUE, sep='\t', skip=1, comment='', quote='')
gain.types = c('copy number gain', 'duplication', 'insertion', 'tandem duplication')
loss.types = c('copy number loss', 'deletion')
clinsv = clinsv %>% filter(clinical_significance=='Pathogenic',
                           variant_call_type %in% c(gain.types, loss.types),
                           chr %in% c(1:22, 'X','Y')) %>%
  mutate(type=ifelse(variant_call_type %in% gain.types, 'INS', 'DEL'),
         start=ifelse(is.na(start), ifelse(is.na(inner_start),
                                           outer_start, inner_start), start),
         end=ifelse(is.na(stop), ifelse(is.na(inner_stop),
                                        outer_stop, inner_stop), stop),
         size=end-start) %>%
  select(chr, start, end, size, type) %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE)

message('Gene list')
gene.list = read.table('nicu-gene-list.tsv', sep='\t', as.is=TRUE, header=TRUE)

message('Simple repeats')
if(!file.exists('simpleRepeat.txt.gz')){
  download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz', 'simpleRepeat.txt.gz')
}
sr = read.table('simpleRepeat.txt.gz')
sr = sr[,c(2,3,4)]
colnames(sr) = c('chrom', 'start', 'end')
sr$chrom = gsub('chr', '', sr$chrom)
sr = sr %>% filter(chrom %in% c(1:22, 'X', 'Y')) %>% makeGRangesFromDataFrame(keep.extra.columns=TRUE)
sr = reduce(sr)

message('DGV catalog')
if(!file.exists('GRCh37_hg19_variants_2020-02-25.txt')){
  download.file('http://dgv.tcag.ca/dgv/docs/GRCh37_hg19_variants_2020-02-25.txt', 'GRCh37_hg19_variants_2020-02-25.txt')
}
dgv = read.table('GRCh37_hg19_variants_2020-02-25.txt', as.is=TRUE, header=TRUE, sep='\t')
dgv = dgv[,c('chr','start','end', 'variantsubtype')]
dgv = makeGRangesFromDataFrame(dgv, keep.extra.columns=TRUE)
dgv.loss = subset(dgv, variantsubtype %in% c('loss', 'deletion', 'gain+loss'))
dgv.gain = subset(dgv, variantsubtype %in% c('duplication', 'gain', 'insertion',
                                             'mobile element insertion',
                                             'novel sequence insertion',
                                             'tandem duplication', 'gain+loss'))


message('CTCF peaks')
## ctcf.df = read.table('./encode-ctcf.tsv', as.is=TRUE, header=TRUE, sep='\t', comment.char='', skip=1, quote='')
## unique(ctcf.df$Biosample.term.name)
## ctcf = lapply(c("thyroid gland", "heart left ventricle", "lung", "liver", "tibial artery",
##                 "cardiac muscle cell", "astrocyte of the cerebellum", "astrocyte",
##                 "kidney", "neural progenitor cell", "neural cell", "pancreas"), function(biosamp){
## ff = ctcf.df %>% filter(Biosample.term.name==biosamp) %>% head(1) %>% 
if(!file.exists('ENCFF415WKV.bed.gz')){
  download.file('https://www.encodeproject.org/files/ENCFF415WKV/@@download/ENCFF415WKV.bed.gz', 'ENCFF415WKV.bed.gz')
}
ctcf = read.table('ENCFF415WKV.bed.gz', as.is=TRUE)
ctcf = with(ctcf, GRanges(gsub('chr', '', V1), IRanges(V2, V3), score=V5))
ctcf = reduce(ctcf)

message('Regulatory regions')
if(!file.exists('ENCFF509DLH.bed.gz')){
  download.file('https://www.encodeproject.org/files/ENCFF509DLH/@@download/ENCFF509DLH.bed.gz', 'ENCFF509DLH.bed.gz')
}
cres = read.table('ENCFF509DLH.bed.gz', as.is=TRUE)
cres = with(cres, GRanges(gsub('chr', '', V1), IRanges(V2, V3), score=V5))
cres = reduce(cres)

message('Conserved regions')
if(!file.exists('phastConsElements100way.txt.gz')){
  download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/phastConsElements100way.txt.gz', 'phastConsElements100way.txt.gz')
}
cons.gr = read.table('phastConsElements100way.txt.gz', as.is=TRUE)
cons.gr = with(cons.gr, GRanges(gsub('chr', '', V2), IRanges(V3, V4)))
cons.gr = reduce(cons.gr)


message('pli scores')
if(!file.exists('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')){
  download.file('https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', 'gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz')
}
pli.df = read.table('gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz', as.is=TRUE, header=TRUE, sep='\t')
pli.df = pli.df %>% select(gene, pLI) %>% unique


message('Cytoband...')
if(!file.exists('cytoBandIdeo.txt.gz')){
  download.file('https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBandIdeo.txt.gz', 'cytoBandIdeo.txt.gz')
}
cyto.df = read.table('cytoBandIdeo.txt.gz', as.is=TRUE, sep='\t')
colnames(cyto.df) = c('chr', 'start', 'end', 'band', 'gieaStain')
cyto.df = cyto.df %>% mutate(chr=gsub('chr', '', chr), arm=substr(band, 1, 1)) %>%
  filter(arm %in% c('p', 'q')) %>% 
  group_by(chr, arm) %>% summarize(start=min(start), end=max(end))
cyto.gr = makeGRangesFromDataFrame(cyto.df, keep.extra.columns=TRUE)


message('Save in RData file...')
save(genc, types.ranked, types.labels, clinsv, gnomad, giab, hpgp, gene.list,
     sr, dgv.loss, dgv.gain, ctcf, cres, cons.gr, pli.df, cyto.gr, file='sv_annotation_database.RData')
