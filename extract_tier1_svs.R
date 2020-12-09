## Tier 1 defined as coding SVs with frequency <1% and not in in-house database

library(sveval)
library(GenomicRanges)
library(dplyr)

winsor <- function(x, u){
  if(any(x>u)) x[x>u] = u
  x
}

message('Load database for annotation...')
load(snakemake@input[['sv_db']])

message('Read VCF...')
svs.gr = readSVvcf(snakemake@input[['vcf']], qual.field='RE', keep.ids=TRUE, other.field='RNAMES')
svs.gr = subset(svs.gr, ac>0)
svs.gr$qual = winsor(svs.gr$qual, 100)
svs.gr = svs.gr[as.character(seqnames(svs.gr)) %in% c(1:22, 'X', 'Y')]

message('Annotate frequency using gnomAD...')
ol.df = svOverlap(svs.gr, gnomad, min.ol=.1, max.ins.dist=100) %>% as.data.frame %>%
  mutate(freq=gnomad$AF[subjectHits]) %>%
  group_by(queryHits) %>%
  summarize(freq=max(freq))
svs.gr$freq = 0
svs.gr$freq[ol.df$queryHits] = ol.df$freq

message('Overlap with GIAB SV catalog from long reads')
ol.gr = svOverlap(svs.gr, giab, min.ol=.5, max.ins.dist=100)
svs.gr$giab = FALSE
svs.gr$giab[ol.gr$queryHits] = TRUE

message('Overlap with SVs from in-house  long-read studies...')
hpgp.nsamps = length(unique(hpgp$sample))
ol.df = svOverlap(svs.gr, hpgp, min.ol=.5, max.ins.dist=100) %>%
  as.data.frame %>%
  mutate(sample=hpgp$sample[subjectHits]) %>%
  group_by(queryHits) %>% summarize(freq=length(unique(sample))/hpgp.nsamps)
svs.gr$hpgp = 0
svs.gr$hpgp[ol.df$queryHits] = ol.df$freq

message('Overlap with coding sequences...')
svs.gr$cds.pc = overlapsAny(svs.gr, subset(genc, type=='CDS' & gene_type=='protein_coding'))

message('Write output TSV...')
## tier 1
## either rare coding
t1.bool = svs.gr$cds.pc & svs.gr$freq<=.01 & !svs.gr$giab & svs.gr$hpgp<=.01
## or large CNVs AF<10%
t1.bool = t1.bool | (svs.gr$cds.pc & svs.gr$freq<=.1 & svs.gr$hpgp<=.1 & svs.gr$size>1e5)
## subset and write
svs.t1 = svs.gr[which(t1.bool)]
svs.t1 %>% as.data.frame %>% write.table(file=snakemake@output[["vcf_t1"]], sep='\t', row.names=FALSE, quote=FALSE)

## other rare SVs
svs.other = svs.gr[which(!t1.bool & svs.gr$freq<=.01 & !svs.gr$giab & svs.gr$hpgp<=.01)]
svs.other %>% as.data.frame %>% write.table(file=snakemake@output[["vcf_other"]], sep='\t', row.names=FALSE, quote=FALSE)
