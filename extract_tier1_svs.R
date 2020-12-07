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
ol.gr = svOverlap(svs.gr, hpgp, min.ol=.5, max.ins.dist=100)
svs.gr$hpgp = FALSE
svs.gr$hpgp[ol.gr$queryHits] = TRUE

message('Overlap with coding sequences...')
svs.gr$cds.pc = overlapsAny(svs.gr, subset(genc, type=='CDS' & gene_type=='protein_coding'))


message('Write output TSV...')
svs.t1 = subset(svs.gr, cds.pc & freq<=.01 & !giab & !hpgp)
svs.t1 %>% as.data.frame %>% write.table(file=snakemake@output[["vcf_t1"]], sep='\t', row.names=FALSE, quote=FALSE)

svs.other = subset(svs.gr, !cds.pc & freq<=.01 & !giab & !hpgp)
svs.other %>% as.data.frame %>% write.table(file=snakemake@output[["vcf_other"]], sep='\t', row.names=FALSE, quote=FALSE)
