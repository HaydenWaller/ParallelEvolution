library(tidyverse)


#Merge count matrices
lk.matrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/matrices/lk.isoform.counts.matrix")
ln.matrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/matrices/ln.isoform.counts.matrix")
lp.matrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/matrices/lp.isoform.counts.matrix")
ls.matrix<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/matrices/ls.isoform.counts.matrix")

lk.lp.matrix<-merge(lk.matrix, lp.matrix, by = "Orthogroup")
lk.lp.ln.matrix<-merge(lk.lp.matrix, ln.matrix, by = "Orthogroup")
lk.lp.ln.ls.matrix<-merge(lk.lp.ln.matrix, ls.matrix, by = "Orthogroup")

write.table(lk.lp.ln.ls.matrix, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/matrices/merged.isoform.counts.matrix")


#Merge TMM EXPR matrices
lk.tmm.expr.matrix<-read.table("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/matrices/lk.isoform.TMM.EXPR.matrix", header = TRUE)
ln.tmm.expr.matrix<-read.table("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/matrices/ln.isoform.TMM.EXPR.matrix", header = TRUE)
lp.tmm.expr.matrix<-read.table("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/matrices/lp.isoform.TMM.EXPR.matrix", header = TRUE)
ls.tmm.expr.matrix<-read.table("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/matrices/ls.isoform.TMM.EXPR.matrix", header = TRUE)

lk.lp.tmm.expr.matrix<-merge(lk.tmm.expr.matrix, lp.tmm.expr.matrix, by = "Orthogroup")
lk.lp.ln.tmm.expr.matrix<-merge(lk.lp.tmm.expr.matrix, ln.tmm.expr.matrix, by = "Orthogroup")
lk.lp.ln.ls.tmm.expr.matrix<-merge(lk.lp.ln.tmm.expr.matrix, ls.tmm.expr.matrix, by = "Orthogroup")

write.table(lk.lp.ln.ls.tmm.expr.matrix, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/matrices/merged.TMM.EXPR.matrix")


#merge pairwise results files with blastp

blast<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/lk_scOrthoProts.pep_blastp.outfmt6")
orthos<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/lk_transcripts_and_orthos.txt")
annot_ortho<-merge(blast, orthos, by="Transcript", all = TRUE)

write_tsv(annot_ortho, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/annots")


annots<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/annots_R.txt")
lk_x_ln<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/DESeq2/merged.isoform.counts.ready.matrix.KOHALENSIS_vs_NEOSPISA.DESeq2.DE_results")
lk_x_lp<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/DESeq2/merged.isoform.counts.ready.matrix.KOHALENSIS_vs_PARANIGRA.DESeq2.DE_results")
lk_x_ls<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/DESeq2/merged.isoform.counts.ready.matrix.KOHALENSIS_vs_SPISA.DESeq2.DE_results")
ln_x_lp<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/DESeq2/merged.isoform.counts.ready.matrix.NEOSPISA_vs_PARANIGRA.DESeq2.DE_results")
ln_x_ls<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/DESeq2/merged.isoform.counts.ready.matrix.NEOSPISA_vs_SPISA.DESeq2.DE_results")
lp_x_ls<-read_tsv("C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/DESeq2/merged.isoform.counts.ready.matrix.PARANIGRA_vs_SPISA.DESeq2.DE_results")

lk_vs_ln.annotated.results<-merge(lk_x_ln, annots, by="Orthogroup", all=TRUE)
lk_vs_lp.annotated.results<-merge(lk_x_lp, annots, by="Orthogroup", all=TRUE)
lk_vs_ls.annotated.results<-merge(lk_x_ls, annots, by="Orthogroup", all=TRUE)
ln_vs_lp.annotated.results<-merge(ln_x_lp, annots, by="Orthogroup", all=TRUE)
ln_vs_ls.annotated.results<-merge(ln_x_ls, annots, by="Orthogroup", all=TRUE)
lp_vs_ls.annotated.results<-merge(lp_x_ls, annots, by="Orthogroup", all=TRUE)


write.table(lk_vs_ln.annotated.results, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/DESeq2/lk_vs_ln.annotated.results")
write.table(lk_vs_lp.annotated.results, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/DESeq2/lk_vs_lp.annotated.results")
write.table(lk_vs_ls.annotated.results, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/DESeq2/lk_vs_ls.annotated.results")
write.table(ln_vs_lp.annotated.results, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/DESeq2/ln_vs_lp.annotated.results")
write.table(ln_vs_ls.annotated.results, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/DESeq2/ln_vs_ls.annotated.results")
write.table(lp_vs_ls.annotated.results, file="C:/Users/hayde/Documents/GRAD_SKEWL/Projects_Experiments/5_species_RNAseq/no_tantalis/DESeq2/lp_vs_ls.annotated.results")

