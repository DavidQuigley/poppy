% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poppy.R
\name{Genome}
\alias{Genome}
\title{Constuctor for object that holds metadata about the genome itself}
\usage{
Genome(
  fn_gene_locs,
  fn_cyto,
  fn_chrom_lengths,
  fn_centromeres,
  fn_ensembl,
  fn_exons
)
}
\arguments{
\item{fn_gene_locs}{TSV, chrom start end NM (refseq_id) score (ignored) strand}

\item{fn_cyto}{TSV, Chromosome ChromStart ChromEnd Band Stain}

\item{fn_chrom_lengths}{TSV, (no header) chrom length}

\item{fn_centromeres}{TSV, chrom chromStart (pos) chromEnd (pos) chrom_idx}

\item{fn_ensembl}{TSV, generated from GTF, ensembl_id symbol chrom start end strand gene_type}

\item{fn_exons}{TSV, symbol (no individual exon identifier) chrom start end strand}
}
\value{
A new Genome object
}
\description{
Constuctor for object that holds metadata about the genome itself
}
