% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNA.R
\name{somatic_add_CNA_copycat}
\alias{somatic_add_CNA_copycat}
\title{Wrapper function to add CNA segments, filter, and adjust for ploidy/purity using CLONET}
\usage{
somatic_add_CNA_copycat(
  genome,
  somatic,
  fn_segs,
  fn_bins,
  fn_pileup_normal,
  fn_pileup_tumor,
  min_segment_percent = 0.005,
  max_percent_centromere_overlap = 0.2,
  purity_override = NA,
  ploidy_override = NA,
  min_required_SNPs = 10,
  sex = "M",
  n_cores = 1,
  verbose = FALSE
)
}
\arguments{
\item{genome}{object created by call to Genome()}

\item{somatic}{somatic data object to modify}

\item{fn_segs}{path to segment file}

\item{fn_bins}{path to binned copy numbers}

\item{fn_pileup_normal}{pileup file created by ASEQ program}

\item{fn_pileup_tumor}{pileup file created by ASEQ program}

\item{min_segment_percent}{minimal size of a segment expressed as a percent of the chromosome size}

\item{max_percent_centromere_overlap}{eliminate segments where total overlap exceeds this percent, 0-1}

\item{purity_override}{use this value instead of estimating purity}

\item{ploidy_override}{use this value instead of estimating ploidy}

\item{min_required_SNPs}{minimum number of SNPs in a segment to make a call}

\item{sex}{M or F}

\item{n_cores}{allow for multithreaded calculations}

\item{verbose}{TRUE or FALSE}
}
\description{
Wrapper function to add CNA segments, filter, and adjust for ploidy/purity using CLONET
}
