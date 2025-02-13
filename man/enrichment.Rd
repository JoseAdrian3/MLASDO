% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment.R
\name{enrichment}
\alias{enrichment}
\title{Enrichment analysis}
\usage{
enrichment(
  importance,
  saving_name,
  text,
  enrichment,
  omic_data_changed,
  omic_data_maintained,
  bool
)
}
\arguments{
\item{importance}{Data frame containing variables and their importance scores in $importance}

\item{saving_name}{String for the saving directory and file prefix.}

\item{text}{String to append to file names.}

\item{enrichment}{String with the enrichment type: "transcriptomic" or "metabolomic".}

\item{omic_data_changed}{Data frame of omics data with the changed condition.}

\item{omic_data_maintained}{Data frame of omics data with the original condition.}

\item{bool}{String indicating condition: "changed" or other.}
}
\value{
Saves enrichment analysis results and elbow data.
}
\description{
This function performs enrichment analysis based on the most important genes
in the dataset using the elbow rule method. Omic may be transcriptomic,
proteomic or metabolomic.
}
