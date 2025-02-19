% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_results.R
\name{annotate_results}
\alias{annotate_results}
\title{Annotation of DE results}
\usage{
annotate_results(input, data, ensembl2id, pathway_dic, circ = F)
}
\arguments{
\item{input}{UI input, includes species}

\item{data}{a DE result data.frame including either ensembl_gene_id as a column (for circRNA data) or row names for linRNA.}

\item{ensembl2id}{a reactive data.frame including gene symbol information}

\item{pathway_dic}{a reactive data.frame with the pathways involved for each gene}

\item{circ}{boolean indicating if the data is of circRNA data structure or linear.}
}
\value{
data.frame
}
\description{
Translates Ensembl gene IDs to gene symbol, information on Transcription factor status and the pathways in which the gene is involved. It adds secondary DE analysis results along with annotation from chosen files.
}
