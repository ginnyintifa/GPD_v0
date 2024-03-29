% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/somatic_extraction_mapping_wrap.R
\name{somatic_extraction_annotation_pos}
\alias{somatic_extraction_annotation_pos}
\title{Extract germline mutations annotated by the snpEff tool for a cancer type of your interest.}
\usage{
somatic_extraction_annotation_pos(mc3_info_file, cancer_type,
  cancer_barcode, output_dir)
}
\arguments{
\item{mc3_info_file}{File name of the mc3 somatic variants from TCGA.}

\item{cancer_type}{TCGA cancer type identifier, e.g. STAD, BRCA.}

\item{output_dir}{The directory you would like to have your output files in.}

\item{sample_cn_id}{A barcode codebook provided by the packge.}
}
\description{
This function generates two files as the final products. First, mutations which result in protein-consequence with
effect on proteins annotated; Second, mutations without protein-consequence.
}
\examples{
somatic_extraction_annotation_pos(mc3_info_file = "/data/ginny/tcga_pancan/important_files/mc3_info.tsv" ,
                                 cancer_type = "STAD",
                                 cancer_barcode = stad_barcode,
                                 output_dir = "/data/ginny/tcga_pancan/STAD_somatic/")
}
