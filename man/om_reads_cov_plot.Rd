% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/om_reseq_basic.R
\name{om_reads_cov_plot}
\alias{om_reads_cov_plot}
\title{Plot Reads coverage distribution}
\usage{
om_reads_cov_plot(coverage_table, sample_limit, max_depth, out_prefix = NULL)
}
\arguments{
\item{coverage_table}{coverage stats table of all samples}

\item{sample_limit}{sample number limit to show detail inform of each sample}

\item{max_depth}{Sequencing depth limit to show in plot}

\item{out_prefix}{output file prefix, default is NULL, don't output file}
}
\description{
Plot bwa mapping coverage summrized by samtools
}
\details{
Mapping coverage of each sample is generated by \code{samtools stats sorted.bam};

Then extracted using \code{grep}

\code{grep ^COV| cut -f 3,4}
}
\examples{

cov_stats <- system.file("extdata", "all_sample.genome.cov.xls", package = "omplotr")

# show coverage summary
om_reads_cov_plot(cov_stats, 5, 100)
# show detail sample information
om_reads_cov_plot(cov_stats, 10, 100)
}
