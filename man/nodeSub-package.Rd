\name{nodeSub-package}
\alias{nodeSub-package}
\alias{nodeSub}
\docType{package}
\title{
	Package providing functions to simulate sequences under different DNA evolution models
}
\description{
    Simulate DNA sequences for the node substitution model.
    In the node substitution model, substitutions accumulate additionally
    during a speciation event, providing a potential mechanistic explanation for
    substitution rate variation. This package provides tools to simulate
    such a process, simulate a reference process with only substitutions along
    the branches, and provides tools to infer phylogenies from alignments. More
    information can be found in Janzen (2021) <doi:10.1093/sysbio/syab085>.

    Version History: \cr
    Version 1.2.7   -  Removed beta calculation due to apTreeshape removal from
                       CRAN \cr
    Version 1.2.3   -  Removed summary statistic tests for CRAN \cr
    Version 1.2.2   -  Changed codedov links in README \cr
    Version 1.2.1   -  Expanded depenency on RPANDA \cr
    Version 1.2     -  Release on CRAN \cr
}

\author{
Thijs Janzen
Maintainer: Thijs Janzen <t.janzen@rug.nl>
}

\references{
Thijs Janzen, Folmer Bokma, Rampal S Etienne, Nucleotide Substitutions during Speciation may Explain Substitution Rate Variation, Systematic Biology, 2021; syab085
}
