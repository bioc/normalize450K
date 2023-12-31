\name{estimateLC}
\alias{estimateLC}

\title{
Estimation of Leukocyte composition for whole blood samples
}
\description{
    Estimate leukocyte composition from whole blood DNA methylation
}
\usage{
estimateLC(eSet)
}

\arguments{
    \item{eSet}{
        A Biobase eSet object as returned from a call of \code{normalize450K}
    }
}
\details{
    Cell proportions are estimated using the algorithm developed by Houseman et al. (2012) by two different models. The first model was trained on a dataset of purified leukocytes (Reinius et al., 2012) and provides predictions for six cell types (granulocytes, monocytes, CD8+ T cells, CD4+ T cells, natural killer cells and CD19+ B cells), the second model was trained on whole blood samples from the LOLIPOP study as described by Heiss et al. (2016) and provides predictions for 4 cell types (neutrophils, eosinophils, lymphocytes, monocytes -- ignore the prediction for basophils). Use this function only for normalized data (with \code{normalize450K(...,tissue='Blood')}.
}
\value{
    Returns the eSet object with cell proportions estimates added to the phenoData slot.
}

\author{
    Jonathan A. Heiss
}

\references{
    Houseman EA, et al. (2012) DNA methylation arrays as surrogate measures of cell mixture distribution. BMC Bioinformatics, doi:\href{http://dx.doi.org/10.1186/1471-2105-13-86}{10.1186/1471-2105-13-86}

    Reinius LE, et al. (2012) Differential DNA methylation in purified human blood cells: implications for cell lineage and studies on disease susceptibility. PloS ONE, doi:\href{http://dx.doi.org/10.1371/journal.pone.0041361}{10.1371/journal.pone.0041361}

    Heiss JA, et al. (2016). Training a model for estimating leukocyte composition using whole-blood DNA methylation and cell counts as reference. Epigenomics, doi:\href{http://dx.doi.org/10.2217/epi-2016-0091}{10.2217/epi-2016-0091}
}


