\name{read_and_normalize450K}
\alias{read450K}
\alias{normalize450K}
\alias{dont_normalize450K}

\title{
Normalization of 450K data by LOESS method
}
\description{
    Read 450K '.idat' files and compute raw or normalized beta-values.
}
\usage{
read450K(idat_files)
normalize450K(intensities,tissue='')
dont_normalize450K(intensities)
}

\arguments{
    \item{idat_files}{
    a character vector containing the paths to the .idat files stripped from the '_Grn.idat' suffix with one entry for each sample (.idat files for green and red intensities have to be in the same folder).
    }
    \item{intensities}{
        List object containing raw signal intensities. Result of calling \code{read450K}.
    }
    \item{tissue}{
        If set to 'Blood', a set of prespecified reference values are used for normalization. This is recommended if you plan to use \code{estimateLC}.
    }
}
\details{
    Function \code{read450K} reads .idat files and returns a list object containing raw signal intensities. \code{dont_normalize450K} returns an ExpressionSet containing beta-values without normalization. \code{normalize450K} performs dye bias correction using the extension controls probes followed by normalization by local regression (Heiss and Brenner, 2015) and returns an ExpressionSet containing beta-values, too.
}
\value{
    For \code{read450K} a list containing the methylated, unmethylated and control signal intensities. For \code{dont_normalize450K} and \code{normalize450K} an ExpressionSet containing beta-values, rows corresponding to CpG sites (named) and columns to samples (in the same order as 'idat_files').
}

\examples{
\dontrun{
 library(minfiData) ## this package includes some .idat files
 library(data.table)

 path <- system.file("extdata",package="minfiData")
 samples = fread(file.path(path, 'SampleSheet.csv'),integer64='character')

 samples[,file:=file.path(path,Sentrix_ID,paste0(Sentrix_ID,'_',Sentrix_Position))]
 ## samples$file is a character vector containing the location of the
 ## .idat files, but without the suffixes "_Red.idat" or "_Grn.idat"

 raw = read450K(samples$file)
 none = dont_normalize450K(raw) ## no normalization
 norm = normalize450K(raw)
}
}
\author{
    Jonathan A. Heiss
}

\note{
    A benchmark comparing the performance of this method with other normalization approaches is provided in the vignette.
}

\references{
    Heiss, Jonathan A., and Hermann Brenner. \emph{Between-array normalization for 450K data}. Frontiers in genetics 6 (2015). doi:\href{http://dx.doi.org/10.3389/fgene.2015.00092}{10.3389/fgene.2015.00092}
}
