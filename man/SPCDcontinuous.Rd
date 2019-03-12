\name{SPCDcontinuous}
\alias{SPCDcontinuous}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Analyze continuous data for the Sequential Parallel Comparison Design(SPCD) study
%%  ~~function to do ... ~~
}
\description{ Patients are randomized in three groups, patients who receive placebo in phase 1 and again in phase 2 of the study, patients who receive placebo in phase 1 and active in phase 2 and patients who receive active therapy in phase 1 and are not included in phase 2. A response criteria is determined and the phase 2 data of patients who respond in phase 1 is eliminated. Each phase is analyzed separately and the results are pooled. The idea of the design is that by removing the placebo responders in the analysis in the second phase the treatment effect will be greater and their will be an increase in power, both from the reuse of patients and from enrichment in the second phase
}
\usage{
SPCDcontinuous(combData, recordID, times, group, transition, nonResponder,outcome,
runmod = runSlopes,w=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{combData}{The data in long form that is to be analyzed.
%%     ~~Describe \code{combData} here~~
}
  \item{recordID}{The name or column number of the variable indicating the patient identifier.
%%     ~~Describe \code{recordID} here~~
}
  \item{times}{The name or column number of the variable indicating observation time.
%%     ~~Describe \code{times} here~~
}
  \item{group}{The name or column number of the variable indicating randomization group as defined above
%%     ~~Describe \code{group} here~~
}
  \item{transition}{ The time that is the end of phase 1 and the beginning of phase 2.
%%     ~~Describe \code{transition} here~~
}
  \item{nonResponder}{The name or column number of the logical variable indicating that the patient should be included in the phase 2 analysis
%%     ~~Describe \code{nonResponder} here~~
}
\item{outcome}{The name or column number of the outcome variable}
  \item{runmod}{A function to analyze the data from each phase. It must output a vector of three numbers, the treatment effect, it's standard error and the ratio of the treatment effect to its standard error. This package includes two examples \code{runSlopes} which uses a random slopes effects model, with a random slope and intercept, where the time-treatment interaction is the measure of treatment effect and \code{runGls} which is a repeated measures model where there is not treatment effect at time 0 and the outcome is the treatment effect at the last observation time.  These models are described in more detail in the documentation of the package \code{LPower} which gives sample size methods for these models. You need to make sure these functions are selecting the correct coefficient if you use them.
}
\item{w}{The weight given to the first phase results, the weight given to the second phase is
\code{1-w}.
}

}
\details{The program uses the number of patients analyzed in each phase as the degree's of freedom for calculating the p-value and uses the Welch–Satterthwaite equation to calculate the p-value of the pooled result.
}
\value{A numeric vector with the treatment effect, its standard error and the ratio for each phase as well as the pooled values.
}
\references{
Fava, M., Evins, A. E., Dorer, D. J., and Schoenfeld, D. A. (2003). The problem of the placebo response in clinical trials for psychiatric disorders: culprits, possible remedies, and a novel study design approach. \emph{Psychotherapy and psychosomatics}, \bold{72},3, 115--127.

Chen, Y. F., Yang, Y., Hung, H. J., & Wang, S. J. (2011). Evaluation of performance of some enrichment designs dealing with high placebo response in psychiatric clinical trials. \emph{Contemporary clinical trials}, \bold{32},4, 592-604.
}
\author{
David A. Schoenfeld, dschoenfeld@mgh.harvard.edu
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{createTestData}},\code{\link{runSlopes}},\code{\link{runGls}}}

\examples{
data1=createTestData(times=0:11,n=c(50,50,50),transition=5,rx1=0,rx2=0)
data1=placeboNonResponder(data1,k=0)
SPCDcontinuous(data1,'ID','time','rx','transition',
  nonResponder='nonResponder',outcome='y',w=.5,runmod=runSlopes)
data1=createTestData(times=c(0,5,11),n=c(50,50,50),
  transition=5,rx1=0,rx2=0)
data1=placeboNonResponder(data1,k=0)
SPCDcontinuous(data1,'ID','time','rx','transition',
  nonResponder='nonResponder',outcome='y',w=.5,runmod=runGls)
}
\keyword{htest}



