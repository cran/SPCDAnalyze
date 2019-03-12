\name{SPCDbinary}
\alias{SPCDbinary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Analyze a Sequential Parallel Comparison Design(SPCD), study with a binary endpoint
%%  ~~function to do ... ~~
}
\description{
 Patients are randomized in three groups, patients who receive placebo in phase 1 and again in phase 2 of the study, patients who receive placebo in phase 1 and active in phase 2 and patients who receive active therapy in phase 1 and are not included in phase 2. A response criteria is determined and the phase 2 data of patients who respond in phase 1 is eliminated. Each phase is analyzed separately and the results are pooled. The idea of the design is that by removing the placebo responders in the analysis in the second phase the treatment effect will be greater and their will be an increase in power, both from the reuse of patients and from enrichment in the second phase
}
\usage{
SPCDbinary(results, w = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{results}{
This is a 3 by 4 matrix where with the first column being the number of responders in phase 1, in the 3 groups, Placebo-Placebo, Placebo-Drug and Drug-Drug. The second column is the number of non-responders, the third, the number of responders in the second phase and the forth the number of non-responders in the second phase.
}
  \item{w}{
The weight using to combine the first and second phases.
}
}

\value{
A numeric vector with the effect, the standard error and the ratio of the effect to the standard error, for each phase and the pooled analysis followed by the three significance levels.
}
\references{
Fava, M., Evins, A. E., Dorer, D. J., and Schoenfeld, D. A. (2003). The problem of the placebo response in clinical trials for psychiatric disorders: culprits, possible remedies, and a novel study design approach. \emph{Psychotherapy and psychosomatics}, \bold{72},3, 115--127.

Tamura, R. N., & Huang, X. (2007). An examination of the efficiency of the sequential parallel design in psychiatric clinical trials. \emph{Clinical Trials}, \bold{4},4, 309-31.
}
\author{
David A. Schoenfeld
}
\seealso{
\code{\link{SPCDcontinuous}}
}
\examples{
SPCDbinary(matrix(c(18,15,3,11,21,12,7,4,20,14,NA,NA),3,4,byrow=TRUE))
}

\keyword{htest}

