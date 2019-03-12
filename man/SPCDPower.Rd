\name{SPCDPower}
\alias{SPCDPower}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the power, sample size for a Sequential Parallel Comparison Design(SPCD) study
}
\description{
Patients are randomized in three groups, patients who receive placebo in phase 1 and again in phase 2 of the study, patients who receive placebo in phase 1 and active in phase 2 and patients who receive active therapy in phase 1 and are not included in phase 2. A response criteria is determined and the phase 2 data of patients who respond in phase 1 is eliminated. Each phase is analyzed separately and the results are pooled. Calculates power or sample size as a function of the alternative hypothesis, posed in terms of response rates or effect sizes, for both binary and continuous outcomes.
}

\usage{
SPCDPower(n=NULL, power=NULL, p, w=0.5, placeboProp=.66, drop = 0, alpha = 0.025,
         effect_size = rep(NULL, 2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{n}{Total sample size of the study, leave as null if you want the sample size computed.
}
  \item{power}{Power of the study, leave as null if you want the power computed.
}
  \item{p}{ A 2 by 2 matrix, \code{matrix(c(Phase1.response.drug,phase1.response.placebo,phase2.response.drug,phase2.response.placebo),2,2)} indicating the alternative hypothesis
}
  \item{w}{Weight for the first phase in the combined test
%%     ~~Describe \code{w} here~~
}
  \item{placeboProp}{
Proportion of patients randomized to placebo in the first phase
}

  \item{drop}{ The proportion of placebo non-responders that drop after the first phase
%%     ~~Describe \code{drop} here~~
}
  \item{alpha}{
Significance level
}
  \item{effect_size}{
This is an alternative method of specifying the alternative. If it is used only p[2,1] needs to be specified. This is useful in the situation where a continuous endpoint is used and treatment response is not defined as the endpoint being greater than a constant.
}
}
\details{
This program considers the situation in which response rates are supplied by the investigator, response is judged as by whether or not a continuous variable is greater than a constant, and the continuous variable is analyzed rather than the response variable. In this case it turns out the effect size for a comparison, where the response rates are \code{p} and \code{q} for placebo and active drug is \code{qnorm(1-p)-qnorm(1-q)}.

}
\value{
A numeric vector with the following fields, sample size \code{n}, Power for for the SPCD when using a dichotomous response outcome, Power for the SPCD using a continuous outcome where response is judged as a continuous variable being greater than a fixed constant, Power for a conventional design for a dichotomous variable and a continuous variable, Power for a SPCD design where the null is rejected if either the first phase or the second phase shows a significant difference. The first value is not corrected for multiple comparisons while the second uses a bonferroni correction.

}
\references{
Fava, M., Evins, A. E., Dorer, D. J., and Schoenfeld, D. A. (2003). The problem of the placebo response in clinical trials for psychiatric disorders: culprits, possible remedies, and a novel study design approach. \emph{Psychotherapy and psychosomatics}, \bold{72},3, 115--127.

Tamura, R. N., & Huang, X. (2007). An examination of the efficiency of the sequential parallel design in psychiatric clinical trials. \emph{Clinical Trials}, \bold{4},4, 309-31.

}
\author{
David A. Schoenfeld dschoenfeld@mgh.harvard.edu
}

\seealso{
\code{\link{SPCDbinary}},\code{\link{SPCDcontinuous}}
}
\examples{
SPCDPower(n=150, power=NULL, p=matrix(c(.6,.3,.5,.3),2,2), w=0.5,
placeboProp=.66, drop = .1, alpha = 0.025,effect_size = rep(NULL, 2))
}
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")

