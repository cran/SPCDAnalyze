\name{runSlopes}
\alias{runSlopes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Runs a random slopes model for use with \code{SPCDcontinuous}
}
\description{
The function \code{SPCDFcontinuous} requires an analysis function which estimates the treatment effect and it's standard error for each phase of the study. This function uses the random slopes model to analyze the data.This is a random effects model with a random slope and intercept for each patient. The treatment effect is the difference in the mean slopes which modeled as the time-treatment interaction.  This function can either be used as is or be used as a model for a user supplied function.
}
\usage{
runSlopes(combData, times = "time")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{combData}{
This is the input file data file, it needs to have fields \code{ID} which is the patient identifier,
\code{time} which is the observation time, \code{y} the outcome variable, \code{trt} a treatment code which is either
\code{Placebo} or \code{Active}
}
  \item{times}{
This is the time variable.The observation time can be different for each patient.  If they are the same then the \code{runGls} function might also be appropriate.
}
}
\details{
This fits the random slopes model
\eqn{y_{t}=\mu+\beta_1 t+ \beta_2 t*I(rx=1)+u+b t+\sigma \epsilon},
where \eqn{u,b,\epsilon} are random variables. Note that a treatment main effect is not included in the model by default, because in a randomized study the treatments should be the same at the baseline visit. This practice may vary.
}
\value{
A numerical vector of length 3 which gives the treatment effect, it's standard are and the z-score which is the ratio of the treatment effect to it's standard error. The function \code{SPCDcontinuous} will use the number of patients as the degrees of freedom to compute p-values.
}
\references{
Laird, N. M., & Ware, J. H. (1982). Random-effects models for longitudinal data. \emph{Biometrics}, 963-974.}

\author{
David A. Schoenfeld dschoenfeld@mgh.harvard.edu
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{runGls}},\code{\link{SPCDcontinuous}}
}
\examples{
data1=createTestData(times=0:11,n=c(50,50,50),transition=5,rx1=0,rx2=0)
data1=placeboNonResponder(data1,k=0)
#Example of use of function outside of SPCDcontinous
data1$trt=ifelse(data1$rx==3,'Active','Placebo')
runSlopes(data1,time='time')

}

\keyword{htest}% use one of  RShowDoc("KEYWORDS")

