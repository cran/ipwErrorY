\name{ipwErrorY-package}
\alias{ipwErrorY-package}
\docType{package}
\title{
Inverse Probability Weighting (IPW) Estimation of Average Treatment Effect (ATE) with Outcome Misclassification
}
\description{
This package is an implementation of the correction methods proposed by Shu and Yi (2017) for the inverse probability weighting (IPW) estimation of average treatment effect (ATE) with misclassified outcomes. Logistic regression model is assumed for treatment model for all implemented correction methods, and is assumed for the outcome model for the implemented doubly robust correction method.
}
\details{
The \code{ipwErrorY} package implements correction methods developed by Shu and Yi (2017) to adjust for outcome misclassification in the estimation of average treatment effect. The function \code{\link[ipwErrorY]{KnownError}} implements the correction method with known outcome misclassification probabilities. The function \code{\link[ipwErrorY]{EstValidation}} implements the optimal linear combination correction method when validation data are available. The function \code{\link[ipwErrorY]{Est2Replicates}} implements the correction method when two independent replicates of the outcome are available. The function \code{\link[ipwErrorY]{KnownErrorDR}} implements the doubly robust correction method with known outcome misclassification probabilities.
}
\author{
Di Shu and Grace Y. Yi

Maintainer: Di Shu <shudi1991@gmail.com>
}
\references{
Shu D, Yi GY. (2017). Causal inference with measurement error in outcomes: bias analysis and estimation methods. \emph{Statistical Methods in Medical Research} <doi:10.1177/0962280217743777>
}
\keyword{ package }

