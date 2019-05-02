\name{ipwErrorY-package}
\alias{ipwErrorY-package}
\docType{package}
\title{
Inverse Probability Weighted (IPW) Estimation of Average Treatment Effect (ATE) with Misclassified Binary Outcome
}
\description{
This package is an implementation of the correction methods proposed by Shu and Yi (2017) for the inverse probability weighted (IPW) estimation of average treatment effect (ATE) with misclassified binary outcomes. Logistic regression model is assumed for treatment model for all implemented correction methods, and is assumed for the outcome model for the implemented doubly robust correction method. Misclassification probability given a true value of the outcome is assumed to be the same for all individuals.
}
\details{
The \code{ipwErrorY} package implements correction methods developed by Shu and Yi (2017) to adjust for misclassification in binary outcomes in the inverse probability weighted estimation of average treatment effect. The function \code{\link[ipwErrorY]{KnownError}} implements the correction method with known outcome misclassification probabilities. The function \code{\link[ipwErrorY]{EstValidation}} implements the optimal linear combination correction method when validation data are available. The function \code{\link[ipwErrorY]{Est2Replicates}} implements the correction method when two independent replicates of the outcome are available. The function \code{\link[ipwErrorY]{KnownErrorDR}} implements the doubly robust correction method with known outcome misclassification probabilities.
}
\author{
Di Shu and Grace Y. Yi

Maintainer: Di Shu <shudi1991@gmail.com>
}
\references{
Shu D, Yi GY. (2017). Causal inference with measurement error in outcomes: bias analysis and estimation methods. \emph{Statistical Methods in Medical Research} <doi:10.1177/0962280217743777>
}
\keyword{ package }

