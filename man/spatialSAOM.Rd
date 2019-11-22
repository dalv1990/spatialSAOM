\name{spatialSAOM}

\title{Spatial Stochastic Actor-Oriented Model}

\description{Fit a SAOM regression model based on a static, structural network and a set of covariates. The typical application is a diffusion model, where the network is defined by geographic adjacency, and the dependent variable is binary or categorical.}

\usage{
spatialSAOM(formula, data, subset, network, diffusion = list(), rateFix = 20, maxRounds = 10, method = "avSim", projname = "SAOM", ...)
}

\arguments{
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. Contains dependent variable and covariates, excluding the diffusion effect.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which spatialSAOM is called.}
  \item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}
  \item{network}{square connection matrix of all observations. The dimension must match the number of observations in the data frame.}
  \item{diffusion}{a list of covariates through which diffusion takes place, if not through the dependent variable.}
  \item{rateFix}{rate at which behavioral changes are fixed to guarantee estimation.}
  \item{maxRounds}{maximum number of iterations of the SAOM algorithm.}
  \item{method}{the method to be used; for fitting, currently only "avAlt" and "avSim" are supported.}
  \item{projname}{name under which temporary output is saved by the SAOM implementation.}
  \item{\dots}{not used.}
}

\value{
  Returns an object of class \code{\link{sienaFit}}.
}

\author{Johan A. Elkink and Thomas U. Grund}

\references{\itemize{
  \item Elkink, Johan A. and Thomas U. Grund. 2019. "Modelling Diffusion through Statistical Network Analysis: A Simulation Study." Arxiv 1903.08648.
  \item Snijders, Tom A.B. and Christian E.G. Steglich. 2015. "Representing Micro–Macro Linkages by Actor-based Dynamic Network Models." Sociological Methods & Research 44(2):222– 271.
  }
}


\keyword{}