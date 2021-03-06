% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dociMoMo.R, R/iMoMo.R
\docType{package}
\name{iMoMo}
\alias{iMoMo}
\title{iMoMo.}
\usage{
iMoMo(
  staticAgeFun = TRUE,
  periodAgeFun = "NP",
  cohortAgeFun = NULL,
  type = c("indirect", "direct"),
  constFun = function(ax, bx, kt, b0x, gc, wxt, ages) list(ax = ax, bx = bx, kt = kt,
    b0x = b0x, gc = gc),
  constFunEst = function(Ax, ax, bx, kt, b0x, gc, wxt, ages) list(Ax = Ax, ax = ax, bx
    = bx, kt = kt, b0x = b0x, gc = gc)
)
}
\arguments{
\item{staticAgeFun}{logical value indicating if a static age function
\eqn{\alpha_x} is to be included.}

\item{periodAgeFun}{a list of length \eqn{N} with the definitions of the
period age modulating parameters \eqn{\beta_x^{(i)}}. Each entry can take
values: \code{"NP"} for non-parametric age terms, \code{"1"} for
\eqn{\beta_x^{(i)}=1} or a predefined parametric function of
age (see details). Set this to \code{NULL} if there are no period terms
in the model.}

\item{cohortAgeFun}{defines the cohort age modulating parameter
\eqn{\beta_x^{(0)}}. It can take values: \code{"NP"} for non-parametric
age terms, \code{"1"} for \eqn{\beta_x^{(0)}=1}, a predefined parametric
function of age (see details) or \code{NULL} if there is no cohort effect.}

\item{type}{defines the type of estimation method to be used.
\code{"indirect"} would estimate the equivalent mortality rate model
and the transform the model into improvement rates. \code{"direct"}
would estimate the model directly on the improvement rate model.}

\item{constFun}{function defining the identifiability constraints of the
model. It must be a function of the form
\code{constFun <- function(ax, bx, kt, b0x, gc, wxt, ages)} taking a set
of fitted model parameters and returning a list
\code{list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)}
of the model parameters with the identifiability constraints applied. If
omitted no identifiability constraints are applied to the model.}

\item{constFunEst}{function defining the identifiability constraints for the
equivalent mortality rate model. It must be a function of the form
\code{constFunEst <- function(Ax, ax, bx, kt, b0x, gc, wxt, ages)} taking a set
of fitted model parameters and returning a list
\code{list(Ax = Ax, ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)}
of the model parameters with the identifiability constraints applied. If
omitted no identifiability constraints are applied to the estimaion model.}
}
\description{
Initialises an iMoMo oject which represents a generalised
Age-Period-Cohort improvement rate mortality model
}
