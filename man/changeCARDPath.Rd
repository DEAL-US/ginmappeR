\name{changeCARDPath}
\alias{changeCARDPath}
\title{Change path to download CARD database}
\description{
Defines a new path to the folder where to download CARD database. By default,
this path is the temporary folder tempdir().
}
\arguments{
  \item{path}{valid route to desired folder. Default value: tempdir()}
}
\value{
Updated path for CARD database to be downloaded.
}
\usage{
changeCARDPath(path=tempdir())
}
