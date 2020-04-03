.onLoad <- function(libname, pkgname) {
.jpackage(pkgname, lib.loc = libname)

defaults <- list(
  bioconductorSELEX.workdir = NULL,
  bioconductorSELEX.verbose = FALSE,
  bioconductorSELEX.maxthreads = NULL
)
toSet <- !(names(defaults) %in% names(options()))
if(any(toSet)) options(defaults[toSet])

workdir = getOption("bioconductorSELEX.workdir")
if(is.null(workdir)) workdir=file.path(tempdir(), "bioconductorSELEX")
verbose = getOption("bioconductorSELEX.verbose")
maxthreads = getOption("bioconductorSELEX.maxthreads")

selex.config(workdir, verbose, maxthreads)
}
