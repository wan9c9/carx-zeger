#!/bin/sh
#R CMD BATCH --slave ./summarize.R
if [ -x pandoc >/dev/null ]; then 
	R -e "library(rmarkdown);render('summarize.Rmd')"
else
  R CMD BATCH --slave ./summarize.R summary.txt
fi

