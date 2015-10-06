#!/bin/sh
R CMD BATCH --no-save --no-restore "--args iRep=${1}" ./simulation_single.R
