#!/bin/sh
#SBATCH --job-name=gcapc

module load R/4.1.0
#R -e "rmarkdown::render('/accounts/campus/koenvdberge/gcapc/211119_gcapcPeakCallingCalderon_scf.Rmd')"
R CMD BATCH --no-save 211119_gcapcPeakCallingCalderon_scf.R gcapc.Rout


#  cluster_kvdb.sh
#
#
#  Created by Koen Van den Berge on 10/9/20.
