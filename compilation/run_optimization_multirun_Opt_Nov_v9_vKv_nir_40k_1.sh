#!/bin/csh -f
#$ -cwd
#  input           = /dev/null
#  output          = /u/scratch/d/danieleb
#$ -o /u/scratch/d/danieleb/NitrOMZ/iNitrOMZ_v6.1/compilation/Opt_Nov_v9_vKv_nir_40k/optimization_run1.joblog
#$ -j y
#$ -l h_data=2.5G,h_rt=48:00:00,arch=intel-E5-2650v4,highp
#$ -pe shared 12
#$ -M dbianchi@atmos.ucla.edu
#$ -m bea 
source /u/local/Modules/default/init/modules.csh
module load matlab
/u/scratch/d/danieleb/NitrOMZ/iNitrOMZ_v6.1/compilation/Opt_Nov_v9_vKv_nir_40k/optimize_cmaes/optimize_cmaes
