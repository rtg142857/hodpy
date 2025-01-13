#!/bin/bash -l

# First argument is the path to the path_config.yml file

#SBATCH --ntasks 1 # The number of cores you need...
#SBATCH -J get_mass_function #Give it something meaningful.
#SBATCH -o logs/out_mass_function
#SBATCH -e logs/err_mass_function
#SBATCH -p cosma8 #or some other partition, e.g. cosma, cosma8, etc.
#SBATCH -A dp004
#SBATCH --exclusive
#SBATCH -t 60
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=tlrt88@durham.ac.uk #PLEASE PUT YOUR EMAIL ADDRESS HERE (without the <>)

module purge
#cosmodesienv main
#conda activate abacus-env
module use /cosma/apps/dp004/dc-mene1/desi/cosmodesiconda/my-desiconda/modulefiles
module load cosmodesiconda/my-desiconda
#source /cosma/home/dp004/dc-mene1/activate-cosmodesiconda

#module load gcc
#module load gsl
#module unload craype-hugepages2M

#module load python/3.10.12

python get_mass_function_single_redshift_flamingo.py $1
