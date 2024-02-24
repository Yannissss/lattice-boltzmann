# Add any `module load` or `export` commands that your code needs to
# compile and run to this file.
echo "Load required modules..."
module purge
module load gcc/10.2
module load llvm/14
module load python/3.10.11-tensorflow2.10
module load maqao
module load intel/compiler
module load openmpi/4.1.4.2-intel-ilp64
echo "Done."