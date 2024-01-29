# Add any `module load` or `export` commands that your code needs to
# compile and run to this file.
module purge
module load gcc/10.2
module load llvm/14
module load python/3.10.11-tensorflow2.10
module load maqao
module load intel/compiler/2024.0.0
module list
echo "Done."