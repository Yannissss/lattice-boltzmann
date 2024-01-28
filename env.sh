# Add any `module load` or `export` commands that your code needs to
# compile and run to this file.
module purge
module load gcc/10.2
module load llvm/14
module list
echo "Done."