# Initialise conda
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate gadcp
hg clone https://currents.soest.hawaii.edu/hgstage/geomag
cd geomag
# Clean the version control stuff
rm -rf .hg .hgignore
# Compile
make magdec