# lfric_hj_bench_code
Scripts to plot results of LFRic runs for hot Jupiters

# Install ANTS and its dependencies

```bash
# Create a conda environment
mamba create -n lfric_ana python=3.12 ipykernel pip aeolus esmpy=8.4 gdal iris matplotlib numba pykdtree
# Download mule and install it
svn co https://code.metoffice.gov.uk/svn/um/mule/trunk@119685 mule_2023.08.1
cd mule_2023.08.1/mule
python setup.py install
# Download and install ANTS
svn co https://code.metoffice.gov.uk/svn/ancil/ants/tags/2.0.0 ants_2.0.0
cd ants_2.0.0
python setup.py install
# Test the installation
python -c 'import ants; print(ants.__file__)'
# Download ants_contrib
svn co https://code.metoffice.gov.uk/svn/ancil/contrib/tags/2.0.0 contrib_2.0.0
```
