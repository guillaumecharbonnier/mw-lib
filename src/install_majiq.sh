# test Majiq

rm -rf /media/gcharbonnier/samsung_m3/miniconda/envs/majiq

conda env create -f /media/gcharbonnier/samsung_m3/mw-lib/src/snakemake/envs/majiq.yaml

conda activate majiq

export HTSLIB_LIBRARY_DIR=/media/gcharbonnier/samsung_m3/miniconda/envs/majiq/libexec/htslib/
export HTSLIB_INCLUDE_DIR=/media/gcharbonnier/samsung_m3/miniconda/envs/majiq/include/htslib


pip install git+https://bitbucket.org/biociphers/majiq_academic.git#egg=majiq


/media/gcharbonnier/samsung_m3/miniconda/envs/majiq/lib/python3.8/site-packages/majiq/run_majiq.py










conda deactivate
python3 -m venv env
source env/bin/activate

export HTSLIB_LIBRARY_DIR=/media/gcharbonnier/samsung_m3/miniconda/envs/majiq/libexec/htslib/
export HTSLIB_INCLUDE_DIR=/media/gcharbonnier/samsung_m3/miniconda/envs/majiq/include/htslib/



export HTSLIB_LIBRARY_DIR=/usr/local/lib

pip install git+https://bitbucket.org/biociphers/majiq_academic.git#egg=majiq







##### AFTER THIS CHUNK, MAJIQ WAS WORKING:

export HTSLIB_LIBRARY_DIR=/usr/local/lib
export HTSLIB_INCLUDE_DIR=/usr/local/include

pip install majiq_academic/

