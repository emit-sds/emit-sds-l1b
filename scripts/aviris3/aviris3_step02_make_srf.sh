# David R Thompson

# This script relies on the caller to have first dark-subtracted and pedestal-shift-corrected all 
# the SRF data from the monochromator sweeps
# It outputs text files, one per base file, corresponding to the SRFs

wavelengthfile=/home/drt/src/emit-sds-l1b-drt/data/aviris3/AVIRIS3_Wavelengths_20230609.txt

fin=/beegfs/scratch/drt/20230608_AVIRIS3_SRF/Col_919_VNIR_darksub_pedestal 
#python ../../utils/makesrf.py --top 100 --bottom 7100  --plot --wavelengths $wavelengthfile $fin #> ${fin}.txt

fin=/beegfs/scratch/drt/20230608_AVIRIS3_SRF/Col_919_SWIR_darksub_pedestal 
#python ../../utils/makesrf.py --top 100 --bottom 22100  --wavelengths $wavelengthfile $fin > ${fin}.txt

fin=/beegfs/scratch/drt/20230608_AVIRIS3_SRF/Col_351_VNIR_darksub_pedestal
python ../../utils/makesrf.py --top 100 --bottom 7100  --wavelengths $wavelengthfile $fin > ${fin}.txt
