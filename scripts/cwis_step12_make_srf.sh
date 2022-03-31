# David R Thompson

# This script relies on the caller to have first dark-subtracted and pedestal-shift-corrected all the SRF data

# Make SRFs
fin=/beegfs/scratch/drt//20220112_CWIS2/20220110_SRF/SRF_VNR_clip_darksub_pedestal
python ../utils/makesrf.py --piecewise_velocity --wavelengths ../data/CWIS_Wavelengths_20220203.txt --bottom 6800 --top_margin 500 --target_index 637 $fin  > ${fin}_637.txt
python ../utils/makesrf.py --piecewise_velocity --wavelengths ../data/CWIS_Wavelengths_20220203.txt --bottom 6800 --top_margin 500 --target_index 360 $fin  > ${fin}_360.txt
python ../utils/makesrf.py --piecewise_velocity --wavelengths ../data/CWIS_Wavelengths_20220203.txt --bottom 6800 --top_margin 500 --target_index 77 $fin  > ${fin}_77.txt
python ../utils/makesrf.py --piecewise_velocity --wavelengths ../data/CWIS_Wavelengths_20220203.txt --bottom 6800 --top_margin 500 --target_index 961 $fin  > ${fin}_961.txt
python ../utils/makesrf.py --piecewise_velocity --wavelengths ../data/CWIS_Wavelengths_20220203.txt --bottom 6800 --top_margin 500 --target_index 1211 $fin  > ${fin}_1211.txt

fin=/beegfs/scratch/drt//20220112_CWIS2/20220110_SRF/SRF_SWIR_clip_darksub_pedestal
python ../utils/makesrf.py --piecewise_velocity --wavelengths ../data/CWIS_Wavelengths_20220203.txt --bottom 21500 --top_margin 500 --target_index 77 $fin  > ${fin}_77.txt
python ../utils/makesrf.py --piecewise_velocity --wavelengths ../data/CWIS_Wavelengths_20220203.txt --bottom 21500 --top_margin 500 --target_index 360 $fin  > ${fin}_360.txt
python ../utils/makesrf.py --piecewise_velocity --wavelengths ../data/CWIS_Wavelengths_20220203.txt --bottom 21500 --top_margin 500 --target_index 637 $fin  > ${fin}_637.txt
python ../utils/makesrf.py --piecewise_velocity --wavelengths ../data/CWIS_Wavelengths_20220203.txt --bottom 21500 --top_margin 500 --target_index 961 $fin  > ${fin}_961.txt
python ../utils/makesrf.py --piecewise_velocity --wavelengths ../data/CWIS_Wavelengths_20220203.txt --bottom 21500 --top_margin 500 --target_index 1211 $fin  > ${fin}_1211.txt


