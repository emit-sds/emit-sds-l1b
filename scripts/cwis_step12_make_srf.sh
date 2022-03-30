# David R Thompson

# This script relies on the caller to have first dark-subtracted and pedestal-shift-corrected all the SRF data

# Make bad pixels
for fin in `ls /beegfs/scratch/drt//20220112_CWIS2/20220110_SRF/*clip*pedestal`; do
   echo $fin
   python ../utils/makesrf.py --wavelengths ../data/CWIS_Wavelengths_20220203.txt --bottom 21500 --top_margin 500 --target_index 637 $fin # > ${fin}.txt
done


