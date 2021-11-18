# David R Thompson

# Make bad pixels
for fin in `ls /beegfs/scratch/drt/20211113_EMIT_SRF/all/20211116*pedestal`; do
python ~/src/emit-sds-l1b/utils/makesrf.py $fin > ${fin}.txt
done


