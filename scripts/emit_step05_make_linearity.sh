# David R Thompson

# collect linearity curves
for fieldpoint in 90 165 240 315 390 465 540 615 690 765 840 915 990 1065 1140 1215; do 
python ../utils/makelinearity.py /beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/*Field${fieldpoint}*clip*badfix /beegfs/scratch/drt/20211115_EMIT_Linearity/Field_${fieldpoint}_clipped
done

# Collect data points for an example plot at field point 315
for fieldpoint in 315; do
python ../utils/makelinearity.py --plot /beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/*Field${fieldpoint}*clip*badfix /beegfs/scratch/drt/20211115_EMIT_Linearity/Field_${fieldpoint}_clipped_test
done
