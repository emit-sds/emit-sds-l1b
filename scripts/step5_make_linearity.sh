# David R Thompson

# collect linearity curves
for fieldpoint in 15 90 165 240 315 390 465 540 615 690 765 840 915 990 1065 1140 1215; do 

  python ~/src/emit-sds-l1b/utils/makelinearity.py /beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/*Field${fieldpoint}*pedestal /beegfs/scratch/drt/20211115_EMIT_Linearity/Field_${fieldpoint}_curves

done