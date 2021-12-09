# David R Thompson

# collect new linearity curves
#for fieldpoint in 15 90 165 240 315 390 465 540 615 690 765 840 915 990 1065 1140 1215; do 
#for fieldpoint in 690 765 840 915 990 1065 1140 1215; do 

# python ~/src/emit-sds-l1b/utils/makelinearity.py --plot /beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/*Field${fieldpoint}*linear /beegfs/scratch/drt/20211115_EMIT_Linearity/Field_${fieldpoint}_linearcheck

#done


# Matador test!

python ../utils/matador.py \
/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/20211117_025005_UTC_LinearitySphere_Field1215_Step14p8mm_PD1289p0candelam2_darksub_pedestal_badfix \
/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/20211117_024003_UTC_LinearitySphere_Field1215_Step18p4mm_PD6p0candelam2_darksub_pedestal_badfix \
/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/20211117_025005_UTC_LinearitySphere_Field1215_Step14p8mm_PD1289p0candelam2_darksub_pedestal_badfix_linear \
/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/20211117_024003_UTC_LinearitySphere_Field1215_Step18p4mm_PD6p0candelam2_darksub_pedestal_badfix_linear
