# David R Thompson
import os, sys, glob

# This script was used when we were fitting linearity measurements to each FPA
# element.  It is no longer used.

#fieldpoints = [15,90,165,240,315,390,465,540,615,690,765,840,915,990,1065,1140,1215]
#for fieldpoint in fieldpoints:
#    cmd = 'python ../utils/fitlinearity.py ' 
#    fp = str(fieldpoint)
#    if fieldpoint != fieldpoints[0]:
#        cmd = cmd + '--draft ../data/EMIT_LinearityMap_20220117 '
#    infiles = glob.glob('/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/*Field'+fp+'*clip*pedestal')
#    infiles.sort()                                                                                                          
#    for infile in infiles:
#       cmd = cmd + ' '+infile
#    cmd = cmd + ' ../data/EMIT_LinearityBasis_20220117 ../data/EMIT_LinearityMap_20220117'
#    print(cmd)
#    os.system(cmd)

