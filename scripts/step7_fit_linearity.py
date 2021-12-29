# David R Thompson
import os, sys, glob

fieldpoints = [15,90,165,240,315,390,465,540,615,690,765,840,915,990,1065,1140,1215]
for fieldpoint in fieldpoints:
    cmd = 'python ~/src/emit-sds-l1b/utils/fitlinearity.py ' 
    fp = str(fieldpoint)
    if fieldpoint != fieldpoints[0]:
        cmd = cmd + '--draft ../data/EMIT_LinearityMap_20211215 '
    infiles = glob.glob('/beegfs/scratch/drt/20211115_EMIT_Linearity/20211117_023623_UTC_LinearitySphere/*Field'+fp+'*pedestal')
    infiles.sort()                                                                                                          
    for infile in infiles:
       cmd = cmd + ' '+infile
    cmd = cmd + ' ../data/EMIT_LinearityBasis_20211215 ../data/EMIT_LinearityMap_20211215'
    print(cmd)
    os.system(cmd)

