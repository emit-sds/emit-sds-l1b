#!/home/drt/src/anaconda3/bin/python
# David R Thompson
import sys, os
print(sys.argv)

for i in range(1,len(sys.argv),2):

    if not '.raw' in sys.argv[i]:
        raise FileNotFoundError('oops')

    if not '.raw' in sys.argv[i+1]:
        raise FileNotFoundError('argh!')
        
    dark = sys.argv[i].replace('.raw','_avg')
    cmd = 'python /home/drt/src/emit-sds-l1b/utils/emit2dark.py %s %s' % (sys.argv[i],dark)
    print(cmd)
    os.system(cmd)

    base = sys.argv[i+1].replace('.raw','_darksub')
    cmd = 'python /home/drt/src/emit-sds-l1b/utils/darksubtract.py %s %s %s' % (sys.argv[i+1],dark, base)
    print(cmd)
    os.system(cmd)
      
    ped = base + '_pedestal'
    cmd = 'python /home/drt/src/emit-sds-l1b/utils/pedestal.py %s %s' % (base,ped)
    print(cmd)
    os.system(cmd)

    badfix = ped + '_badfix'
    cmd = 'python /home/drt/src/emit-sds-l1b/utils/fixbad.py %s /home/drt/src/emit-sds-l1b/data/EMIT_Bad_Elements_20211229 %s' % (ped,badfix)
    print(cmd)
    os.system(cmd)

    osffix = badfix + '_osffix'
    cmd = 'python /home/drt/src/emit-sds-l1b/utils/fixosf.py %s %s' % (badfix,osffix)
    print(cmd)
    os.system(cmd)

    linear = osffix + '_linear'
    cmd = 'python /home/drt/src/emit-sds-l1b/utils/fixlinearity.py %s ~/src/emit-sds-l1b/data/EMIT_LinearityBasis_20211215 ~/src/emit-sds-l1b/data/EMIT_LinearityMap_20211215 %s' % (osffix,linear)
    print(cmd)
    os.system(cmd)

    scatterfix = linear + '_scatterfix'
    cmd = 'python /home/drt/src/emit-sds-l1b/utils/fixscatter.py %s ~/src/emit-sds-l1b/data/EMIT_SpatialScatter_20211226 ~/src/emit-sds-l1b/data/EMIT_SpectralScatter_20211226 %s' % (linear,scatterfix)
    print(cmd)
    os.system(cmd)
  
    ghostfix = scatterfix + '_ghostfix'
    cmd = 'python /home/drt/src/emit-sds-l1b/utils/fixghostpoints.py %s ~/src/emit-sds-l1b/data/EMIT_GhostMap_20211228.txt %s' % (scatterfix,ghostfix)
    print(cmd)
    os.system(cmd)
  
