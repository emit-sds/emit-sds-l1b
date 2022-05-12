#!/home/drt/src/anaconda3/bin/python
# David R Thompson
import sys, os, os.path
print(sys.argv)

base_dir = os.path.split(os.path.abspath(__file__))[0]+'/../'

for i in range(1,len(sys.argv),2):

    if not '.img' in sys.argv[i]:
        raise FileNotFoundError('oops')

    if not '.img' in sys.argv[i+1]:
        raise FileNotFoundError('argh!')
        
    strip = sys.argv[i].replace('.img','_strip')
    cmd = 'python '+base_dir+'/utils/strip_bad.py %s %s' % (sys.argv[i],strip)
    print(cmd)
    os.system(cmd)

    shift = strip + '_shift'
    cmd = 'python '+base_dir+'/utils/leftshift.py %s %s' % (strip, shift)
    print(cmd)
    os.system(cmd)

    dark = shift + '_avg'
    cmd = 'python '+base_dir+'/utils/emit2dark.py %s %s' % (shift, dark)
    print(cmd)
    os.system(cmd)

    strip = sys.argv[i+1].replace('.img','_strip')
    cmd = 'python '+base_dir+'/utils/strip_bad.py %s %s' % (sys.argv[i+1], strip)
    print(cmd)
    os.system(cmd)

    shift = strip + '_shift'
    cmd = 'python '+base_dir+'/utils/leftshift.py %s %s' % (strip, shift)
    print(cmd)
    os.system(cmd)

    base = shift + '_darksub'
    cmd = 'python '+base_dir+'/utils/darksubtract.py %s %s %s' % (shift, dark, base)
    print(cmd)
    os.system(cmd)
      
    ped = base + '_pedestal'
    cmd = 'python '+base_dir+'/utils/pedestal.py %s %s' % (base,ped)
    print(cmd)
    os.system(cmd)

    badfix = ped + '_badfix'
    cmd = ('python '+base_dir+'/utils/fixbad.py %s '+base_dir+'/data/EMIT_BadElements_20220117 %s') % (ped,badfix)
    print(cmd)
    os.system(cmd)

    osffix = badfix + '_osffix'
    cmd = 'python '+base_dir+'/utils/fixosf.py %s %s' % (badfix,osffix)
    print(cmd)
    os.system(cmd)

    linear = osffix + '_linear'
    cmd = ('python '+base_dir+'/utils/fixlinearity.py %s '+base_dir+'/data/EMIT_LinearityBasis_20220117 '+base_dir+'/data/EMIT_LinearityMap_20220117 %s') % (osffix,linear)
    print(cmd)
    os.system(cmd)

    scatterfix = linear + '_scatterfix'
    cmd = ('python '+base_dir+'/utils/fixscatter.py %s '+base_dir+'/data/EMIT_SpatialScatter_20220117 '+base_dir+'/data/EMIT_SpectralScatter_20220117 %s') % (linear,scatterfix)
    print(cmd)
    os.system(cmd)
  
    ghostfix = scatterfix + '_ghostfix'
    cmd = 'python '+base_dir+'/utils/fixghostraster.py %s ~/src/emit-sds-l1b/data/EMIT_GhostMap_20220117.json %s' % (scatterfix,ghostfix)
    print(cmd)
    os.system(cmd)
  
