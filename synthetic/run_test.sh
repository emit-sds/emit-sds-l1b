# David R Thompson

# Create synthetic data
 python synthesize_L1A.py emit_benchmarks_ang20170323t202244.json
 python synthesize_L1A.py emit_benchmarks_ang20160218t083205.json
 python synthesize_L1A.py emit_benchmarks_ang20170328t202059.json 

 # Subset
 gdal_translate -of envi -co interleave=bil -srcwin 0 10000 1280 1000  \
                  /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000000_oOOOOO_sSSS_L0_raw_bBBB_vVV.img \
                  /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000003_oOOOOO_sSSS_L0_raw_bBBB_vVV.img
 gdal_translate -of envi -co interleave=bip -srcwin 0 10000 1280 1000  \
                  /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000000_oOOOOO_sSSS_L0_loc_bBBB_vVV.img \
                  /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000003_oOOOOO_sSSS_L0_loc_bBBB_vVV.img
 gdal_translate -of envi -co interleave=bip -srcwin 0 10000 1280 1000  \
                  /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000000_oOOOOO_sSSS_L0_obs_bBBB_vVV.img \
                  /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000003_oOOOOO_sSSS_L0_obs_bBBB_vVV.img
 cp /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000000_oOOOOO_sSSS_L0_dark_bBBB_vVV.img \
   /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000003_oOOOOO_sSSS_L0_dark_bBBB_vVV.img 
 cp /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000000_oOOOOO_sSSS_L0_dark_bBBB_vVV.img.hdr \
   /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000003_oOOOOO_sSSS_L0_dark_bBBB_vVV.img.hdr 

 gdal_translate -of envi -co interleave=bil -srcwin 0 10000 1280 0300  \
                  /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000000_oOOOOO_sSSS_L0_raw_bBBB_vVV.img \
                  /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000004_oOOOOO_sSSS_L0_raw_bBBB_vVV.img
 gdal_translate -of envi -co interleave=bip -srcwin 0 10000 1280 0300  \
                  /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000000_oOOOOO_sSSS_L0_loc_bBBB_vVV.img \
                  /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000004_oOOOOO_sSSS_L0_loc_bBBB_vVV.img
 gdal_translate -of envi -co interleave=bip -srcwin 0 10000 1280 0300  \
                  /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000000_oOOOOO_sSSS_L0_obs_bBBB_vVV.img \
                  /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000004_oOOOOO_sSSS_L0_obs_bBBB_vVV.img
 cp /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000000_oOOOOO_sSSS_L0_dark_bBBB_vVV.img \
   /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000004_oOOOOO_sSSS_L0_dark_bBBB_vVV.img 
 cp /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000000_oOOOOO_sSSS_L0_dark_bBBB_vVV.img.hdr \
   /beegfs/store/emit/drt/pge_validation_data/synthetic/emit20210110T000004_oOOOOO_sSSS_L0_dark_bBBB_vVV.img.hdr 

 # Perform calibration
 python ../emitrdn.py config_ang20160218t083205.json 
 python ../emitrdn.py config_ang20170328t202059.json   
 python ../emitrdn.py config_ang20170323t202244.json 
 python ../emitrdn.py config_ang20170328t202059-1000.json 
 python ../emitrdn.py config_ang20170328t202059-300.json 
