# David R Thompson

# Create synthetic data
python synthesize_L1A.py emit_benchmarks_ang20170323t202244.json
python synthesize_L1A.py emit_benchmarks_ang20160218t083205.json
python synthesize_L1A.py emit_benchmarks_ang20170328t202059.json
                    
# Perform calibration
python ../emitrdn.py config_ang20160218t083205.json 
python ../emitrdn.py config_ang20170328t202059.json   
python ../emitrdn.py config_ang20170323t202244.json 
