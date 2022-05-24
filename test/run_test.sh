# David R Thompson

python ../emitrdn.py emit20220305t002601_subset_raw emit20220305t002444_subset_raw ../config/EMIT_20220504.json testout

if (diff reference testout)
then
   echo "PASS"
else
   echo "FAIL"
fi


