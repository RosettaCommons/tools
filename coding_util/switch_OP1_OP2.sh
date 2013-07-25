#for i in cc hh params txt; do 
#for i in params; do 
 sed -i 's/OP1/OPX/g' $1
 sed -i 's/OP2/OP1/g' $1
 sed -i 's/OPX/OP2/g' $1
#done

#find ./ -name '*.old' | xargs rm
