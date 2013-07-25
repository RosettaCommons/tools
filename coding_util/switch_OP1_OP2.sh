# if you are on a mac, change "-i" to "-i .old" below.
for i in "$@"; do 
 sed -i 's/OP1/OPX/g' $1
 sed -i 's/OP2/OP1/g' $1
 sed -i 's/OPX/OP2/g' $1
done


