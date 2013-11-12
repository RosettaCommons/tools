# if you are on a mac, change "-i" to "-i .old" below.
for i in "$@"; do 

sed -i  .old "s/H21/H2X/g" $i
sed -i  .old "s/H22/H21/g" $i
sed -i  .old "s/H2X/H22/g" $i

done

