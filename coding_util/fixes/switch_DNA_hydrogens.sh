# if you are on a mac, change "-i" to "-i .old" below.
for i in "$@"; do 
sed -i   "s/H2''/H2XX/g" $i
sed -i   "s/ H2'/H2''/g" $i
sed -i   "s/H2XX/ H2'/g" $i

sed -i   "s/H5''/H5XX/g" $i
sed -i   "s/ H5'/H5''/g" $i
sed -i   "s/H5XX/ H5'/g" $i

sed -i   "s/H41/H4X/g" $i
sed -i   "s/H42/H41/g" $i
sed -i   "s/H4X/H42/g" $i

done

