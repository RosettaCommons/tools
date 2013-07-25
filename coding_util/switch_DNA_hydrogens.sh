for i in "$@"; do 
sed -i .old  "s/H2''/H2XX/g" $i
sed -i .old  "s/ H2'/H2''/g" $i
sed -i .old  "s/H2XX/ H2'/g" $i

sed -i .old  "s/H5''/H5XX/g" $i
sed -i .old  "s/ H5'/H5''/g" $i
sed -i .old  "s/H5XX/ H5'/g" $i

sed -i .old  "s/H41/H4X/g" $i
sed -i .old  "s/H42/H41/g" $i
sed -i .old  "s/H4X/H42/g" $i

done

