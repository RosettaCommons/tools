# if you are on a mac, change "-i" to "-i .old" below.
for i in "$@"; do 

sed -i .old  "s/H41/H4X/g" $i
sed -i .old  "s/H42/H41/g" $i
sed -i .old "s/H4X/H42/g" $i

done

