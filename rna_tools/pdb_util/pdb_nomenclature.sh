for i in "$@"; do 
    sed -i  's/O1P/OP2/g' $i
    sed -i  's/O2P/OP1/g' $i
    sed -i  "s/C1\*/C1\'/g" $i
    sed -i  "s/C2\*/C2\'/g" $i
    sed -i  "s/C3\*/C3\'/g" $i
    sed -i  "s/C4\*/C4\'/g" $i
    sed -i  "s/C5\*/C5\'/g" $i
    sed -i  "s/O2\*/O2\'/g" $i
    sed -i  "s/O3\*/O3\'/g" $i
    sed -i  "s/O4\*/O4\'/g" $i
    sed -i  "s/O5\*/O5\'/g" $i
    sed -i  "s/H1\*/H1\'/g" $i
    sed -i  "s/1H2\*/ H2\'/g" $i
    sed -i  "s/H3\*/H3\'/g" $i
    sed -i  "s/H4\*/H4\'/g" $i
    sed -i  "s/1H5\*/ H5\'/g" $i
    sed -i  "s/2H5\*/H5\'\'/g" $i
    sed -i  "s/2HO\*/HO2\'/g" $i

    sed -i  "s/C5M/C7 /g" $i
    sed -i  "s/1H5M/ H71/g" $i
    sed -i  "s/2H5M/ H72/g" $i
    sed -i  "s/3H5M/ H73/g" $i

    sed -i  "s/1H4 / H41/g" $i
    sed -i  "s/2H4 / H42/g" $i
    sed -i  "s/1H4/H41/g" $i
    sed -i  "s/2H4/H42/g" $i

    sed -i  "s/H21\*/H2''/g" $i
    sed -i  "s/H21'/H2''/g" $i
    
    sed -i  "s/2H2 / H21/g" $i
    sed -i  "s/1H2 / H22/g" $i
    sed -i  "s/2H2/H21/g" $i
    sed -i  "s/1H2/H22/g" $i
    
    sed -i "s/1H6 / H61/g" $i
    sed -i "s/2H6 / H62/g" $i
    sed -i "s/1H6/H61/g" $i
    sed -i "s/2H6/H62/g" $i

done

