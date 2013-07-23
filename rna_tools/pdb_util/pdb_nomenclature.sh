#for i in cc hh params txt; do 
for i in pdb; do 
    find ./ -name \*.$i | xargs sed -i 's/O1P/OP2/g'
    find ./ -name \*.$i | xargs sed -i 's/O2P/OP1/g'
    find ./ -type f -name \*.$i | xargs sed -i "s/C1\*/C1\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/C2\*/C2\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/C3\*/C3\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/C4\*/C4\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/C5\*/C5\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/O2\*/O2\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/O3\*/O3\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/O4\*/O4\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/O5\*/O5\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/H1\*/H1\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/1H2\*/ H2\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/H3\*/H3\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/H4\*/H4\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/1H5\*/ H5\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/2H5\*/H5\'\'/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/2HO\*/HO2\'/g"

    find ./ -type f -name \*.$i | xargs sed -i "s/H21\*/H21\'/g"

    find ./ -type f -name \*.$i | xargs sed -i "s/1H4 / H41/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/2H4 / H42/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/1H4/H41/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/2H4/H42/g"
    
    find ./ -type f -name \*.$i | xargs sed -i "s/2H2 / H21/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/1H2 / H22/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/2H2/H21/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/1H2/H22/g"
    
    find ./ -type f -name \*.$i | xargs sed -i "s/1H6 / H61/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/2H6 / H62/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/1H6/H61/g"
    find ./ -type f -name \*.$i | xargs sed -i "s/2H6/H62/g"

    find ./ -type f -name \*.$i | xargs sed -i "s/H22'/H21'/g"

done

#find ./ -name '*.old' | xargs rm
