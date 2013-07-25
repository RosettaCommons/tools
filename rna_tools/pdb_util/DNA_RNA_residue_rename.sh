for i in "$@"; do 
    sed -i  .old 's/"  A"/" DA"/g' $i
    sed -i  .old 's/"  C"/" DC"/g' $i
    sed -i  .old 's/"  G"/" DG"/g' $i
    sed -i  .old 's/"  T"/" DT"/g' $i
    sed -i  .old 's/" rA"/"  A"/g' $i
    sed -i  .old 's/" rC"/"  C"/g' $i
    sed -i  .old 's/" rG"/"  G"/g' $i
    sed -i  .old 's/" rU"/"  U"/g' $i
    rm $i.old
done

