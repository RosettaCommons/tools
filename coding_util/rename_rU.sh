for i in "$@"; do 
    sed -i  .old 's/" rU"/"  U"/g' $i
    rm $i.old
done

