for i in "$@"; do 
    sed -i  's/"  A"/" DA"/g' $i
    sed -i  's/"  C"/" DC"/g' $i
    sed -i  's/"  G"/" DG"/g' $i
    sed -i  's/"  T"/" DT"/g' $i
    sed -i  's/" rA"/"  A"/g' $i
    sed -i  's/" rC"/"  C"/g' $i
    sed -i  's/" rG"/"  G"/g' $i
    sed -i  's/" rU"/"  U"/g' $i
    sed -i 's/   A a/  DA a/g' $i
    sed -i 's/   C c/  DC c/g' $i
    sed -i 's/   G g/  DG g/g' $i
    sed -i 's/   T t/  DT t/g' $i
    sed -i 's/  rA a/   A a/g' $i
    sed -i 's/  rC c/   C c/g' $i
    sed -i 's/  rG g/   G g/g' $i
    sed -i 's/  rU u/   U u/g' $i
    #rm $i.old
done

