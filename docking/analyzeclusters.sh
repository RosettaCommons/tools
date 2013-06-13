ls cluster*pdb > clusterlist
rms2.pl clusterlist > clusters.rms
while read pdb
do
  findcontacts.py $pdb CI
  pdbid=$(echo $pdb | cut -c 1-8)
  sort -nr +7 $pdbid.contacts > tmp
  mv tmp $pdbid.contacts
done < clusterlist