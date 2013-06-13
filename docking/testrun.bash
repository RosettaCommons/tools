#Mike Daily, 10/23/03
#CAUTION; THIS IS NOT THE CURE-ALL TESTRUN SCRIPT
#but it has most flags you will need for a testrun
#If you need a special flag (eg norepack, fab, randomize, dock_pert) you must provide this YOURSELF

pdb=$1

#Use the fix disulfide flag if necessary 
#The ppk.bash script should create one if necessary
#Then you can copy it over from startfiles in your prepacking directory.

if [ -f $pdb.fixdisulf ]
then
    disulfflags="-fix_disulf ../$pdb.fixdisulf -use_disulf_logfile $pdb.disulflog -norepack_disulf"
else
    disulfflags=""
fi

rrun.sh gnu $pdb aa -dock_mcm $disulfflags -chat -ex34 -ex1 -norepack2 -randomize1 -fake_native
