#!/bin/bash

#Analysis script for repertoire metrics. Builds a large text file with SASA, hSASA, Charge repertoire metrics associated. Paths will need to be altered for proper processing. A list of repertoire-wide gene assignemts of each antibody isalso highly recommended for completeness. Gene assignment inclusion requires uncommenting lines 97-101.


if [ -z "$1" ]
then
    echo Usage: $0 repertoire_directory [force]
    echo H3 jobs processing
    exit
fi
#path
repertoire=$1
cd $repertoire
rm Antibody_Id
rm Antibody_repertoire_by_H3
rm Antibody_repertoire_output
rm CDR_charge
rm CDR_SASA
rm CDR_SASA_HP
rm CPUs
rm Decoy
rm H1
rm H1_charge
rm H1_SASA
rm H1_SASA_HP
rm H2
rm H2_charge
rm H2_SASA
rm H2_SASA_HP
rm H3
rm H3_charge
rm H3_SASA
rm H3_SASA_HP
rm L1
rm L1_charge
rm L1_SASA
rm L1_SASA_HP
rm L2
rm L2_charge
rm L2_SASA
rm L2_SASA_HP
rm L3
rm L3_charge
rm L3_SASA
rm L3_SASA_HP
rm IGHD
rm IGHJ
rm IGHV
rm IGLJ
rm IGLV

dir=$(ls -v -dp *| grep '/' | sed 's:/::')

#loop over repertoire
for d in $dir; do                                      
    cd $d
    pwd
    echo $d >> ../Antibody_Id
    if [ ! -s "remodel_h3.fasc.sort.top10" ]
    then
	#sort score_flie by chainbreak, score
	head -2 remodel_h3.fasc | tail -1 > remodel_h3.fasc.sort
	sort -n -k  2 remodel_h3.fasc |grep -vE 'SEQUENCE|H1_RMS' >> remodel_h3.fasc.sort
	head -2 remodel_h3.fasc | tail -1 > remodel_h3.fasc.chainbreak
	sort -n -k 11 remodel_h3.fasc |grep -vE 'SEQUENCE|H1_RMS' >> remodel_h3.fasc.chainbreak
       
        #top 10 models
	head -11 remodel_h3.fasc.sort | awk '{print $1, $2, $42}' > top10.scores
	head -11 remodel_h3.fasc.sort | tail | awk '{print $42}' > top10.models
	head -11 remodel_h3.fasc.sort > remodel_h3.fasc.sort.top10
    else
	cat remodel_h3.fasc.sort.top10 |  awk '{print $1, $2, $42}' > top10.scores
	cat remodel_h3.fasc.sort.top10 |  tail | awk '{print $42}' > top10.models
    fi

    if [ ! -d top10 ] 
    then
	mkdir top10
	i=0
	for x in $(cat top10.models); do
	    (( i += 1 ))
	    cp pdbs/$x.pdb.gz top10/$d-m$i.pdb.gz
	done
	gzip -d top10/*.pdb.gz
    fi	

    #Grab H, L length data and paratope data
    if [ -s "remodel_h3.fasc.sort.top10" ] && [ -s "grafting.log" ]
    then
	grep "H3 detected" grafting.log | awk '{print $5}' >> ../H3
	grep "H1 detected" grafting.log | awk '{print $5}' >> ../H1
	grep "H2 detected" grafting.log | awk '{print $5}' >> ../H2
	grep "L1 detected" grafting.log | awk '{print $5}' >> ../L1
	grep "L2 detected" grafting.log | awk '{print $5}' >> ../L2
	grep "L3 detected" grafting.log | awk '{print $5}' >> ../L3
	#gep "IGH" $d"HL.fasta" | awk -F " " '{ for(j=1;j<=NF;j++) if ($j == "[V:") print $(j+1) }' >> ../IGHV
	#gep "IGH" $d"HL.fasta" | awk -F " " '{ for(j=1;j<=NF;j++) if ($j == "[D:") print $(j+1) }' >> ../IGHD
	#grep "IGH" $d"HL.fasta" | awk -F " " '{ for(j=1;j<=NF;j++) if ($j == "[J:") print $(j+1) }' >> ../IGHJ
	#grep -E "IGL|IGK" $d"HL.fasta" | awk -F " " '{ for(j=1;j<=NF;j++) if ($j == "[V:") print $(j+1) }' >> ../IGLV
	#grep -E "IGL|IGK" $d"HL.fasta" | awk -F " " '{ for(j=1;j<=NF;j++) if ($j == "[J:") print $(j+1) }' >> ../IGLJ
	
	if [ ! -s "metrics_out.log" ]
	then
	    /Rosetta/main/source/bin/antibody_CR_paratope.linuxiccrelease -in:file:s top10/*-m1.pdb >metrics_out.log
	    fi
	    
    else
	echo --- >> ../H3
	echo --- >> ../H1
	echo --- >> ../H2
	echo --- >> ../L1
	echo --- >> ../L2
	echo --- >> ../L3
	echo --- >> ../IGHV
        echo --- >> ../IGHD
        echo --- >> ../IGHJ
        echo --- >> ../IGLV
        echo --- >> ../IGLJ
      	fi
    
    #CPU time
    list=$(ls outerr/* | grep -v "out_0")
    seconds=0
    decoys=0
    for x in $list; do
	if grep -q "no more batches to process..." "$x"
        then
            time=$(tail $x | sed -n '9p' | awk '{print $10}')
            decoy_num=$(tail $x | sed -n '9p' | awk '{print $6}')
            seconds=$(( seconds+time ))
            decoys=$(( $decoys+$decoy_num ))
	fi
    done

    cpu=$( echo "scale=3; $seconds/3600" | bc )
    if [ $decoys!=0 ]
    then
	abdecoytime=$( echo "scale=3; $seconds/$decoys" | bc )
    else 
	abdecoytime= ---
	cpu= ---
    fi
    echo $cpu >> ../CPUs
    echo $cpu" cpu"
    echo $abdecoytime >> ../Decoy
    echo $abdecoytime " decoy"


    if [ -s "metrics_out.log" ] && grep -q "Paratope charge:" metrics_out.log
    then
	grep "Paratope_SASA" metrics_out.log | awk '{print $4}' >> ../CDR_SASA
	grep "Paratope_hSASA" metrics_out.log | awk '{print $4}' >> ../CDR_SASA_HP
	grep "Paratope charge:"  metrics_out.log | awk '{print $5}' >> ../CDR_charge
        grep "Loop H1: CDR_SASA" metrics_out.log | awk '{print $6}' >> ../H1_SASA
	grep "Loop H2: CDR_SASA" metrics_out.log | awk '{print $6}' >> ../H2_SASA
	grep "Loop H3: CDR_SASA" metrics_out.log | awk '{print $6}' >> ../H3_SASA
	grep "Loop L1: CDR_SASA" metrics_out.log | awk '{print $6}' >> ../L1_SASA
	grep "Loop L2: CDR_SASA" metrics_out.log | awk '{print $6}' >> ../L2_SASA
	grep "Loop L3: CDR_SASA" metrics_out.log | awk '{print $6}' >> ../L3_SASA
	grep "Loop H1: CDR_SASA" metrics_out.log | awk '{print $8}' >> ../H1_SASA_HP
	grep "Loop H2: CDR_SASA" metrics_out.log | awk '{print $8}' >> ../H2_SASA_HP
        grep "Loop H3: CDR_SASA" metrics_out.log | awk '{print $8}' >> ../H3_SASA_HP
        grep "Loop L1: CDR_SASA" metrics_out.log | awk '{print $8}' >> ../L1_SASA_HP
        grep "Loop L2: CDR_SASA" metrics_out.log | awk '{print $8}' >> ../L2_SASA_HP
        grep "Loop L3: CDR_SASA" metrics_out.log | awk '{print $8}' >> ../L3_SASA_HP
	grep "Loop H1:" metrics_out.log | grep "anions" | awk '{print $10}' >> ../H1_charge
        grep "Loop H2:" metrics_out.log | grep "anions" | awk '{print $10}' >> ../H2_charge
        grep "Loop H3:" metrics_out.log | grep "anions" | awk '{print $10}' >> ../H3_charge
        grep "Loop L1:" metrics_out.log | grep "anions" | awk '{print $10}' >> ../L1_charge
        grep "Loop L2:" metrics_out.log | grep "anions" | awk '{print $10}' >> ../L2_charge
        grep "Loop L3:" metrics_out.log | grep "anions" | awk '{print $10}' >> ../L3_charge
    else
	echo --- >> ../CDR_SASA
	echo --- >> ../CDR_SASA_HP
	echo --- >> ../CDR_charge
	echo --- >> ../H1_SASA
        echo --- >> ../H2_SASA
        echo --- >> ../H3_SASA
        echo --- >> ../L1_SASA
        echo --- >> ../L2_SASA
        echo --- >> ../L3_SASA
        echo --- >> ../H1_SASA_HP
        echo --- >> ../H2_SASA_HP
        echo --- >> ../H3_SASA_HP
        echo --- >> ../L1_SASA_HP
        echo --- >> ../L2_SASA_HP
        echo --- >> ../L3_SASA_HP
        echo --- >> ../H1_charge
        echo --- >> ../H2_charge
        echo --- >> ../H3_charge
        echo --- >> ../L1_charge
        echo --- >> ../L2_charge
        echo --- >> ../L3_charge
        fi
    cd ../
done

sed -i "1iAntibody" Antibody_Id
sed -i "1iCDR_SASA" CDR_SASA
sed -i "1iCDR_Hydrophobic_SASA" CDR_SASA_HP
sed -i "1iCDR_Charge" CDR_charge
sed -i "1iH1_Length" H1
sed -i "1iH2_Length" H2
sed -i "1iH3_Length" H3
sed -i "1iL1_Length" L1
sed -i "1iL2_Length" L2
sed -i "1iL3_Length" L3
sed -i "1iH1_SASA" H1_SASA
sed -i "1iH2_SASA" H2_SASA
sed -i "1iH3_SASA" H3_SASA
sed -i "1iL1_SASA" L1_SASA
sed -i "1iL2_SASA" L2_SASA
sed -i "1iL3_SASA" L3_SASA
sed -i "1iH1_HP_SASA" H1_SASA_HP
sed -i "1iH2_HP_SASA" H2_SASA_HP
sed -i "1iH3_HP_SASA" H3_SASA_HP
sed -i "1iL1_HP_SASA" L1_SASA_HP
sed -i "1iL2_HP_SASA" L2_SASA_HP
sed -i "1iL3_HP_SASA" L3_SASA_HP
sed -i "1iH1_Charge" H1_charge
sed -i "1iH2_Charge" H2_charge
sed -i "1iH3_Charge" H3_charge
sed -i "1iL1_Charge" L1_charge
sed -i "1iL2_Charge" L2_charge
sed -i "1iL3_Charge" L3_charge
sed -i "1iIGHV" IGHV
sed -i "1iIGHD" IGHD
sed -i "1iIGHJ" IGHJ
sed -i "1iIGLV" IGLV
sed -i "1iIGLJ" IGLJ
sed -i "1iCPU_hrs" CPUs
sed -i "1iDecoy_time" Decoy

paste Antibody_Id H1 H2 H3 L1 L2 L3 CDR_SASA H1_SASA H2_SASA H3_SASA L1_SASA L2_SASA L3_SASA CDR_SASA_HP H1_SASA_HP H2_SASA_HP H3_SASA_HP L1_SASA_HP L2_SASA_HP L3_SASA_HP CDR_charge H1_charge H2_charge H3_charge L1_charge L2_charge L3_charge IGHV IGHD IGHJ IGLV IGLJ CPUs Decoy | column -t > Antibody_output

#sort -n -k4 Antibody_repertoire_output >Antibody_repertoire_by_H3 