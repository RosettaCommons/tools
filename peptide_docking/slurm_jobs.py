#!/usr/bin/python3

import os
import subprocess

import pfpd_const as pfpd

##########################################################################################
"""WE use SLURM workload manager. If you use something else - change the commands below"""
##########################################################################################

SBATCH_HEADER = '#!/bin/sh\n' \
               '#SBATCH --ntasks={ntasks}\n' \
               '#SBATCH --time=50:00:00\n' \
               '#SBATCH --get-user-env\n' \
               '#SBATCH --mem-per-cpu=1600m\n'

RUN_PIPER = ['sbatch', 'run_piper']
RUN_EXTRACT_DECOYS = ['sbatch', 'run_extract_decoys']
RUN_PREP_FPD_INPUTS = ['sbatch', 'run_prepare_fpd_inputs']
RUN_REFINEMENT = ['sbatch', 'run_refinement']
RUN_EXTRACT_TOP_MODEL = ['sbatch', 'extract_model']
RUN_RESCORING = ['sbatch', 'rescoring']
RUN_CLUSTERING = ['sbatch', 'run_clustering']


def send_piper_job(jobs_list, rec_name, lig_name, ppk_receptor, refinement_dir):
    """run PIPER docking, extract top 250 decoys and prepare input for refinement"""
    ####################################################################################
    # This function will execute:
    #
    # PIPER_BIN/piper -vv -c1.0 -k4 --msur_k=1.0 --maskr=1.0 -T FFTW_EXHAUSTIVE -R {decoys} -t 1
    # -p PIPER_DIR/prms/atoms.prms. -f PIPER_DIR/prms/coeffs.prms
    # -r PIPER_DIR/prms/rots.prm (receptor) (ligand)
    #
    # for f in `awk '{print $1}' ft.000.00 | head -250`; do if [ ! -f {f}.pdb ];
    # then PIPER_DIR/apply_ftresult.py -i $f ft.000.00 PIPER_DIR/prms/rots.prm (ligand) --out-prefix $f;fi;done
    #
    # piper_run=`pwd | awk -F'/' '{print $NF}'`
    # for f in `ls [0-9]*.pdb`;do cat (receptor.ppk) ${f} | grep ATOM > REFINEMENT_DIR/${piper_run}.${f};
    # gzip REFINEMENT_DIR/${piper_run}.${f};done
    #############################################

    with open('run_piper', 'w') as piper:
        piper.write(SBATCH_HEADER.format(ntasks=1))
        piper.write(pfpd.PIPER_DOCKING.format(decoys=pfpd.N_ROTS, r=rec_name, l=lig_name))
        piper.write(pfpd.EXTRACT_DECOYS % (pfpd.PIPER_MODELS_EXTRACTED, lig_name))
        piper.write(pfpd.PREP_FPD_INPUTS % (ppk_receptor, refinement_dir, refinement_dir))

        # run PIPER docking, extract top 250 decoys and prepare input for refinement
    piper_job = str(subprocess.check_output(RUN_PIPER))
    jobs_list.append(''.join(d for d in piper_job if d.isdigit()))

    return jobs_list


def send_fpd_job(fpd_command, jobs_list, refinement_dir):
    ###########################################
    # This function will execute:
    #
    # ls *gz >input_list
    # mpirun ROSETTA_BIN/FlexPepDocking.mpiserialization.linuxgccrelease -database ROSETTA_DB @flags
    ############################################

    with open(os.path.join(refinement_dir, 'run_refinement'), 'w') as refinement:
        refinement.write(SBATCH_HEADER.format(ntasks=300))
        refinement.write(fpd_command)

    all_prep_jobs = ''

    for job_id in jobs_list:
        all_prep_jobs += ':' + job_id  # there are multiple jobs that need to be separated by semicolon
    RUN_REFINEMENT.insert(1, '--dependency=aftercorr%s' % all_prep_jobs)
    run_refinement = str(subprocess.check_output(RUN_REFINEMENT))

    refinement_id = ''.join(d for d in run_refinement if d.isdigit())

    return refinement_id


def send_clustering_jobs(clustering, refinement_id, clustering_dir, refinement_dir, native=True, rescoring=None):
    ####################################################################################
    # This function will execute:
    #
    # python PFPD_SCRIPTS/extract_top_model.py
    # python PFPD_SCRIPTS/rescoring.py {sc_func} {rec}'
    # python PFPD_SCRIPTS/clustering.py 2.0 {native} {decoys}'
    ##########################################################

    with open(os.path.join(clustering_dir, 'run_clustering'), 'w') as cluster:
        cluster.write(SBATCH_HEADER.format(ntasks=1))
        cluster.write(clustering)

    if not native:
        with open(os.path.join(refinement_dir, 'rescoring'), 'w') as rescore:
            rescore.write(SBATCH_HEADER.format(ntasks=30))
            rescore.write(pfpd.EXTRACT_MODEL)
            rescore.write(rescoring)

        RUN_RESCORING.insert(1, '--dependency=aftercorr:%s' % refinement_id)
        rescoring = str(subprocess.check_output(RUN_RESCORING))

        rescoring_id = ''.join(d for d in rescoring if d.isdigit())
        os.chdir(clustering_dir)
        RUN_CLUSTERING.insert(1, '--dependency=aftercorr:%s' % rescoring_id)
        subprocess.call(RUN_CLUSTERING)
    else:

        os.chdir(clustering_dir)
        RUN_CLUSTERING.insert(1, '--dependency=aftercorr:%s' % refinement_id)
        subprocess.call(RUN_CLUSTERING)

