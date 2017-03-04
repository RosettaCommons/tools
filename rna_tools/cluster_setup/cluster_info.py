#!/usr/bin/env python

from os.path import expanduser,expandvars,basename
import subprocess

user_name = basename( expanduser('~') )
biox3_user_name = expandvars( '$BIOX3_USER_NAME' )
if biox3_user_name == '$BIOX3_USER_NAME': biox3_user_name = user_name
sherlock_user_name = expandvars( '$SHERLOCK_USER_NAME' )
if not len(sherlock_user_name): sherlock_user_name = user_name
comet_user_name = expandvars( '$COMET_USER_NAME' )
if not len(comet_user_name): comet_user_name = user_name
xsede_user_name = expandvars( '$XSEDE_USER_NAME' )
xsede_dir_number = expandvars( '$XSEDE_DIR_NUMBER' )

def cluster_check( cluster_in ):
    clusterlist = [ 'syd','niau','seth','bes','hapy','apep','gebb','ptah','yah','isis','yah','maat','nut','fin','dig','biox2','biox2_scratch','biox3','biox3_scratch','vanlang_scratch','ade','ade.stanford.edu','steele','steele_scratch','tg-condor','tg-condor_scratch','abe','ncsa','abe_scratch','ade_scratch','vanlang','kwipapat','kwip','lovejoy','tsuname','lovejoy_scratch','backup','lonestar','ranger','lonestar_work','lonestar_scratch','trestles','stampede','stampede_scratch','sherlock', 'comet', 'sherlock','sherlock_scratch' ];

    cluster = cluster_in
    if cluster not in clusterlist:
        print 'Hey, '+cluster+' is not a known cluster.'
        cluster = 'unknown'

    cluster_dir = ''

    old_teragrid_user_name = 'dasr' # defunct

    if cluster == 'biox2': cluster = 'biox2.stanford.edu'
    if cluster == 'ade': cluster = '%s@ade.stanford.edu' % user_name
    if cluster == 'steele': cluster = '%s@tg-steele.purdue.teragrid.org' % old_teragrid_user_name
    if cluster == 'tg-condor': cluster ='%s@tg-condor.purdue.teragrid.org'  % old_teragrid_user_name
    if cluster == 'abe': cluster = '%s@login-abe.ncsa.teragrid.org' % xsede_user_name
    if cluster == 'ncsa': cluster ='%s@tg-login.ncsa.teragrid.org' % xsede_user_name
    if cluster == 'lonestar': cluster ='%s@lonestar.tacc.xsede.org' % xsede_user_name
    if cluster == 'trestles': cluster ='%s@trestles.sdsc.edu' % xsede_user_name
    if cluster == 'ranger': cluster ='%s@ranger.tacc.xsede.org' % xsede_user_name

    if cluster == 'backup':
        cluster = ''
        cluster_dir = '/Volumes/RhijuBackup/%s/' % user_name

    if cluster == 'steele_scratch': # defunct
        cluster = 'dasr@tg-steele.purdue.teragrid.org'
        cluster_dir = '/scratch/scratch95/d/%s/' % old_teragrid_user_name

    if cluster == 'stampede':
        cluster = '%s@stampede.tacc.xsede.org' % xsede_user_name
        cluster_dir = '/work/%s/%s/' % (xsede_dir_number, xsede_user_name)

    if cluster == 'stampede_scratch':
        cluster = '%s@stampede.tacc.xsede.org' % xsede_user_name
        cluster_dir = '/scratch/%s/%s/' % (xsede_dir_number, xsede_user_name)

    if cluster == 'lonestar_work':
        cluster ='%s@lonestar.tacc.xsede.org' % xsede_user_name
        cluster_dir = '/work/%s/%s/' % (xsede_dir_number, xsede_user_name)

    if cluster == 'lonestar_scratch':
        cluster ='%s@lonestar.tacc.xsede.org' % xsede_user_name
        cluster_dir = '/scratch/%s/%s/' % (xsede_dir_number, xsede_user_name)

    if cluster == 'tg-condor_scratch':
        cluster = 'dasr@tg-condor.purdue.teragrid.org'
        cluster_dir = '/scratch/scratch95/d/dasr/'

    if cluster == 'biox2_scratch':
        cluster = 'biox2.stanford.edu'
        cluster_dir = '/scratch/users/%s/' % user_name

    if cluster == 'biox3':
        cluster = '%s@biox3.stanford.edu' % biox3_user_name
        cluster_dir = '/home/%s/' % biox3_user_name

    if cluster == 'biox3_scratch':
        cluster = '%s@biox3.stanford.edu' % biox3_user_name
        cluster_dir = '/scratch/users/%s/' % biox3_user_name

    if cluster == 'sherlock':
        cluster = '%s@sherlock.stanford.edu' % biox3_user_name
        cluster_dir = '/home/%s/' % biox3_user_name

    if cluster == 'sherlock_scratch':
        cluster = '%s@sherlock.stanford.edu' % biox3_user_name
        cluster_dir = '/scratch/users/%s/' % biox3_user_name

    if cluster == 'ade_scratch':
        cluster = 'ade.stanford.edu'
        cluster_dir = '/scr/%s/' % user_name

    if cluster == 'ade':
        cluster = 'ade.stanford.edu'
        cluster_dir = '/home/%s/' % user_name

    if cluster == 'abe_scratch':
        cluster = '%s@login-abe.ncsa.teragrid.org' % xsede_user_name
        cluster_dir = '/scratch/users/%s/' % xsede_user_name

    if cluster == 'kwipapat':  cluster = 'kwipapat@biox3.stanford.edu'
    if cluster == 'kwip':      cluster = 'kwipapat@biox3.stanford.edu'
    if cluster == 'tsuname':   cluster = 'tsuname@biox3.stanford.edu'

    if cluster == 'sherlock':
        cluster = '%s@sherlock.stanford.edu' % sherlock_user_name
        cluster_dir = '/home/%s/' % sherlock_user_name

    if cluster == 'comet':
        cluster = '%s@comet.sdsc.xsede.org' % comet_user_name
        cluster_dir = '/home/%s/' % comet_user_name

    return (cluster,cluster_dir)


def strip_home_dirname( clusterdir ):
    clusterdir = clusterdir.replace('/Users/%s/' % user_name,'')
    clusterdir = clusterdir.replace('Dropbox/','')
    clusterdir = clusterdir.replace('/scratch/users/%s/' % user_name,'')
    clusterdir = clusterdir.replace('/work/%s/' % user_name,'')
    clusterdir = clusterdir.replace('/home/%s/' % user_name,'')
    clusterdir = clusterdir.replace('/home1/%s/%s/' % ( xsede_dir_number, xsede_user_name ),'')
    clusterdir = clusterdir.replace('/work/%s/%s/' % (xsede_dir_number, xsede_user_name ),'')
    return clusterdir
