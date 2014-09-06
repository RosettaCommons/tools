#!/usr/bin/python

from os.path import expanduser,expandvars,basename

user_name = basename( expanduser('~') )
xsede_user_name = expandvars( '$XSEDE_USER_NAME' )
xsede_dir_number = expandvars( '$XSEDE_DIR_NUMBER' )

def cluster_check( cluster_in ):
    clusterlist = [ 'syd','niau','seth','bes','hapy','apep','gebb','ptah','yah','isis','yah','maat','nut','fin','dig','biox2','biox2_scratch','biox3','biox3_scratch','vanlang_scratch','ade','ade.stanford.edu','steele','steele_scratch','tg-condor','tg-condor_scratch','abe','ncsa','abe_scratch','ade_scratch','vanlang','kwipapat','kwip','lovejoy','tsuname','lovejoy_scratch','backup','lonestar','ranger','lonestar_work','lonestar_scratch','trestles','stampede','stampede_scratch' ];

    cluster = cluster_in
    if cluster not in clusterlist:
        print 'Hey, '+cluster+' is not a known cluster.'
        cluster = 'unknown'

    cluster_dir = ''

    old_teragrid_user_name = 'dasr' # defunct

    if cluster == 'biox2': cluster = 'biox2.stanford.edu'
    if cluster == 'biox3': cluster = 'biox3.stanford.edu'
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

    if cluster == 'biox3_scratch':
        cluster = 'biox3.stanford.edu'
        cluster_dir = '/scratch/users/%s/' % user_name

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

    return (cluster,cluster_dir)

