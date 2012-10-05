#!/usr/bin/env phenix.python

import os.path
import imp

file_path = os.path.split( os.path.abspath(__file__) ) [0]
imp.load_source('erraser_util', file_path + '/erraser_util.py')
imp.load_source('measure_params', file_path + '/measure_params.py')

from measure_params import compute_torsion
from erraser_util import *
import glob

prefix = parse_options( sys.argv, 'prefix', '' )
res = parse_options( sys.argv, 'res', '' )
outfile = parse_options( sys.argv, 'out', 'erraser_single_res.analysis' )

data = []

####################################
def find_chi_angle_std_pdb( input_pdb, res_find ) :

    def find_coord( coords_list, atm_name ) :
        for i in coords_list :
            if atm_name == i[0] :
                return i[1]
        return []

    data_list = []
    current_res = ''
    current_res_name = ''
    res_coords = []
    for line in open(input_pdb) :
        if len(line) > 6 and line[0:4] == 'ATOM' :
            res_id = line[21:27].replace(' ', '')
            res_name = line[16:20]
            if current_res != res_id :
                if current_res != '' :
                    data_list.append( [current_res, current_res_name, res_coords] )
                res_coords = []
                current_res = res_id
                current_res_name = res_name
            atm_name = line[12:16].replace(' ', '')
            coord = []
            coord.append( float( line[30:38] ) )
            coord.append( float( line[38:46] ) )
            coord.append( float( line[46:54] ) )
            res_coords.append( [atm_name, coord] )

    purines = ['   A', ' ADE', ' A  ', '  rA', '   G', ' GUA', ' G  ', '  rG']
    pyrimidines = ['   U', ' URI', ' U  ', '  rU', '   C', ' CYT', ' C  ', '  rC']
    for res_coords in data_list :
        if res_find == res_coords[0] :
            res_name = res_coords[1]
            if res_name in purines :
                atom1 = find_coord(res_coords[2], "O4'")
                atom2 = find_coord(res_coords[2], "C1'")
                atom3 = find_coord(res_coords[2], "N9")
                atom4 = find_coord(res_coords[2], "C4")
            elif res_name in pyrimidines :
                atom1 = find_coord(res_coords[2], "O4'")
                atom2 = find_coord(res_coords[2], "C1'")
                atom3 = find_coord(res_coords[2], "N1")
                atom4 = find_coord(res_coords[2], "C2")
            else :
                continue
            if atom1 == [] or atom2 == [] or atom3 == [] or atom4 == [] :
                return 'NA'
            else :
                return compute_torsion(atom1, atom2, atom3, atom4)
    #return 'NA'
###############################
def syn_anti_check( chi ) :
    if chi <= 140 and chi > -40 :
        return " syn"
    else :
        return "anti"
####################################
def extract_info( pdb, res, name ) :
    subprocess_call( 'phenix.rna_validate outliers_only=False %s > rna_validate.temp' % pdb )
    subprocess_call( 'phenix.clashscore %s > clash.temp' % pdb )
    clash = 0
    suite_i = '__'
    suite_i_plus1 = '__'
    pucker_out = 'OK'
    for line in open('clash.temp') :
        if "clashscore =" in line :
            clash = float( line.split() [-1] )
            clash_out = '%.2f' % clash
            break

    current_entry = ''
    is_break_next = False
    for line in open('rna_validate.temp') :
        if "Pucker" in line :
            current_entry = 'pucker'
        elif "Bond" in line :
            current_entry = 'bond'
        elif "Angle" in line :
            current_entry = 'angle'
        elif "Suite" in line :
            current_entry = 'suite'

        if line[0] != '#' and ' :' in line :
            curr_res = line[5:11].replace(' ', '')
            if current_entry == 'pucker' :
                if curr_res == res :
                    pucker = float(line.split(':')[1])
                    pucker_out = '%4.0f/!!' % pucker
            elif current_entry == 'suite' :
                curr_suite = line[12:14]
                if curr_res == res :
                    suite_i = curr_suite
                    is_break_next = True
                    continue
                if is_break_next :
                    suite_i_plus1 = curr_suite
                    break
    chi = find_chi_angle_std_pdb( pdb, res )
    if chi == 'NA' :
        chi_out = chi
    else :
        chi_out = '%4.0f' % chi + '/' + syn_anti_check( chi )

    score = 0
    for line in open('scores.out') :
        if name in line :
            score = float(line.split() [-1])
    return [clash_out, suite_i, suite_i_plus1, pucker_out, chi_out, score]

####################################

pdb_list = sorted( glob.glob( prefix + '_*.pdb') )
pdb_id = 0
start_min_score = 0
for i in pdb_list :
    if 'start_min' in i :
        name = 'start_min'
        data.append( [name, extract_info(i, res, name)] )
        start_min_score = data[-1][1][-1]
    else :
        name = 'model_%d' % pdb_id
        pdb_id += 1
        data.append( [name, extract_info(i, res, name)] )


out = open(outfile, 'w')
out.write( "Changes Introduced by ERRASER Single-Res\n" )
out.write( "========================================================================\n" )
out.write( "           Clashscore  Suite_i  Suite_i+1   Pucker        Chi      Score\n" )
# "start_min      999.99       !!         !!  -100/!!  -100/anti  -999999.9"
for i in data :
    score = '%9.1f' % (i[1][-1] - start_min_score)
    out.write( '%9s%12s%9s%11s%9s%11s%11s\n' % (i[0], i[1][0], i[1][1], i[1][2], i[1][3], i[1][4], score) )
out.write( "========================================================================\n\n" )
out.write( "# ERRASER rebuilds two suites (i and i+1) and one pucker for each run.  \n" )
out.write( "# For outlier puckers (!!) the corresponding delta angle are shown.     \n" )
out.write( "# Chi intervals for syn: (-40, 140]. anti: (-180, -40] and (140, 180].  \n" )

remove('*.temp')
