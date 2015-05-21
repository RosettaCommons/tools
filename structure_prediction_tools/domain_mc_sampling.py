#!/usr/local/bin/python2.7
from argparse import ArgumentParser
from multiprocessing import Pool
from os import popen, system
from os.path import basename, exists
from sys import exit, stderr, stdout
import random
from pprint import pprint
import operator
import math
from rosetta_util import read_silent_header
import cPickle as pickle

def density_score_reader( density_score_file, nullscore ):
    '''
    read into density score
    '''

    Dict = {}
    dens_top = {}

    with open( density_score_file ) as f:
        '''format :
            density_score after_rotation_domain.189.1TYQ_nohet_1.pdb.pdb rmsd'''

        for l in f:
            ls = l.split()
            placement     = ls[1]
            rmsd          = float( ls[2] )
            density_score = float( ls[0] ) * 0.1
            domain_id     = placement.split(".")[2].split("_")[2]

            if domain_id in Dict.keys():
                Dict[ domain_id ][ placement ] = ( density_score, rmsd )
                dens_top[ domain_id ] = max(density_score+nullscore, dens_top[ domain_id ] )
            else:
                Dict[ domain_id ] = { placement : ( density_score, rmsd ) }
                #Dict[ domain_id ][ "null" ] = ( nullscore, 0.0 ) # store a null placement for each domain id
                dens_top[ domain_id ] = density_score + nullscore

        for domain_id in Dict.keys():
            Dict[ domain_id ][ "null" ] = ( dens_top[domain_id], 100.0 )

    return Dict


def clash_score_reader( scorefile ):
    Dict = {}

    stderr.write( "reading %s" % scorefile )
    '''format A B clash_counts rmsd:
        after_rotation_domain.189.1TYQ_nohet_1.pdb.pdb after_rotation_domain.189.1TYQ_nohet_1.pdb.pdb clash_counts rmsd'''

    with open( scorefile ) as f:
        for l in f:
            ls = l.split()

            i_domain    = ls[0]
            j_domain    = ls[1]
            clash_score = float(ls[2])

            if i_domain in Dict.keys():
                Dict[ i_domain ][ j_domain ] = clash_score
            else:
                Dict[ i_domain ] = { j_domain : clash_score }

            # should go to both direction
            if j_domain in Dict.keys():
                Dict[ j_domain ][ i_domain ] = clash_score
            else:
                Dict[ j_domain ] = { i_domain : clash_score }

        stderr.write( "\ndone reading score table\n" )

    return Dict


def closab_score_reader( scorefile ):
    Dict = {}
    stderr.write( "reading %s" % scorefile )

    with open( scorefile ) as f:
        for l in f:
            ls = l.split()

            i_domain    = ls[0]
            j_domain    = ls[1]
            closab_score = float(ls[2])

            if i_domain in Dict.keys():
                Dict[ i_domain ][ j_domain ] = closab_score
            else:
                Dict[ i_domain ] = { j_domain : closab_score }

            # should go to both direction
            if j_domain in Dict.keys():
                Dict[ j_domain ][ i_domain ] = closab_score
            else:
                Dict[ j_domain ] = { i_domain : closab_score }

        stderr.write( "\ndone reading score table\n" )

    return Dict


def random_pick( list ):
    assert len(list) > 0
    item = list[ random.randrange( 0, len(list) ) ]
    return item


def normalization( Dict ):
    ''' normalized boltzmann Dict = { placed_domain : prob } '''

    normalized_Dict = {}
    factor = 1 / sum([ tuple[1][0] for tuple in Dict.items() ])

    for placed_domain in Dict.keys():
        normalized_Dict[ placed_domain ] = ( Dict[ placed_domain ][0]*factor, Dict[ placed_domain ][1] )

    if args.debug: print "normalized-factor:", factor
    if args.debug: pprint( normalized_Dict )

    return normalized_Dict


def boltzmann_prob( score, temperature ):
    '''
    my boltzmann
    '''
    #print temperature
    boltz_factor = ( -score ) / temperature
    boltzmann_prob = math.exp( boltz_factor )

    return boltzmann_prob


def read_scorefiles( args ):
    ## read density score
    dens_dict = density_score_reader( args.densityscore, args.nullscore )
    '''
    if not exists("density_score_Dict.pickle"):
        density_score_Dict = density_score_reader( args.densityscore, args.nullscore )
        pickle.dump( density_score_Dict, open("density_score_Dict.pickle", "w") )
    else:
        pkl = open("density_score_Dict.pickle", "rb" )
        density_score_Dict = pickle.load( pkl )
    '''

    ## read clash score
    if not exists("clash_score_Dict.pickle"):
        clash_dict = clash_score_reader( args.clashscore )
        pickle.dump( clash_dict, open("clash_score_Dict.pickle", "w") )
    else:
        pkl = open("clash_score_Dict.pickle", "rb" )
        clash_dict = pickle.load( pkl )

    ## read closab score
    if not exists("closab_score_Dict.pickle"):
        closab_dict = closab_score_reader( args.closabscore )
        pickle.dump( closab_dict, open("closab_score_Dict.pickle", "w") )
    else:
        pkl = open("closab_score_Dict.pickle", "rb" )
        closab_dict = pickle.load( pkl )

    #if args.debug:
        #pprint( density_score_Dict )
        #pprint( clash_score_Dict )

    return dens_dict, clash_dict, closab_dict



if __name__=='__main__':
    parser = ArgumentParser()
    parser.add_argument("-d", "--densityscore", required=True, help="")
    parser.add_argument("-dw", "--densityscore_weight", default=1.0, type=float, help="")

    parser.add_argument("-c", "--clashscore", required=True, help="")
    parser.add_argument("-cw", "--clashscore_weight", default=1.0, type=float, help="")

    parser.add_argument("-s", "--closabscore", required=True, help="")
    parser.add_argument("-sw", "--closabscore_weight", default=1.0, type=float, help="")

    parser.add_argument("-n", "--nullscore", default=-100.0, type=float, help="")
    parser.add_argument("--n_null_domain_allowed", default=2, type=int, help="")
    parser.add_argument("--log", required=True, help="")

    parser.add_argument("--sa_start_temp", default=50.0, type=float, help="")
    parser.add_argument("--sa_end_temp", default=1.0, type=float, help="")
    parser.add_argument("--sa_nsteps", default=25, type=int, help="")
    parser.add_argument("--mc_nsteps", default=10000, type=int, help="")

    #parser.add_argument("-n", "--norandom", type=int, help="")
    parser.add_argument("--models", default=1, type=int, help="")
    parser.add_argument("--seed", default=-1, type=int, help="")
    parser.add_argument("--debug", default=False, action="store_true", help="")
    args = parser.parse_args()

    density_score_Dict, clash_score_Dict, closab_score_Dict = read_scorefiles( args )
    if (args.seed>=0):
        print "Using constant random number seed"
        random.seed(args.seed)

    '''
    after the first round, every position has a domain assigned, now we can do Monte Carlo sampling
    domain_poses_Dict: is the highest-hierarchy strage
    '''
    logfile = open( args.log, "w" )
    temp_scale = ( args.sa_end_temp/args.sa_start_temp )**(1.0/args.sa_nsteps)
    print "time scale: ", temp_scale

    for each_model in range( args.models ):
        domain_poses_Dict = {}

        '''
        fill in random domains for each position for the first round - probably should do low-density score placement for each domain
        variable naming level:
            domain_id
            placed_domain
        '''
        ## store those numbersj
        domain_list = density_score_Dict.keys() # 1 2 3 4 5 6 7

        # assign intial placements by randomly pick one placement for each domain
        assign_dict = {}
        for domain_id in domain_list:
            placements_list = density_score_Dict[ domain_id ].keys()
            assign_dict[domain_id] = random_pick( placements_list )

        for domain_id in domain_list:
            placed_domain = assign_dict[domain_id]
            placed_domain_density_score = density_score_Dict[ domain_id ][ placed_domain ][0]
            placed_domain_rmsd = density_score_Dict[ domain_id ][ placed_domain ][1]
            selected_probility = 0.05
            total_score = 0.0 #for this domain

            all_closab_scores, all_clash_scores = ( 0.0, )*2

            for this_domain_id in domain_list:
                if this_domain_id == domain_id: continue #itself
                placed_domain_to_compare = assign_dict[ this_domain_id ]

                # clash score
                try:
                    clash_score = clash_score_Dict[ placed_domain ][ placed_domain_to_compare ]
                except KeyError:
                    if placed_domain_to_compare == "null" or placed_domain == "null":
                        clash_score = 0.0
                    else:
                        stderr.write("ERROR: domain %s/%s and %s/%s doesn't have clash score\n" %( domain_id, placed_domain, this_domain_id, placed_domain_to_compare ))
                        exit()
                all_clash_scores += clash_score*args.clashscore_weight

                # closab score
                if placed_domain_to_compare == "null" or placed_domain == "null":
                    closab_score = 0 # non-adjacent domains
                else:
                    closab_score = 0 # non-adjacent domains

                all_closab_scores += closab_score*args.closabscore_weight

            # end of evaluating all candidate placements at one domain
            total_score = placed_domain_density_score + all_clash_scores + all_closab_scores

            domain_poses_Dict[ domain_id ] = ( placed_domain,
                                             ( selected_probility, total_score ),
                                             ( placed_domain_density_score, placed_domain_rmsd ) )
        print domain_poses_Dict
        #exit()

        #print domain_poses_Dict
        mc_temp = args.sa_start_temp
        for sa_step in range( args.sa_nsteps ):
            # each temperature schedule
            for mc_step in range( args.mc_nsteps ):
                if args.debug:
                    print "model: %s  cycle: %s" %( each_model, mc_step )

                # pick a random domain_id to start optmize with
                domain_id = random_pick( domain_list )
                #print domain_id

                boltzmann_rank_prob_Dict = {}

                ''' assign a boltzmann prob for each domain at that position '''
                for placed_domain in density_score_Dict[ domain_id ].keys():
                    #print placed_domain
                    #calc allscore for each placement of this domain
                    density_score, all_closab_scores, all_clash_scores = ( 0.0, )*3

                    # one-body score:
                    try:
                        # this could be a null placement
                        density_score_Tuple = density_score_Dict[ domain_id ][ placed_domain ] # ( density_score, rmsd )
                        density_score = density_score_Tuple[0]
                        rmsd = density_score_Tuple[1]
                        #print density_score_Tuple

                    except KeyError:
                        stderr.write("ERROR: domain %s/%s doesn't have density score\n" %( domain_id, placed_domain ))
                        exit()

                    # two-body score: the domain to the rest
                    ''' calculate placed_domain to the selected ones for two body terms '''
                    # the rest of other domain_ids besides placed_domain
                    # calculate the one to consider and ones that have been selected
                    for this_domain_id in domain_list: #domain_poses_Dict.keys():
                        #print this_domain_id
                        if this_domain_id == domain_id: continue #itself
                        placed_domain_to_compare = domain_poses_Dict[ this_domain_id ][0]

                        # clash score
                        try:
                            clash_score = clash_score_Dict[ placed_domain ][ placed_domain_to_compare ]
                        except KeyError:
                            if placed_domain_to_compare == "null" or placed_domain == "null":
                                clash_score = 0.0
                            else:
                                stderr.write("ERROR: domain %s/%s and %s/%s doesn't have clash score\n" %( domain_id, placed_domain, this_domain_id, placed_domain_to_compare ))
                                exit()
                        all_clash_scores += clash_score*args.clashscore_weight

                        # closab score
                        try:
                            closab_score = closab_score_Dict[ placed_domain ][ placed_domain_to_compare ]
                        except KeyError:
                            if placed_domain_to_compare == "null" or placed_domain == "null":
                                closab_score = 0 # non-adjacent domains
                            else:
                                stderr.write("ERROR: domain %s/%s and %s/%s doesn't have closab score\n" %( domain_id, placed_domain, this_domain_id, placed_domain_to_compare ))
                                exit()
                        all_closab_scores += closab_score*args.closabscore_weight

                    # end of evaluating all candidate placements at one domain
                    total_scores = density_score + all_clash_scores + all_closab_scores
                    boltzmann_rank_prob_Dict[ placed_domain ] = ( boltzmann_prob( total_scores, mc_temp ), total_scores )

                    #after_rotation_domain.623.1TYQ_nohet_3.pdb': (0.9578719689882209, -70.2713)
                    #if (placed_domain=="after_rotation_domain.623.1TYQ_nohet_3.pdb"):
                    #    if (domain_poses_Dict[ this_domain_id ][0] == "after_rotation_domain.1606.1TYQ_nohet_6.pdb"):
                    #        print all_clash_scores
                    #        if (total_scores == -70.2713):
                    #            print density_score, all_clash_scores, all_closab_scores
                    #            print sa_step, mc_step
                    #            exit()

                    if args.debug:
                        print "domain %s/%s\trmsd: %5.4f\ttotal_score: %5.4f\tdensity_score: %5.4f\tclash_score: %5.4f\tclosab_Score: %5.4f\n\nboltzmann: %s" %(
                            domain_id,
                            placed_domain,
                            rmsd,
                            total_scores,
                            density_score,
                            all_clash_scores,
                            all_closab_scores,
                            boltzmann_rank_prob_Dict[ placed_domain ])

                ''' random_prob from (0,1) - candidate_placed_domain at that position '''
                # boltzmann_Dict = ( prob, total_scores )
                try:
                    normalized_boltzmann_rank_prob_Dict = normalization( boltzmann_rank_prob_Dict )
                except ZeroDivisionError: # this means clash dominates, need to optimize other domains
                    continue

                random_prob = random.uniform( 0, 1 )
                for candidate_placed_domain in normalized_boltzmann_rank_prob_Dict.keys():
                    random_prob = random_prob - normalized_boltzmann_rank_prob_Dict[ candidate_placed_domain ][0]
                    if args.debug:
                        print candidate_placed_domain, normalized_boltzmann_rank_prob_Dict[ candidate_placed_domain ]

                    selected_placed_domain = candidate_placed_domain
                    if random_prob < 0:
                        break

                if args.debug:
                    print "This domain: %s %s %s got picked!!\n" %( selected_placed_domain, normalized_boltzmann_rank_prob_Dict[ selected_placed_domain ], density_score_Dict[ domain_id ][ selected_placed_domain ] )

                domain_poses_Dict[ domain_id ] = ( selected_placed_domain,
                                                  normalized_boltzmann_rank_prob_Dict[ selected_placed_domain ],
                                                  density_score_Dict[ domain_id ][ selected_placed_domain ] )

                ## theoretically other domains need to be re-assigned but, due to the mc protocol here, it is not necessary,
                ## when the domain with wrong two-body score is selected, score will be updated
                ## if it is not selected, wrong score won't affect bolzman weight any way, so we only need to recalculate energy at the end

                #exit()
                if args.debug:
                    pprint( domain_poses_Dict )

            mc_temp = temp_scale*mc_temp
            #print mc_temp

        #redo all score
        for domain_id in domain_list:
            placed_domain = domain_poses_Dict[domain_id][0]
            placed_domain_density_score = density_score_Dict[ domain_id ][ placed_domain ][0]
            placed_domain_rmsd = density_score_Dict[ domain_id ][ placed_domain ][1]
            selected_probility = domain_poses_Dict[domain_id][1][0]
            total_score = 0.0 #for this domain

            all_closab_scores, all_clash_scores = ( 0.0, )*2

            for this_domain_id in domain_list:
                if this_domain_id == domain_id: continue #itself
                placed_domain_to_compare = domain_poses_Dict[this_domain_id][0]

                # clash score
                try:
                    clash_score = clash_score_Dict[ placed_domain ][ placed_domain_to_compare ]
                except KeyError:
                    if placed_domain_to_compare == "null" or placed_domain == "null":
                        clash_score = 0.0
                    else:
                        stderr.write("ERROR: domain %s/%s and %s/%s doesn't have clash score\n" %( domain_id, placed_domain, this_domain_id, placed_domain_to_compare ))
                        exit()
                all_clash_scores += clash_score*args.clashscore_weight

                # closab score
                if placed_domain_to_compare == "null" or placed_domain == "null":
                    closab_score = 0 # non-adjacent domains
                else:
                    closab_score = 0 # non-adjacent domains

                all_closab_scores += closab_score*args.closabscore_weight

            # end of evaluating all candidate placements at one domain
            total_score = placed_domain_density_score + all_clash_scores + all_closab_scores

            domain_poses_Dict[ domain_id ] = ( placed_domain,
                                             ( selected_probility, total_score ),
                                             ( placed_domain_density_score, placed_domain_rmsd ) )


        print args.log, "model_", each_model
        pprint( domain_poses_Dict )
        total_score         = 0.0
        total_rmsd_addup    = 0.0
        total_density_score = 0.0
        low_rmsd            = 0.0
        n_null_domain = 0

        outlines = ""
        for domain_id in domain_list: #domain_poses_Dict.keys():
            placed_domain = domain_poses_Dict[ domain_id ][0]
            allscore      = domain_poses_Dict[ domain_id ][1][1] #fix a bug here
            rmsd          = domain_poses_Dict[ domain_id ][2][1]
            density_score = domain_poses_Dict[ domain_id ][2][0]
            pair_score = allscore-density_score

            outlines += "model: %d,  domain_id: %s,  allscore %f,  rmsd: %5.4f,   %s,  density_score: %f\n" \
                        %( each_model, domain_id, allscore, rmsd, placed_domain, density_score )

            total_rmsd_addup    += rmsd
            total_density_score += density_score
            total_score         += (density_score + pair_score/2.0)

            if placed_domain == "null":
                n_null_domain += 1

        # control n_null_domain allowed
        if n_null_domain > args.n_null_domain_allowed:
            continue

        outlines +=  "model: %d,  total_score: %f, total_rmsd: %f, total_density_score: %f, average_rmsd: %f\n" \
                     %( each_model, total_score, total_rmsd_addup, total_density_score, total_rmsd_addup/len(domain_poses_Dict.keys()) )

        logfile.write( outlines )
        logfile.flush()
    logfile.close()

