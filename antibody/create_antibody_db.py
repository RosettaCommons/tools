#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   antibody.py
## @brief  Script to create antibody.db file from old info/* files
## @author Sergey Lyskov

import os, sys, time, os.path, urllib2, httplib, socket, json


def read_formated_text_file(file_name):
    ''' return dict with first field as key and value equal to dict(field_name=value)
    '''
    info, legend = {}, ''
    for l in file(file_name):
        if l.startswith('# '): legend = l[2:].split()
        elif len(l)>8: info[l.split()[0]] =  dict( zip(legend, l.split() ) )

    return info, legend


def write_formated_text_file(file_name, data, keys):
    ''' Create plain text file ising data dict and fill it with keys columns
    '''

    with file(file_name, 'w') as f:

        column_lengths = {}
        for k in keys:
            l = len('# '+k) if k==keys[0] else len(k)  # adjusting for legend prefix
            for p in data:
                if k in data[p]: l = max(l, len(data[p][k]) )
            column_lengths[k] = l


        f.write(' '.join(['{0:^{1}}'.format('# '+k if k==keys[0] else k, column_lengths[k]) for k in keys]) + '\n')

        sorted_pdbs = sorted(data.keys()) #, key=lambda x: x.partition('_chothia.pdb')[0][3:] )

        for p in sorted_pdbs:
            value = data.get(k, '-')

            f.write( ' '.join([ '{0:^{1}}'.format(data[p].get(k, '-'), column_lengths[k]) for k in keys]) + '\n')



def donwload_fasta_files(pdb_codes):
    ''' Donwload fasta files from PDB web sites and return map pdb_code:file_name
    '''
    fasta_dir = '._temp_'
    if not os.path.isdir(fasta_dir): os.mkdir( fasta_dir )

    files = {}

    for pdb in pdb_codes:
        url = 'http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=FASTA&compression=NO&structureId='+pdb.capitalize()
        file_name = fasta_dir + '/' + pdb + '.fasta'

        if not os.path.isfile(file_name):

            print 'Donwloading {} from {}...'.format(file_name, url)

            while True:
                try:
                    u = urllib2.urlopen( urllib2.Request(url) )
                    with file(file_name, 'w') as f: f.write(u.read())
                    break

                except (urllib2.HTTPError, urllib2.URLError, httplib.BadStatusLine, socket.error) as e: pass

                print 'ERROR: Download failed! Going to sleep and retry...'
                time.sleep(16)



        files[pdb] = file_name

    return files



def parse_fasta_file(file_name):
    ''' Parse FASTA file and return list of sequences it contain
    '''
    sequences = []
    with file(file_name) as f:
        for line in f:
            if line.startswith('>'): sequences.append('')
            else: sequences[-1] += line[:-1]  # removing end-line

    return sequences


def detect_fr_components(sequence, cdr1, cdr2, cdr3, fr):
    ''' Tries to tedect individual framework component and return them as tuple. Return empty string on failure
    '''
    fr1_end = sequence.find(cdr1)
    fr2_end = sequence.find(cdr2)
    fr3_end = sequence.find(cdr3)

    if fr1_end < 0  or  fr2_end < 0  or  fr3_end < 0: return '', '', '', ''

    fr2_begin = fr1_end + len(cdr1)
    fr3_begin = fr2_end + len(cdr2)
    fr4_begin = fr3_end + len(cdr3)

    fr2 = sequence[fr2_begin:fr2_end]
    fr3 = sequence[fr3_begin:fr3_end]

    print 'sequence:', sequence
    print 'cdr1:', cdr1
    print 'cdr2:', cdr2
    print 'cdr3:', cdr3

    print 'fr:', fr
    print 'fr2:', fr2
    print 'fr3:', fr3

    fr1, maybe_fr2, _ = fr.partition(fr2)
    _, maybe_fr3, fr4 = fr.partition(fr3)

    if maybe_fr2 != fr2: print 'detect_fr_components failed! maybe_fr2 != fr2'; #sys.exit(1)
    if maybe_fr3 != fr3: print 'detect_fr_components failed! maybe_fr3 != fr3'; #sys.exit(1)

    if fr1+fr2+fr3+fr4 != fr: print 'detect_fr_components failed! fr1+fr2+fr3+fr4 != fr'; #sys.exit(1)

    return fr1, fr2, fr3, fr4


def detect_components(cdr_info):
    fasta_files = donwload_fasta_files( cdr_info.keys() )
    #print 'FASTA:', parse_fasta_file('._temp_/1a7q.fasta')

    for pdb in cdr_info:
        #print 'PDB:', pdb
        sequences = parse_fasta_file( fasta_files[pdb] )

        h1 = cdr_info[pdb]['H1']
        h2 = cdr_info[pdb]['H2']
        h3 = cdr_info[pdb]['H3']
        frh = cdr_info[pdb]['FRH'] if 'FRH' in cdr_info[pdb] else '-'

        l1 = cdr_info[pdb]['L1']
        l2 = cdr_info[pdb]['L2']
        l3 = cdr_info[pdb]['L3']
        frl = cdr_info[pdb]['FRL'] if 'FRL' in cdr_info[pdb] else '-'

        if h1 == '-': h1=''
        if h2 == '-': h2=''
        if h3 == '-': h3=''
        if frh == '-': frh=''

        if l1 == '-': l1=''
        if l2 == '-': l2=''
        if l3 == '-': l3=''
        if frl == '-': frl=''

        heavy = ''
        if h1 or h2 or h3:
            for s in sequences:
                if s.find(h1) >= 0  and  s.find(h2) >= 0  and  s.find(h3) >= 0: heavy=s; break

        light = ''
        if l1 or l2 or l3:
            for s in sequences:
                if s.find(l1) >= 0  and  s.find(l2) >= 0  and  s.find(l3) >= 0: light=s; break

        # frh1, frh2, frh3, frh4 = '', '', '', ''
        # if h1 and h2 and h3 and heavy and frh:  frh1, frh2, frh3, frh4 = detect_fr_components(heavy, h1, h2, h3, frh)

        # frl1, frl2, frl3, frl4 = '', '', '', ''
        # if l1 and l2 and l3 and light and frl:  frl1, frl2, frhl, frl4 = detect_fr_components(light, l1, l2, l3, frl)

        heavy = heavy or '-'
        light = light or '-'

        cdr_info[pdb]['HEAVY'] = heavy
        cdr_info[pdb]['LIGHT'] = light



def create_antibody_db():
    cdr_info, cdr_legend = read_formated_text_file('info/cdr_info')


    frh_info, frh_legend = read_formated_text_file('info/frh_info')
    frl_info, frl_legend = read_formated_text_file('info/frl_info')

    # altering data representation
    for p in cdr_info.keys():
        if p in frh_info: cdr_info[p]['FRH'] = frh_info[p]['FRH']
        if p in frl_info: cdr_info[p]['FRL'] = frl_info[p]['FRL']

        cdr_info[p]['StructSource'] = cdr_info[p]['StructSouce']

        new_key = p.partition('_chothia.pdb')[0][3:]
        cdr_info[new_key] = cdr_info[p]
        del cdr_info[p]

        cdr_info[new_key]['pdb'] = new_key



        for k in cdr_info[new_key]: cdr_info[new_key][k] = '-' if cdr_info[new_key][k]=='none' else cdr_info[new_key][k]

    detect_components(cdr_info)


    for k,nk in dict(H1='h1', H2='h2', H3='h3', L1='l1', L2='l2', L3='l3', FRH='frh', FRL='frl', HEAVY='heavy', LIGHT='light').items():
        for p in cdr_info:
            if k in cdr_info[p]:
                cdr_info[p][nk] = cdr_info[p][k]
                del cdr_info[p][k]

    write_formated_text_file('info/antibody.info', cdr_info, 'pdb resolution BioType date LightType StructSource h1 h2 h3 l1 l2 l3 frh frl heavy light'.split())



def read_fasta_file(file_name):
    """ return array of sequences from FASTA formatted file

    file_name -- relative or absolute path to FASTA file
    """
    seqArray=[]
    seqArrayNames=[]
    seq=""
    for l in file(file_name):
        if l.startswith('>'):
            if not "" == seq:
	        seqArray.append(seq)
		seq=""
            seqArrayNames.append(l.rstrip())
        else:
            seq=seq+l.rstrip()
    seqArray.append(seq)
    return seqArray


# Generating data set for Rosetta CDR unit tests
def create_cdr_test_data():
    tests = [d for d in os.listdir('test') if d != '.svn']

    data = {}

    for t in tests:
        #print 'CDR test target: {}...'.format(t)

        r = json.load( file('test/%s/%s.json' % (t, t) ) )
        for k in r: r[k.lower()] = r.pop(k)

        r['heavy_chain_sequence'] = read_fasta_file('test/%s/query_h.fasta' % t)[0]
        r['light_chain_sequence'] = read_fasta_file('test/%s/query_l.fasta' % t)[0]

        for chain in [ 'L', 'H' ]:
            sequence = []
            numbering = []
            for line in file( 'test/{}/numbering_{}.txt'.format(t, chain) ):
                aa, num = line.split()
                sequence.append(aa)
                numbering.append(num)

            key = dict(H='heavy', L='light')[chain]

            r['trimmed_' + key + '_sequence'] = ''.join(sequence)
            r['trimmed_' + key + '_numbering'] = ' '.join(numbering)


        data[ t.lower() ] = r

    with file('info/cdr-test-data.json', 'w') as f: json.dump(data, f, sort_keys=True, indent=2)



def main(args):
    ''' Script to create antibody.db file from old info/* files and cdr-test-data.json file for unit tests
    '''
    create_antibody_db()

    create_cdr_test_data()


if __name__ == "__main__": main(sys.argv)
