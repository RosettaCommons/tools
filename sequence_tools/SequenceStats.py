#Author: Jared Adolf-Bryfogle

import sys
from restype_definitions import definitions

class SequenceStats:
    """
    Class for getting data from an array of strings of sequences (one letter code) of equal length.
    """

    def __init__(self, sequence_list,):
        self.sequence_list = sequence_list
        if len(self.sequence_list)==0:
            print "No sequences found."
            return
        self.stats = []; #Vector of a map of each amino acid
        self.freq = []; #Same as stats

        self.aas = definitions().get_all_one_letter_codes()
        self.init_data_map()
        self.compute_stats()

    def set_sequences(self, sequence_list):
        """
        Set a sequence list
        """

        self.sequence_list = sequence_list
        self.compute_stats()

    def get_probability(self, position, aa):
        '''
        Get probability of the current position (starting from 0) and aa
        '''

        try:
            if not self.stats[position].has_key(aa):
                return 0
        except IndexError:
            sys.exit("Position not found in stats.  This is bad. We have a problem")

        return self.stats[position][aa]

    def get_consensus_sequence(self):
        consensus = ""
        for position in range(0, len(self.stats)):
            pos = position
            aa_letter = ""
            for aa in self.stats[pos]:
                prob = self.stats[pos][aa]
                if prob > .9:
                    aa_letter = aa.upper()
                    break
                elif prob > .2:
                    aa_letter = aa.lower()
                    break
                else:
                    aa_letter = '-'
                    continue

            consensus = consensus+aa_letter

        return consensus

    def get_frequency(self, position, aa):
        return self.freq[position][aa]

    ##Private Methods##
    def compute_stats(self):
        """
        Compute frequency and probability (0-1) for each position for each amino acid
        """
        if len(self.sequence_list)==0:return

        #print "Total Unique CDR Sequences: "+repr(len(self.sequence_list))
        for position in range(0, len(self.sequence_list[0])):

            totals = self.return_initialized_total_map()
            for seq in self.sequence_list:

                totals[seq[position]]+=1
            for aa in totals:

                self.stats[position][aa] = totals[aa]/float(len(self.sequence_list))
                self.freq[position][aa] = totals[aa]

    def init_data_map(self):
        """
        Sets all probabilities 0 and appends each map to the stats vector
        """

        for i in range(1, len(self.sequence_list[0])+1):
            prob = dict()
            freq = dict()
            for aa in self.aas:
                prob[aa] = 0
                freq[aa] = 0
            self.stats.append(prob)
            self.freq.append(freq)

    def return_initialized_total_map(self):
        totals = dict()
        for aa in self.aas:
            totals[aa] = 0

        return totals

