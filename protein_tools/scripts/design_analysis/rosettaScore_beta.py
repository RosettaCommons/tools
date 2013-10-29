import sys
import warnings
from Bio.PDB import PDBExceptions
from Bio.PDB import PDBParser as parser
import amino_acids
warnings.simplefilter('ignore', PDBExceptions.PDBConstructionWarning)

if sys.version_info < (2, 7):
    raise Exception("You must use python2.7 to call rosettaScore_beta class")


exception = ["CYD", "TYS", "CYS"]


def _get_score_table(path):
    score_table = []
    table = False
    with open(path) as openfile:
        for count, line in enumerate(openfile, start=1):
            try:
                if line.split()[0] == "#BEGIN_POSE_ENERGIES_TABLE":
                    table = True
                    continue
                if table and line.split()[0] != "#END_POSE_ENERGIES_TABLE":
                    score_table.append(line)
            except IndexError:
                print "an index error occured on line {0}, probably can't split this line".format(count)
    if score_table:
        return score_table
    else:
        raise Exception("The pdb {0} you are trying to pass does not contain a score file".format(path))


def get_table(path):
    raw_table = []
    infile = open(path, 'r')
    table = False
    for line in infile:
        line_split = line.split()
        try:
            if line_split[0] == "END_POSE_ENERGIES_TABLE":
                break
            if line_split[0] == "#BEGIN_POSE_ENERGIES_TABLE":
                table = True
                raw_table.append(line)
            elif table:
                raw_table.append(line)
        except IndexError:
            print "pdb file probably contains a blank line"
            continue
    infile.close()
    return raw_table


class ScoreRecord:
    def __init__(self, name, scores, chainid, resid, native):
        self.name = name
        self.scores = scores
        self.resid = resid
        self.chain = chainid
        self.native = native


class ScorePoseRecord:
    def __init__(self, name, poseid, scores, native):
        self.name = name
        self.scores = scores
        self.poseid = poseid
        self.native = native


class ScoreTable:
    def __init__(self, path, already_open=False):
        # get the energy table by a function that will only return the lines in the rosetta score table
        energy_table = _get_score_table(path)
        header = energy_table[0].split()[1:]
        weights = energy_table[1].split()[1:]
        pose = energy_table[2].split()[1:]
        residue_energies = energy_table[3:]
        # map pose_id_to_res_id
        self.pose_id_mapped_to_res_id = {}
        # we can hold on to two dictionaries, weights and records.
        self.weights = {}
        self.records = {}
        self.pose = {}
        if already_open:
            pp = path.get_residues()
        else:
            pp = parser().get_structure(path.split(".")[-4:-1], path).get_residues()
        for residue, line in zip(pp, residue_energies):
            name_pose_res = line.split("_")[0]
            name_pdb_res = residue.get_resname()
            if name_pose_res != name_pdb_res and name_pdb_res not in exception:
                raise Exception("%s,%spose name and res name do not match" % (name_pose_res, name_pdb_res))

            resid = int(residue.get_id()[1])
            poseid = int(line.split()[0].split("_")[-1])
            chainid = residue.get_parent().get_id()
            scores = line.split()[1:]
            scoredictionary = {}
            for term, score in zip(header, scores):
                scoredictionary[term] = float(score)
            self.records[(resid, chainid)] = ScoreRecord(name_pose_res, scoredictionary, chainid, resid, name_pose_res)
            self.records[poseid] = ScorePoseRecord(name_pose_res, poseid, scoredictionary, name_pdb_res)
            self.pose_id_mapped_to_res_id[(chainid, resid)] = poseid

        for term, weight in zip(header, weights):
            self.weights[term] = weight

        for term, score in zip(header, pose):
            self.pose[term] = float(score)

    def get_score(self, poseid=None, term=None, chain=None, pdbres=None):
        if chain and pdbres:
            try:
                self.score_record = self.records[(int(pdbres), chain)]
            except KeyError:
                print "No residue with the chain %s,%i exist in this pdb file - your results may be empty" % (chain, pdbres)
        elif poseid:
            try:
                self.score_record = self.records[poseid]
            except KeyError:
                print "This %i is not a viable pose number" % poseid
        else:
            raise Exception("You must put in a pose id, or a residue id and a chain id so the score can be determined")
        if term:
            return (self.score_record.scores[term], self.score_record.name)
        else:
            return (self.score_record.scores, self.score_record.name)

    def get_all_score_terms(self, chain, pdbres):
        score = self.get_score(chain=chain, pdbres=pdbres)[0]
        return score

    def get_all_score_terms_by_pose_id(self, pose_id):
        score = self.get_score(poseid=pose_id)[0]
        return score

    def get_pose_all_scores(self):
        return self.pose

    def get_weight(self, term):
        return self.weights[term]

    def get_pose_score(self, term):
        return self.pose[term]

    def get_pose_id_from_chain_id(self, chain, pdbres):
        return self.pose_id_mapped_to_res_id[(chain, int(pdbres))]
if __name__ == '__main__':
    if sys.version_info < (2, 7):
        raise Exception("You must use python2.7 to call rosettaScore_beta class")

    file = "/blue/crowelab/home/willisjr/healthy_donor/illumina/pg_threading_with_complex/per_ddg_models/21252:127887_input_complex_0001_0001.pdb"
    print ScoreTable(file).get_score(chain="G", pdbres=30, term="fa_atr")
    print ScoreTable(file).get_score(poseid=471, term="fa_atr")
    print ScoreTable(file).get_weight("fa_atr")
    print ScoreTable(file).get_pose_score("total")
