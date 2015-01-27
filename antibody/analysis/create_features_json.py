#!/usr/bin/python
# Author: Jared Adolf-Bryfogle

#This script will create either cluster features or antibody features json for use in Features R script.
#I am just really sick of doing this by hand.
#Example Cmd-line:  python create_features_json.py --database databases/baseline_comparison.txt --scripts cluster

from optparse import OptionParser, IndentedHelpFormatter
import os
import sys
import json
from collections import defaultdict




def setup_baseline_scripts_and_formats(json_dict, type):
    current_path = os.path.split(os.path.abspath(__file__))[0]

    FILE = open(current_path+"/"+type+"_features.json")
    data = json.load(FILE)
    FILE.close()
    append_scripts_formats_to_json_dict(data, json_dict)

def append_scripts_formats_to_json_dict(data, json_dict):
    json_dict["sample_source_comparisons"][0]["output_formats"].extend(data["sample_source_comparisons"][0]["output_formats"])
    json_dict["sample_source_comparisons"][0]["analysis_scripts"].extend(data["sample_source_comparisons"][0]["analysis_scripts"])

def initialize_json_dict(out_dir):

    formats = defaultdict()
    formats["output_formats"] = []
    formats["analysis_scripts"] = []
    formats["sample_sources"] = []

    json_dict = defaultdict()
    json_dict["sample_source_comparisons"] = []
    json_dict["sample_source_comparisons"].append(formats)
    json_dict["output_dir"] = out_dir

    return json_dict


def add_sample_source_comparisons(current_json_dict, json_file_path = None):
    if json_file_path:
        if not os.path.exists(json_file_path):
            sys.exit("Passed Json file path does not exist!!")
        FILE = open(json_file_path, 'r')
        json_dict = json.load(FILE)
        FILE.close()
        current_json_dict["sample_source_comparisons"].extend(json_dict["sample_source_comparisons"])

def add_sample_source(json_dict, sample_source_dict):
    json_dict["sample_source_comparisons"][0]["sample_sources"].append(sample_source_dict)


def write_json_for_single_recovery_experiment(db_path_exp, db_path_natives, exp_id, out_path):
    """
    Create a JSON file for recovery of a single experiment.
    """
    if not os.path.exists(db_path_exp):
        sys.exit("SQLITE3 database path does not exist!!")

    json_dict = initialize_json_dict()
    json_dict["output_dir"] = out_path
    setup_baseline_scripts_and_formats(json_dict, "cluster")

    #Ref HAS to go first here!
    add_sample_source(json_dict, _create_json_info_dic(db_path_natives, "ref", True))
    add_sample_source(json_dict, _create_json_info_dic(db_path_exp, exp_id, False))


    if not os.path.exists(out_path):
        os.mkdir(out_path)
    if not os.path.exists(out_path+"/jsons"):
        os.mkdir(out_path+"/jsons")
    OUTFILE = open(out_path+"/jsons/"+"cluster_features."+exp_id+".json", 'w')
    json.dump(json_dict, OUTFILE, indent=1)
    OUTFILE.close()

    print "Complete"


########################################################################################################################
def _create_json_info_dic(db_path, id, ref=False):
    """
    Create a simple dict that holds info to create the json.
    """
    databases = defaultdict()
    databases["database_path"] = db_path
    databases["id"] = id
    databases["reference"] = ref

    return databases



########################################################################################################################
####### Class based implementation
########################################################################################################################

class JsonCreator:
    """
    Basic implementation of a simple JsonCreator to create Jsons.  Could be expanded to not load jsons with pre-set scripts.
    A nicer implementation would be a GUI for running the FeaturesReporter scripts.
    """
    def __init__(self, out_path, script_type):
        self.script_types = ["antibody", "interface", "cluster", "antibody_min_hbond_analysis", "antibody_no_hbond_analysis"]
        if not script_type in self.script_types:
            sys.exit(script_type +" unrecognized.  Available JSON script types are: "+repr(self.script_types))

        self.json_dict = initialize_json_dict(out_path)
        setup_baseline_scripts_and_formats(self.json_dict, script_type)

    def add_sample_source_info(self, db_path, id, ref = False):
        info = _create_json_info_dic(db_path, id, ref)
        add_sample_source(self.json_dict, info)

    def add_current_sample_sources_to_json_dict(self, json_file_path):
        """
        Combine a JSON with data held in this class
        """
        if os.path.exists(json_file_path):
            add_sample_source_comparisons(self.json_dict, json_file_path)

    def save_json(self, out_path = "local_json.txt"):
        OUTFILE = open(out_path, 'w')


        json.dump(self.json_dict, OUTFILE, indent=1)
        OUTFILE.close()
        self.json_path = out_path

        print "Json written to: "+out_path
        print "Json path set to JsonCreator"

    def run_json(self):
        run_features_json(self.json_path)


def run_features_json(json_path):
    """
    Convenience function
    Run compare_sample_sources with json path.  Must have compare_sample_sources.R in PATH
    """
    r_cmd = "compare_sample_sources.R --config "+json_path
    print "Running: "+r_cmd
    os.system(r_cmd)

if __name__ == "__main__":


    parser = OptionParser()

    parser.add_option("--databases", "-d",
        help = "Path to database file.  2 - 3 columns.  "
               "db_name, short_name, ref keyword if the reference database. One entry per line. Comments ok. "
    )

    parser.add_option("--script", "-s",
        help = "Script type.  Options are cluster, antibody, interface.  Will setup the appropriate output formats and R scripts"
    )

    parser.add_option("--db_path", "-p",
        help = "Path to databases.  Default is pwd/databases",
        default=os.getcwd()+"/databases")

    parser.add_option("--out_path", "-o",
        help = "Where to put the result of the analysis scripts.  Currently unsupported by the features framework.",
        default=os.getcwd())

    parser.add_option("--out_name", "-n",
        help = "Output file name of json file",
        default = "local_json.json" )

    parser.add_option("--add_comparison_to_this_json", "-a",
        help = "Add all this data as another sample source comparison to passed in json file path. Untested")

    (options, args) = parser.parse_args(sys.argv[1:])

    if not os.path.exists(options.databases):
        sys.exit("Database file does not exist!!")

    if not os.path.exists(options.db_path):
        sys.exit("SQLITE3 database path does not exist!!")

    script_types = ["cluster", "antibody", "interface"]

    if not options.scripts in script_types:
        sys.exit("Unrecognized script type. Options are: "+repr(script_types))


    json_creator = JsonCreator(options.out_path, options.script)

    INFILE = open(options.databases, "r")
    for line in INFILE:

        line = line.strip()
        if not line: continue
        if line[0] == "#" or line[0] == " ": continue

        lineSP = line.split()
        ref = False
        if len(lineSP) == 3 and lineSP[2] == "ref": ref = True

        db_path = options.db_path+"/"+lineSP[0]; id= lineSP[1]

        json_creator.add_sample_source_info(db_path, id, ref)
        #print repr(json_dict)

    INFILE.close()

    json_creator.add_current_sample_sources_to_json_dict(options.add_comparison_to_this_json)
    json_creator.save_json(options.out_path)

    print "Complete"
