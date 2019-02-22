#!/bin/sh
antibody.macosclangrelease -fasta $1.truncated.fasta -antibody:grafting_database ~/Rosetta/tools/antibody-update -no_relax -antibody:exclude_pdb $1
