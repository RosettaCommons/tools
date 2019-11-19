# This file is meant to list all the #includes that should be left intact.

import re


class DontRemoveInclude:
    def __init__(self):

        self.files_to_leave_intact = set()
        self.regexes_to_leave_intact = []

        self.includes_to_leave_intact = set()
        self.includes_to_leave_be = set()
        self.initialize_keep_intact_lists()

        self.surrogates = {}
        self.initialize_surrogates_for_files()


    # returns True if no headers should be removed
    # from a particular file
    def leave_file_intact(self, fname):
        if fname in self.files_to_leave_intact:
            return True
        if fname in self.surrogates:
            return False
        for regex in self.regexes_to_leave_intact:
            if regex.search(fname):
                return True
        return False

    def attempt_include_removal_for_file(self, fname):
        return not self.leave_file_intact(fname)

    # returns True if included_file_name shold not be removed from source_file
    def leave_include_in_file(self, included_file_name, source_file):
        if (included_file_name, source_file) in self.includes_to_leave_intact:
            return True
        elif included_file_name in self.includes_to_leave_be:
            return True
        return False

    def initialize_keep_intact_lists(self):
        self.initialize_files_to_leave_intact()
        self.initialize_regexes_to_leave_intact()
        self.initialize_includes_to_leave_intact()
        self.initialize_includes_to_leave_be()

    def initialize_files_to_leave_intact(self):
        intact = set()
        intact.add("protocols/forge/methods/pose_mod.hh")
        intact.add("core/id/DOF_ID_Map.Pose.hh")
        intact.add("protocols/viewer/viewers.hh")
        intact.add("protocols/viewer/viewers.cc")
        intact.add("protocols/jd2/JobDistributorFactory.cc")
        intact.add("protocols/jd2/MPIWorkPoolJobDistributor.hh")
        intact.add("protocols/jd2/MPIWorkPoolJobDistributor.cc")
        intact.add("protocols/jd2/MPIWorkPartitionJobDistributor.hh")
        intact.add("protocols/jd2/MPIWorkPartitionJobDistributor.cc")
        intact.add("protocols/jd2/MPIFileBuffer.cc")
        intact.add("protocols/jd2/archive/MPIArchiveJobDistributor.hh")
        intact.add("protocols/jd2/archive/MPIArchiveJobDistributor.cc")
        intact.add("core/util/prof.hh")
        intact.add("core/scoring/etable/etrie/CountPairDataGeneric.hh")
        intact.add("devel/InvKinLigLoopDesign/ints.hh")
        intact.add("utility/vectorL.hh")
        intact.add("utility/vector1.hh")
        intact.add("utility/vector0.hh")
        self.files_to_leave_intact = intact

    def initialize_regexes_to_leave_intact(self):
        self.regexes_to_leave_intact.append(re.compile("\.tmpl\.hh$"))
        self.regexes_to_leave_intact.append(re.compile("\.impl\.hh$"))
        self.regexes_to_leave_intact.append(re.compile("\.ipp$"))

    def initialize_includes_to_leave_intact(self):
        pairs = [
            ("cassert", "utility/backtrace.hh"),
            ("utility/backtrace.hh", "utility/assert.hh"), # for debug assert
            ("utility/exit.hh", "utility/to_string.hh"),
            ("core/types.hh", "core/pack/rotamer_set/RotamerSetsBase.hh"),
            ("core/types.hh", "core/grid/CartGrid.hh"),
            ("core/types.hh", "core/pack/rotamer_trials.hh"),
            ("core/types.hh", "core/pack_basic/RotamerSetsBase.hh"),
            ("core/types.hh", "protocols/dna/DNAParameters.hh"),
            ("core/types.hh", "protocols/pockets/PocketGrid.hh"),
            ("core/types.hh", "protocols/scoring/methods/pcs2/PcsDataLanthanide.hh"),
            ("utility/exit.hh", "protocols/moves/DataMap.hh"),
            ("string", "protocols/moves/DataMap.hh"),
            ("utility/basic_sys_util.hh", "core/init.cc"),
            ("string", "protocols/ligand_docking/ligand_options/ProtocolOption.hh"),
            ("utility/exit.hh", "core/kinematics/Edge.hh"),
            ("cassert", "core/kinematics/Edge.hh"),
            ("string", "protocols/frags/VallData.hh"),
            ("string", "protocols/filters/Filter.fwd.hh"),
            (
                "protocols/jd2/JobDistributorFactory.fwd.hh",
                "protocols/jd2/archive/MPIArchiveJobDistributor.hh",
            ),
            ("cstdlib", "protocols/jumping/PairingsList.hh"),
            ("cstring", "core/init.cc"),
            ("string", "core/scoring/geometric_solvation/DatabaseOccSolEne.hh"),
            ("cstring", "core/io/serialization/serialize_pose.hh"),
            ("cstring", "core/scoring/packstat/io.cc"),
            ("utility", "devel/InvKinLigLoopDesign/Ints.hh"),
            ("string", "devel/InvKinLigLoopDesign/Ints.hh"),
            ("string", "core/scoring/AtomVDW.hh"),
            ("string", "core/scoring/constraints/FuncFactory.hh"),
            ("cassert", "devel/InvKinLigLoopDesign/std_extra.hh"),
            ("cassert", "devel/simple_options/option.cc"),
            ("cstring", "devel/simple_options/option.cc"),
            ("utility/keys/Key2Tuple.hh", "core/chemical/ResidueType.hh"),
            ("utility/keys/Key3Tuple.hh", "core/chemical/ResidueType.hh"),
            ("utility/keys/Key4Tuple.hh", "core/chemical/ResidueType.hh"),
            ("core/chemical/Adduct.hh", "core/chemical/ResidueType.hh"),
            ("core/chemical/Adduct.hh", "core/chemical/ResidueTypeBase.hh"),
            ("core/id/DOF_ID.hh","core/id/DOF_ID_Map.hh"),
            (
                "core/kinematics/ShortestPathInFoldTree.hh",
                "protocols/abinitio/MaxSeqSepConstraintSet.hh",
            ),
            ("core/kinematics/Jump.hh", "core/io/silent/ProteinSilentStruct.hh"),
            ("core/conformation/Residue.hh", "core/pose/util.tmpl.hh"),
            ("core/conformation/Conformation.hh", "core/pose/util.tmpl.hh"),
            ("core/pack/interaction_graph/LinearMemoryInteractionGraph.hh", "core/pack/interaction_graph/SurfaceInteractionGraph.hh"),
            ("core/pack/interaction_graph/LinearMemoryInteractionGraph.hh", "core/pack/interaction_graph/HPatchInteractionGraph.hh"),
            ("core/pack/task/PackerTask.hh", "core/pack/interaction_graph/HPatchInteractionGraph.hh"),
            ("core/pack/rotamer_set/RotamerSet.hh", "core/pack/interaction_graph/HPatchInteractionGraph.hh"),
            ("core/pack/rotamer_set/RotamerSets.hh", "core/pack/interaction_graph/HPatchInteractionGraph.hh"),
            ("core/scoring/TenANeighborGraph.hh", "core/pack/interaction_graph/HPatchInteractionGraph.hh"),
            ("core/pack/interaction_graph/LinearMemoryInteractionGraph.hh", "core/pack/interaction_graph/NPDHBondInteractionGraph.hh"),
            ("core/conformation/Residue.hh", "core/pack/interaction_graph/NPDHBondInteractionGraph.hh"),
            ("core/pack/task/PackerTask.hh", "core/pack/interaction_graph/NPDHBondInteractionGraph.hh"),
            ("core/pack/rotamer_set/RotamerSet.hh", "core/pack/interaction_graph/NPDHBondInteractionGraph.hh"),
            ("core/pack/rotamer_set/RotamerSets.hh", "core/pack/interaction_graph/NPDHBondInteractionGraph.hh"),
            ("core/pose/Pose.hh", "core/pack/interaction_graph/NPDHBondInteractionGraph.hh"),
            ("core/scoring/ScoreFunction.hh", "core/pack/interaction_graph/NPDHBondInteractionGraph.hh"),
            ("core/scoring/Energies.hh", "core/pack/interaction_graph/NPDHBondInteractionGraph.hh"),
            ("core/scoring/hbonds/NPDHBondSet.hh", "core/pack/interaction_graph/NPDHBondInteractionGraph.hh"),
            ("core/scoring/hbonds/HBondOptions.hh", "core/pack/interaction_graph/NPDHBondInteractionGraph.hh"),
            ("utility/graph/Graph.hh", "core/pack/interaction_graph/NPDHBondInteractionGraph.hh"),
            ("core/scoring/EnergyMap.hh", "core/scoring/etable/BaseEtableEnergy.hh"),
            ("core/conformation/Residue.hh", "core/scoring/etable/BaseEtableEnergy.hh"),
            ("core/conformation/Atom.hh","core/scoring/etable/atom_pair_energy_inline.hh"),
            ("core/conformation/Residue.hh","core/scoring/etable/atom_pair_energy_inline.hh"),
            ("core/chemical/AtomType.hh","core/scoring/etable/atom_pair_energy_inline.hh"),
            ("numeric/xyzVector.hh","core/scoring/etable/atom_pair_energy_inline.hh"),
            ("core/conformation/Atom.hh","core/scoring/etable/count_pair/CountPair1B.hh"),
            ("core/conformation/Residue.hh","core/scoring/etable/count_pair/CountPair1B.hh"),
            ("core/conformation/Atom.hh", "core/scoring/etable/count_pair/CountPair2B.hh"),
            ("core/conformation/Residue.hh", "core/scoring/etable/count_pair/CountPair2B.hh"),
            ("core/conformation/Residue.hh", "core/scoring/etable/count_pair/CountPairIntraRes.hh"),
            ("core/pose/util.tmpl.hh", "core/scoring/packing/surf_vol.cc"),
            ("core/scoring/trie/TrieCountPairBase.hh", "core/scoring/trie/RotamerTrie.hh"),
            ("core/scoring/hbonds/hbtrie/HBCPData.hh", "core/scoring/trie/RotamerTrie.hh"),
            ("core/scoring/lkball/lkbtrie/LKBAtom.hh", "core/scoring/trie/RotamerTrie.hh"),
            ("core/scoring/trie/CPDataCorrespondence.hh", "core/scoring/trie/trie.functions.hh"),
            ("core/scoring/trie/RotamerTrieBase.hh", "core/scoring/trie/trie.functions.hh"),
            ("core/conformation/Residue.hh", "core/scoring/trie/trie.functions.hh"),
            ("core/conformation/RotamerSetBase.hh", "core/scoring/trie/trie.functions.hh"),
            ("protocols/antibody/grafting/util.hh", "protocols/antibody/grafting/cdr_detection.hh"),
            ("protocols/antibody/grafting/util.hh", "protocols/antibody/grafting/json_based_cdr_detection.hh"),
            ("utility/json_utilities.hh", "protocols/antibody/grafting/json_based_cdr_detection.hh"),
            ("protocols/antibody/grafting/util.hh", "protocols/antibody/grafting/scs_multi_template.hh"),
            ("utility/json_utilities.hh", "protocols/ddg/CartesianddG.hh"),
            ("core/io/silent/SilentStruct.hh", "protocols/evaluation/PoseEvaluator.hh"),
            ("core/conformation/Residue.hh", "protocols/features/RotamerFeatures.hh"),
            ("core/kinematics/FoldTree.hh", "protocols/forge/methods/fold_tree_functions.hh"),
            ("core/conformation/Residue.hh", "protocols/indexed_structure_store/FragmentLookup.hh"),
            ("core/pose/Pose.hh", "protocols/indexed_structure_store/FragmentLookup.hh"),
            ("core/pack/task/TaskFactory.hh", "protocols/loop_modeling/refiners/packing_helper.hh"),
            ("core/pack/task/PackerTask.hh", "protocols/loop_modeling/refiners/packing_helper.hh"),
            ("core/pose/Pose.hh", "protocols/loop_modeling/refiners/packing_helper.hh"),
            ("utility/tag/Tag.hh", "protocols/loop_modeling/utilities/rosetta_scripts.hh"),
            ("protocols/noesy_assign/CrossPeak.hh", "protocols/noesy_assign/NoesyModule.impl.hh"),
            ("protocols/noesy_assign/CrossPeakList.hh", "protocols/noesy_assign/NoesyModule.impl.hh"),
            ("protocols/fldsgn/topology/SS_Info2.hh", "protocols/sic_dock/scores/MotifHashRigidScore.cc"),
            ("protocols/antibody/grafting/util.hh", "protocols/tcr/TCRloopRefine.hh"),
            ("protocols/antibody/grafting/util.hh", "protocols/tcr/TCRseqInfo.hh"),
            ("protocols/antibody/grafting/util.hh", "protocols/tcr/util.hh"),
            ("core/io/silent/SilentStruct.hh", "protocols/toolbox/Cluster.impl.hh"),
            ("core/conformation/Residue.hh", "protocols/toolbox/DecoySetEvaluation.impl.hh"),
            ("core/pose/Pose.hh", "protocols/toolbox/PyReturnValuePolicyTest.hh"),
            ("core/scoring/ScoreFunction.hh", "protocols/toolbox/PyReturnValuePolicyTest.hh"),
            ("core/conformation/parametric/ParametersSet.hh", "protocols/viewer/viewers.cc"), # I'm not sure about this one
            (
                "protocols/jd2/MPIWorkPoolJobDistributor.hh",
                "apps/public/rosetta_scripts/rosetta_scripts.cc",
            ),
            (
                "protocols/jd2/MPIFileBufJobDistributor.hh",
                "apps/public/rosetta_scripts/rosetta_scripts.cc",
            ),
            (
                "protocols/jd2/BOINCJobDistributor.hh",
                "apps/public/rosetta_scripts/rosetta_scripts.cc",
            ),
            (
                "protocols/jd2/BOINCJobDistributor.hh",
                "protocols/jd2/JobDistributorFactory.cc",
            ),
            (
                "protocols/jd2/MPIWorkPartitionJobDistributor.hh",
                "protocols/jd2/JobDistributorFactory.cc",
            ),
            (
                "protocols/jd2/MPIWorkPoolJobDistributor.hh",
                "protocols/jd2/JobDistributorFactory.cc",
            ),
            (
                "protocols/jd2/archive/MPIArchiveJobDistributor.hh",
                "protocols/jd2/JobDistributorFactory.cc",
            ),
            ("ObjexxFCL/FArray1D.hh", "core/kinematics/DomainMap.hh"),
            (
                "numeric/deriv/dihedral_deriv.hh",
                "core/scoring/constraints/DihedralConstraint.cc",
            ),
            ("iostream", "devel/domain_assembly/DomainAssemblyReader.cc"),
            (
                "core/scoring/mm/MMBondAngleResidueTypeParam.hh",
                "protocols/branch_angle/BranchAngleOptimizer.cc",
            ),
            ("iostream", "protocols/genetic_algorithm/Entity.cc"),
            ("iostream", "protocols/pack_daemon/MultistateFitnessFunction.cc"),
            ("iostream", "protocols/smanager/smanager.cc"),
            ("iostream", "protocols/wum/WorkUnitManager.cc"),
            ("iostream", "protocols/optimize_weights/Arithmetic.cc"),
            ("iostream", "core/scoring/dssp/StrandPairing.cc"),
            ("core/id/DOF_ID_Range.hh", "protocols/moves/ThermodynamicMover.hh"),
            ("core/id/TorsionID_Range.hh", "protocols/moves/ThermodynamicMover.hh"),
            ("iostream", "core/pack/interaction_graph/DensePDInteractionGraph.cc"),
            (
                "core/optimization/MinimizerOptions.hh",
                "apps/benchmark/Minimizer.bench.hh",
            ),
            ("iostream", "apps/pilot/will/test_string.cc"),
            (
                "core/pack/dunbrack/SingleResidueDunbrackLibrary.tmpl.hh",
                "core/pack/dunbrack/RotamerLibrary.cc",
            ),
            ("core/fragment/FrameIteratorWorker_.hh", "core/fragment/FrameIterator.hh"),
            ("core/fragment/Frame.hh", "core/fragment/FragCache.hh"),
            ("core/fragment/FragData.hh", "core/fragment/picking_old/vall/util.hh"),
            ("protocols/docking/DockingProtocol.hh", "apps/benchmark/performance/Docking.bench.hh"),
            ("core/chemical/Atom.hh", "core/chemical/ResidueGraphTypes.hh"),
            ("boost/graph/filtered_graph.hpp", "core/chemical/ResidueGraphTypes.hh"),
            ("core/chemical/AtomTypeSet.fwd.hh", "core/chemical/ResidueGraphTypes.hh"),
            ("boost/graph/adjacency_list.hpp", "core/chemical/ResidueGraphTypes.hh"),
            ("utility", "core/chemical/ResidueGraphTypes.hh"),
            ("iostream", "core/io/nmr/ParaIon.cc"), # but why would it remove it at all??
            ("iostream", "core/scoring/func/SumFunc.cc"), # maybe unnecessary?
            ("iostream", "core/scoring/etable/count_pair/CountPairFunction.cc"),
            ("ostream", "devel/denovo_protein_design/SSClass.hh"),
            ("utility/stream_util.hh", "../test/core/io/mmtf_IO.cxxtest.hh"),
            ("utility/stream_util.hh", "../test/protocols/fldsgn/SheetConstraintGenerator.cxxtest.hh"),
            ("utility/stream_util.hh", "../test/protocols/denovo_design/components/DivideAndConquerorTests.cxxtest.hh"),
            ("numeric/VoxelGrid.impl.hh", "../test/core/scoring/nmr/NMRDummySpinlabelEnsemble.cxxtest.hh"),
        ]

        for pair in pairs:
            self.includes_to_leave_intact.add(pair)

    def initialize_includes_to_leave_be(self):
        self.includes_to_leave_be.add("utility/vector1.hh")
        self.includes_to_leave_be.add("utility/vector0.hh")

    def initialize_surrogates_for_files(self):
        surrogates = [
            ("protocols/loops/Loops.tmpl.hh", "protocols/simple_task_operations/RestrictToLoops.cc"),
            ("protocols/canonical_sampling/BiasEnergy.tmpl.hh", "protocols/canonical_sampling/BiasEnergy.cc"),
            ("protocols/moves/MonteCarlo.tmpl.hh", "protocols/viewer/viewers.cc"),
            ("numeric/interpolation/InterpolatedPotential.tmpl.hh", "core/scoring/rna/data/RNA_DMS_Potential.cc"),
            ("numeric/interpolation/spline/PolycubicSpline.tmpl.hh", "protocols/features/RotamerFeatures.cc"),
            ("numeric/angle.functions.hh", "core/chemical/ring/RingConformerSet.cc"),
            ("core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh", "core/pack/dunbrack/RotamerLibrary.cc"),
            ("core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.tmpl.hh", "core/pack/dunbrack/RotamerLibrary.cc"),
            ("core/pack/dunbrack/RotamericSingleResidueDunbrackLibraryParser.tmpl.hh", "core/pack/dunbrack/SingleResidueDunbrackLibrary.cc"),
            ("core/pose/util.tmpl.hh", "protocols/helical_bundle/PerturbBundle.cc"),
            ("core/pose/util.tmpl.hh", "protocols/forge/remodel/RemodelMover.cc"),
            ("core/io/silent/ProteinSilentStruct.tmpl.hh", "protocols/loophash/WorkUnit_LoopHash.cc"),
            #("core/io/silent/ProteinSilentStruct.tmpl.hh", "protocols/mpi_refinement/WorkUnit_Relax.cc"),
            #("core/io/silent/ProteinSilentStruct.tmpl.hh", "apps/pilot/mike/loophash.cc"),
            ("core/scoring/rms_util.tmpl.hh", "protocols/protein_interface_design/filters/RmsdSimpleFilter.cc"),
            ("core/scoring/etable/BaseEtableEnergy.tmpl.hh", "core/scoring/etable/EtableEnergy.cc"),
            ("core/scoring/NeighborList.tmpl.hh", "core/energy_methods/StackElecEnergy.cc"),
            ("utility/FixedSizeLexicographicalIterator.tmpl.hh", "protocols/match/MatchSet.cc"),
            ("protocols/ligand_docking/ligand_dock_impl.hh", "apps/public/ligand_docking/ligand_dock.cc"),
            ("protocols/noesy_assign/CrossPeakList.impl.hh", "protocols/noesy_assign/NoesyModule.cc"),
            ("protocols/noesy_assign/NoesyModule.impl.hh", "protocols/noesy_assign/NoesyModule.cc"),
            ("protocols/toolbox/Cluster.impl.hh", "apps/public/analysis/fast_clustering.cc"),
            ("protocols/toolbox/DecoySetEvaluation.impl.hh", "protocols/toolbox/DecoySetEvaluation.cc"),
            ("numeric/alignment/rmsd_calc.impl.hh", "numeric/alignment/rmsd_calc.cc"),
            ("numeric/VoxelGrid.impl.hh", "protocols/nmr/pre/PREEnergy.cc"),
        ]
        
        for surrpair in surrogates:
            if surrpair[0] not in self.surrogates:
                self.surrogates[surrpair[0]] = []
            self.surrogates[surrpair[0]].append(surrpair[1])
