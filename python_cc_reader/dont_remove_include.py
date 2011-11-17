# This file is meant to list all the #includes that should be left intact.

import re

class DontRemoveInclude :
   def __init__( self ) :
      self.files_to_leave_intact = set()
      self.regexes_to_leave_intact = []
      self.includes_to_leave_intact = set()
      self.includes_to_leave_be = set()
      self.initialize_keep_intact_lists()

   # returns True if no headers should be removed
   # from a particular file
   def leave_file_intact( self, fname ) :
      if fname in self.files_to_leave_intact :
         return True
      for regex in self.regexes_to_leave_intact :
         if regex.search( fname ) :
            return True
      return False

   def attempt_include_removal_for_file( self, fname ) :
     return not self.leave_file_intact( fname )

   # returns True if included_file_name shold not be removed from source_file
   def leave_include_in_file( self, included_file_name, source_file ) :
      if ( included_file_name, source_file ) in self.includes_to_leave_intact :
         return True
      elif included_file_name in self.includes_to_leave_be :
         return True
      return False

   def initialize_keep_intact_lists( self ) :
      self.initialize_files_to_leave_intact()
      self.initialize_regexes_to_leave_intact()
      self.initialize_includes_to_leave_intact()
      self.initialize_includes_to_leave_be()

   def initialize_files_to_leave_intact( self ) :
      intact = set()
      intact.add( "protocols/forge/methods/pose_mod.hh" )
      intact.add( "core/id/DOF_ID_Map.Pose.hh" )
      intact.add( "protocols/viewer/viewers.hh" )
      intact.add( "protocols/viewer/viewers.cc" )
      intact.add( "protocols/jd2/JobDistributorFactory.cc" )
      intact.add( "protocols/jd2/MPIWorkPoolJobDistributor.hh" )
      intact.add( "protocols/jd2/MPIWorkPoolJobDistributor.cc" )
      intact.add( "protocols/jd2/MPIWorkPartitionJobDistributor.hh" )
      intact.add( "protocols/jd2/MPIWorkPartitionJobDistributor.cc" )
      intact.add( "protocols/jd2/MPIFileBuffer.cc" )
      intact.add( "protocols/jd2/archive/MPIArchiveJobDistributor.hh" )
      intact.add( "protocols/jd2/archive/MPIArchiveJobDistributor.cc" )
      intact.add( "core/util/prof.hh" )
      intact.add( "core/scoring/etable/etrie/CountPairDataGeneric.hh" )
      intact.add( "devel/InvKinLigLoopDesign/ints.hh" )
      self.files_to_leave_intact = intact

   def initialize_regexes_to_leave_intact( self ) :
      self.regexes_to_leave_intact.append( re.compile( "\.tmpl\.hh$" ) )
      self.regexes_to_leave_intact.append( re.compile( "\.ipp$" ) )

   def initialize_includes_to_leave_intact( self ) :
      pairs = [ ( "core/types.hh", "core/pack/rotamer_set/RotamerSetsBase.hh" ), \
                ( "utility/exit.hh", "protocols/moves/DataMap.hh" ), \
                ( "string", "protocols/moves/DataMap.hh" ), \
                ( "utility/basic_sys_util.hh", "core/init.cc" ), \
                ( "string", "protocols/ligand_docking/ligand_options/ProtocolOption.hh" ), \
                ( "utility/exit.hh", "core/kinematics/Edge.hh" ), \
                ( "cassert", "core/kinematics/Edge.hh" ), \
                ( "string", "protocols/frags/VallData.hh" ), \
                ( "string", "protocols/filters/Filter.fwd.hh" ), \
                ( "protocols/jd2/JobDistributorFactory.fwd.hh", "protocols/jd2/archive/MPIArchiveJobDistributor.hh" ), \
                ( "cstdlib", "protocols/jumping/PairingsList.hh" ), \
                ( "cstring", "core/init.cc" ), \
                ( "string", "core/scoring/geometric_solvation/DatabaseOccSolEne.hh" ), \
                ( "cstring", "core/io/serialization/serialize_pose.hh" ), \
                ( "cstring", "core/scoring/packstat/io.cc" ), \
                ( "utility", "devel/InvKinLigLoopDesign/Ints.hh" ), \
                ( "string", "devel/InvKinLigLoopDesign/Ints.hh" ), \
                ( "string", "core/scoring/AtomVDW.hh" ), \
                ( "string", "core/scoring/constraints/FuncFactory.hh" ), \
                ( "cassert", "devel/InvKinLigLoopDesign/std_extra.hh" ), \
                ( "cassert", "devel/simple_options/option.cc" ), \
                ( "cstring", "devel/simple_options/option.cc" ), \
                ( "utility/keys/Key2Tuple.hh", "core/chemical/ResidueType.hh" ), \
                ( "core/kinematics/ShortestPathInFoldTree.hh", "protocols/abinitio/MaxSeqSepConstraintSet.hh" ), \
                ( "protocols/jd2/MPIWorkPoolJobDistributor.hh", "apps/public/rosetta_scripts/rosetta_scripts.cc" ), \
                ( "protocols/jd2/MPIFileBufJobDistributor.hh", "apps/public/rosetta_scripts/rosetta_scripts.cc" ), \
                ( "protocols/jd2/BOINCJobDistributor.hh", "apps/public/rosetta_scripts/rosetta_scripts.cc" ), \
                ( "protocols/jd2/BOINCJobDistributor.hh", "protocols/jd2/JobDistributorFactory.cc" ), \
                ( "protocols/jd2/MPIWorkPartitionJobDistributor.hh", "protocols/jd2/JobDistributorFactory.cc" ), \
                ( "protocols/jd2/MPIWorkPoolJobDistributor.hh", "protocols/jd2/JobDistributorFactory.cc" ), \
                ( "protocols/jd2/archive/MPIArchiveJobDistributor.hh", "protocols/jd2/JobDistributorFactory.cc" ), \
                ( "ObjexxFCL/FArray1D.hh", "core/kinematics/DomainMap.hh" ), \
                ( "numeric/deriv/dihedral_deriv.hh", "core/scoring/constraints/DihedralConstraint.cc" ),
                ( "iostream" "devel/domain_assembly/DomainAssemblyReader.cc", ),
                ( "core/scoring/mm/MMBondAngleResidueTypeParam.hh" "protocols/branch_angle/BranchAngleOptimizer.cc", ),
                ( "iostream" "protocols/genetic_algorithm/Entity.cc", ),
                ( "iostream" "protocols/pack_daemon/MultistateFitnessFunction.cc", ),
                ( "iostream" "protocols/smanager/smanager.cc", ),
                ( "iostream" "protocols/wum/WorkUnitManager.cc", ),
                ( "iostream" "protocols/optimize_weights/Arithmetic.cc", ),
                ( "iostream" "core/scoring/dssp/StrandPairing.cc", ),
                ( "core/id/DOF_ID_Range.hh" "protocols/moves/ThermodynamicMover.hh", ),
                ( "core/id/TorsionID_Range.hh" "protocols/moves/ThermodynamicMover.hh", ),
                ( "iostream" "core/pack/interaction_graph/DensePDInteractionGraph.cc", ),
                ( "core/optimization/MinimizerOptions.hh" "apps/benchmark/Minimizer.bench.hh", ),
                ( "iostream" "apps/pilot/will/test_string.cc", )]

      for pair in pairs :
         self.includes_to_leave_intact.add( pair )

   def initialize_includes_to_leave_be( self ) :
      self.includes_to_leave_be.add( "utility/vector1.hh" )
      self.includes_to_leave_be.add( "utility/vector0.hh" )
