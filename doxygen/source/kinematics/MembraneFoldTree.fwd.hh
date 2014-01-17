// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 core/membrane/kinematics/MembraneFoldTree.fwd.hh
///
/// @brief 	 Membrane Fold Tree
/// @details Sets up a fold tree with default membrane topology. This includes setting a membrane residue
///          as the root of the foldtree and adding jump edges from each virtual embedding residue
///          to the first residue in its respective chain. This class alsom maintains an extension from the edge
///          list that maps the embedding residues to their respective chains.
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_kinematics_MembraneFoldTree_fwd_hh
#define INCLUDED_core_membrane_kimenatics_MembraneFoldTree_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh> 

namespace core {
namespace membrane {
namespace kinematics {

    /// @brief Membrane Fold Tree
    /// @details Constructs and maintains defualt foldtree topology and jumps for membrane proteins
    class MembraneFoldTree;
    typedef utility::pointer::owning_ptr< MembraneFoldTree > MembraneFoldTreeOP;
    typedef utility::pointer::owning_ptr< MembraneFoldTree const > MembraneFoldTreeCOP;
    
} // core
} // membrane
} // kinematics

#endif // INCLUDED_core_membrane_kinematics_MmebraneFoldTree_fwd_hh

