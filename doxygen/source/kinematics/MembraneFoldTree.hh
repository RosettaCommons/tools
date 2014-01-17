// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 	 core/membrane/kinematics/MembraneFoldTree.hh
///
/// @brief 	 Membrane Fold Tree
/// @details Sets up a fold tree with default membrane topology. This includes setting a membrane residue
///          as the root of the foldtree and adding jump edges from each virtual embedding residue
///          to the first residue in its respective chain. This class alsom maintains an extension from the edge
///          list that maps the embedding residues to their respective chains.
///
/// @note    Subclass FoldTree
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_kinematics_MembraneFoldTree_hh
#define INCLUDED_core_membrane_kimenatics_MembraneFoldTree_hh

// Unit Headers
#include <core/membrane/kinematics/MembraneFoldTree.fwd.hh>

// Project Headers
#include <core/membrane/geometry/MembraneResidueFactory.hh>
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

// Utility Headers
#include <core/types.hh>

#include <utility/vector1.hh>

#include <basic/Tracer.hh>

// // C++ Headers
#include <string>
#include <vector>
#include <iostream>

namespace core {
namespace membrane {
namespace kinematics {
            
    /// @brief   Membrane Fold Tree
    /// @details Constructs and maintains defualt foldtree topology and jumps for membrane proteins
    class MembraneFoldTree : public core::kinematics::FoldTree {
    
    public: // methods
        
        /// @brief    Default Constructor for Membrane FoldTree
        /// @details  Constructs and maintains a correct foldtree topology for membrane proteins
        ///
        /// @param      peptide_edges
        ///                 list of peptide edges to construct (peptide edges are per chain)
        /// @param      membrane_edges
        ///                 list of membrane edges to construct (also per chain)
        /// @param      root
        ///                 designated root
        ///
        /// @note   Precondition: peptide_edges.size() == membrane_edges.size()
        /// @note   Precondition: root >= 0
        /// @note   PostCondition: check_foldtree() == true
        /// @note   PostCondition: check_mp_foldtree() == true
        MembraneFoldTree(
                         utility::vector1< std::pair< int, int > > peptide_edges,
                         utility::vector1< std::pair< int, int > > membrane_edges,
                         int root
                         );
        
        /// @brief Destructor
        ~MembraneFoldTree();
        
        /// @brief Copy Constructor
        MembraneFoldTree( MembraneFoldTree const & src );
        
        /// @brief Assignment Operator
        MembraneFoldTree & operator=( MembraneFoldTree const & src );
        
        /// @brief      Check Membrane FoldTree
        /// @details    Check that the membrane foldtree maintains a valid membrane topology
        bool check_mpfoldtree();
        
        /// @brief      Get the root of the membrane foldtree
        /// @details    Returns the index of the memrbane residue which is the root of the fold tree
        int get_membrane_root();
        
        /// @brief      Get a map of each embedding residue to its respective chain
        /// @details    Note, the nth jump edge from membrane residue to embres is also the chain #
        std::map< int, int > get_embres_map();
        
    private: // standard constructor
        
        /// @brief      Emtpy Constructor
        /// @details    Extends standard constructor - creates a valid
        ///             fold tree but an invalid membrane fold tree
        MembraneFoldTree();
        
    private: // methods
        
        /// @brief      Build a membrane foldtree
        /// @details    Construct a membrane foldtree from a jump map and chain map
        void build_foldtree(
                            utility::vector1< std::pair< int, int > > peptide_edges,
                            utility::vector1< std::pair< int, int > > membrane_edges
                            );
        
        /// @brief Copy Data
        /// @details Copy Constructor Helper Function
        void copy_data( MembraneFoldTree src, MembraneFoldTree copy );
        
    private: // data
        
        // Store a reference to the membrane residue
        int membrane_root_;
        
        // Store a mapping from jump number to embedding residue index
        std::map< int, int> embres_map_;
        
    }; // class MembraneFoldTree
    
} // kinematics
} // membrane
} // core

#endif // INCLUDED_core_membrane_kinematics_MmebraneFoldTree_hh

