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
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_kinematics_MembraneFoldTree_hh
#define INCLUDED_core_membrane_kimenatics_MembraneFoldTree_hh

// Unit Header
#include <core/membrane/kinematics/MembraneFoldTree.hh>

// Project Headers
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

// Utility Headers
#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

// // C++ Headers
#include <string>
#include <vector>
#include <iostream>

using namespace core;
using namespace core::membrane::kinematics;

static basic::Tracer TR( "core.membrane.kinematics.MembraneFoldTree" );

namespace core {
namespace membrane {
namespace kinematics {
    
    /// Constructors ////////////////////
    
    /// @brief      Emtpy Constructor
    /// @details    Extends standard constructor - creates a valid
    ///             fold tree but an invalid membrane fold tree
    /// @note       Private
    MembraneFoldTree::MembraneFoldTree() :
        FoldTree(),
        membrane_root_(1)
    {}
    
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
    /// @note   Preconditions: nchains != 0
    /// @note   Precondition: root >= 0
    /// @note   PostCondition: check_foldtree() == true
    /// @note   PostCondition: check_mp_foldtree() == true
    MembraneFoldTree::MembraneFoldTree(
                     utility::vector1< std::pair< int, int > > peptide_edges,
                     utility::vector1< std::pair< int, int > > membrane_edges,
                     int root
            ) :
        FoldTree(),
        membrane_root_(1)
    {
        using namespace core::membrane::util;
        
        // Check Preconditions
        if ( root <= 0 ) {
            TR << "Exception type 1" << std::endl;
            throw new EXCN_MembraneFoldTree("Root out of bounds - cannot construct membrane foldtree!");
        }
        
        if ( peptide_edges.size() == 0 ) {
            TR << "Exception type 2" << std::endl;
            throw new EXCN_MembraneFoldTree("Cannot construct foldtree with an empty list of peptide edges");
        }
        
        if ( membrane_edges.size() == 0 ) {
            TR << "Exception type 3" << std::endl;
            throw new EXCN_MembraneFoldTree("Cannot construct foldtree with an empty list of membrane edges");
        }
        
        if ( peptide_edges.size() != membrane_edges.size() ) {
            TR << "Exception type 1" << std::endl;
            throw new EXCN_MembraneFoldTree("Must have one peptide edge in the map for every membrane edge in a correct membrane foldtree topology");
        }
        
        // Set the Root
        membrane_root_ = root;
        
        // build foldtree from inputs
        build_foldtree( peptide_edges, membrane_edges );
        
    }
    
    /// @brief Destructor
    MembraneFoldTree::~MembraneFoldTree() {}
    
    /// @brief Copy Constructor
    MembraneFoldTree::MembraneFoldTree( MembraneFoldTree const & src ) {
        copy_data( *this, src );
    }
    
    /// @brief Assignment Operator
    MembraneFoldTree & MembraneFoldTree::operator=( MembraneFoldTree const & src ) {
        
        // Abort self-assignment.
        if ( this == &src ) {
            return *this;
        }
        
        copy_data( *this, src );
        return *this;
    }
    
    //// Public Methods ////////////////////////////////////
    
    /// @brief      Check Membrane FoldTree
    /// @details    Check that the membrane foldtree maintains a valid membrane topology
    bool MembraneFoldTree::check_mpfoldtree() {
        
        // Check that we still have a valid foldtree
        if ( ! check_fold_tree() ) { return false; }
        
        // Check membrane root is root
        if ( ! is_root( membrane_root_ ) ) { return false; }
        
        return true;
        
    }
    
    /// @brief      Get the root of the membrane foldtree
    /// @details    Returns the index of the memrbane residue which is the root of the fold tree
    int MembraneFoldTree::get_membrane_root() {
        return membrane_root_;
    }
    
    /// @brief      Get a map of each embedding residue to its respective chain
    /// @details    Note, the nth jump edge from membrane residue to embres is also the chain #
    std::map< int, int > MembraneFoldTree::get_embres_map() {
        return embres_map_;
    }
    
    /// Private Methods /////////////////
    
    /// @brief      Build a membrane foldtree
    /// @details    Construct a membrane foldtree from a jump map and chain map
    void
    MembraneFoldTree::build_foldtree(
                                     utility::vector1< std::pair< int, int > > peptide_edges,
                                     utility::vector1< std::pair< int, int > > membrane_edges
                                     ) {
        
        // Get initial info
        int jump_counter = 1;
        int num_chains = (int) peptide_edges.size();
        
        // Add the membrane jump
        add_edge(1, membrane_root_, jump_counter);
        jump_counter++;
        
        // Construct peptide edges and appropriate chain jump edges
        for ( int i = 1; i <= num_chains; ++i ) {
            
            // Get info for peptide edge
            int const start_peptide = peptide_edges[i].first;
            int const stop_peptide = peptide_edges[i].second;
            
            // Add Edges
            add_edge( start_peptide, stop_peptide, -1 ); // peptide edge
            
            if ( i != num_chains ) {
                add_edge( stop_peptide, peptide_edges[i+1].first, jump_counter ); // temporarily changing chain conectivity so it is consecutive and maybe a correct tree??
                jump_counter++;
            }
        }
        
        // Construct membrane edges
        for ( int i = 1; i <= num_chains; ++i ) {
            
            // Get info for membrane edge
            int const start_membrane = membrane_edges[i].first;
            int const stop_membrane = membrane_edges[i].second;
            
            // Add edge
            embres_map_.insert( std::pair< int, int >( i, stop_membrane ) );
            add_edge( start_membrane, stop_membrane, jump_counter );
            jump_counter++;
        }
        
        //reorder( membrane_root_ );
    }
    
    /// @brief Copy Data
    /// @details Copy Constructor Helper Function
    void
    MembraneFoldTree::copy_data( MembraneFoldTree src, MembraneFoldTree copy ) {
        src.membrane_root_ = copy.membrane_root_;
        src.embres_map_ = copy.embres_map_;
    }
    
} // kinematics
} // membrane
} // core

#endif // INCLUDED_core_membrane_kinematics_MmebraneFoldTree_hh

