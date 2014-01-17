// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/properties/MultiChainSpanningTopology.hh
///
/// @brief      MultiChain Container for Spanning Topology Data
/// @details    Cacheable data object to the pose - acts as a container for spanning topology data to be
///             stored per-chain either continuously or noncontinuously.
///
/// @note       Last Modified: 1/2/14
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_properties_MultiChainSpanningTopology_hh
#define INCLUDED_core_membrane_properties_MultiChainSpanningTopology_hh

// Unit headers
#include <core/membrane/properties/MultiChainSpanningTopology.fwd.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>
#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <utility/vector0.hh>

// Platform headers
#include <platform/types.hh>

// C++ Headers
#include <string>
#include <cstdlib>
#include <map>

using namespace core;

namespace core {
namespace membrane {
namespace properties {
    
    /// @brief      MultiChain Container for Spanning Topology Data
    /// @details    Cacheable data object to the pose - acts as a container for spanning topology data to be
    ///             stored per-chain either continuously or noncontinuously.
    class MultiChainSpanningTopology : public basic::datacache::CacheableData {
        
    public: // constructors
        
        /// @brief Constructor
        MultiChainSpanningTopology();
        
        /// @brief Conpy Constructor
        MultiChainSpanningTopology( MultiChainSpanningTopology const & src );
        
        /// @brief Inherited Clone Method (Cacheable Data)
        basic::datacache::CacheableDataOP
        clone() const
        {
            return new MultiChainSpanningTopology( *this );
        }
        
        /// @brief Destructor
        ~MultiChainSpanningTopology();
        
    public: // getter/setter
        
        /// @brief Get Multi Chain Spanning Topology
        /// @details Maps one chain (int) to a spanning topology (SpanningTopology&)
        std::map< int, SpanningTopology& > get_multichain_topology();
        core::Size get_num_chains();
        
    public: // methods
        
        /// @brief Add By Chain
        /// @details Add a spanning topology by a given chain
        void add_by_chain( std::pair< int, SpanningTopology& > pair );
        
    private: // helper methods
        
        /// @brief Copy Data
        /// @details Copy Constructor Helper function
        void copy_data( MultiChainSpanningTopology src, MultiChainSpanningTopology copy );
        
    private: // data
        
        // Map Spanning topology to chain
        std::map< int, SpanningTopology& > multichain_topology_;
        
    }; // class MultiChainSpanningTopology
    
    /// @brief Add to Pose Cache Methods
    MultiChainSpanningTopology const & MultiChainSpanningTopology_from_pose( pose::Pose const & pose );
    MultiChainSpanningTopology & nonconst_MultiChainSpanningTopology_from_pose( pose::Pose & pose );
    
} // properties
} // membrane
} // core

#endif // INCLUDED_core_membrane_properties_MultiChainSpanningTopology_hh





