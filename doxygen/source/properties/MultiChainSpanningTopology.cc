// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/properties/MultiChainSpanningTopology.cc
///
/// @brief      MultiChain Container for Spanning Topology Data
/// @details    Cacheable data object to the pose - acts as a container for spanning topology data to be
///             stored per-chain either continuously or noncontinuously.
///
/// @note       Last Modified: 1/2/14
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_properties_MultiChainSpanningTopology_cc
#define INCLUDED_core_membrane_properties_MultiChainSpanningTopology_cc

// Unit headers
#include <core/membrane/properties/MultiChainSpanningTopology.hh>

// Project Headers
#include <core/membrane/properties/SpanningTopology.hh>
#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <utility/vector0.hh>

// Platform headers
#include <platform/types.hh>

// C++ Headers
#include <string>

using namespace core;

namespace core {
namespace membrane {
namespace properties {
        
    /// @brief Constructor
    MultiChainSpanningTopology::MultiChainSpanningTopology() :
        basic::datacache::CacheableData()
    {}
    
    /// @brief Conpy Constructor
    MultiChainSpanningTopology::MultiChainSpanningTopology( MultiChainSpanningTopology const & src ) :
        basic::datacache::CacheableData()
    {
        copy_data( src, *this );
    }

    /// @brief Destructor
    MultiChainSpanningTopology::~MultiChainSpanningTopology() {}

    /// @brief Get Multi Chain Spanning Topology
    std::map< int, SpanningTopology& >
    MultiChainSpanningTopology::get_multichain_topology() { return multichain_topology_; }
    
    /// @brief Get the numebr of chains stored in the obj
    core::Size
    MultiChainSpanningTopology::get_num_chains() { return multichain_topology_.size(); }
    
    /// @brief Add By Chain
    /// @details Add a new spanning topology object by a given chain (int)
    void
    MultiChainSpanningTopology::add_by_chain( std::pair< int, SpanningTopology& > pair ) {
        multichain_topology_.insert( pair );
    }
    
    /// @brief Copy Data
    /// @details Copy Constructor helper function
    void
    MultiChainSpanningTopology::copy_data( MultiChainSpanningTopology src, MultiChainSpanningTopology copy ) {
        src.multichain_topology_ = copy.multichain_topology_;
    }
    
    /// @brief Add to Pose Cache (Const)
    /// @details Add a const multichain spanning topology object to the pose cache
    MultiChainSpanningTopology const & MultiChainSpanningTopology_from_pose( pose::Pose const & pose ) {
        
        assert( pose.data().has( core::pose::datacache::CacheableDataType::MULTICHAIN_SPANNING_TOPOLOGY ) );
        return *( static_cast< MultiChainSpanningTopology const * >( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MULTICHAIN_SPANNING_TOPOLOGY )() ));
    }
    
    /// @brief Add to Pose Cache (nonconst)
    /// @details Add a nonconst multichain spanning topology object to the pose cache
    MultiChainSpanningTopology & nonconst_MultiChainSpanningTopology_from_pose( pose::Pose & pose ) {
        
        if ( pose.data().has( core::pose::datacache::CacheableDataType::MULTICHAIN_SPANNING_TOPOLOGY ) ) {
            return *( static_cast< MultiChainSpanningTopology * >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTICHAIN_SPANNING_TOPOLOGY )() ));
        }
        // else
        MultiChainSpanningTopologyOP topology = new MultiChainSpanningTopology;
        pose.data().set( core::pose::datacache::CacheableDataType::MULTICHAIN_SPANNING_TOPOLOGY, topology );
        return *topology;
        
    }
    
} // properties
} // membrane
} // core

#endif // INCLUDED_core_membrane_properties_MultiChainSpanningTopology_hh





