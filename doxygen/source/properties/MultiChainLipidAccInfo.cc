// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/properties/MultiChainLipidAccInfo.cc
///
/// @brief      MultiChain Container for Lipid accessibility info
/// @details    Cacheable data object to the pose - stores lipid accessibility data by chain
///             for continuous and noncontinuous chains
///
/// @note       Last Modified: 1/2/14
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_properties_MultiChainLipidAccInfo_cc
#define INCLUDED_core_membrane_properties_MultiChainLipidAccInfo_cc

// Unit headers
#include <core/membrane/properties/MultiChainLipidAccInfo.hh>

// Project Headers
#include <core/membrane/properties/LipidAccInfo.hh>
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
    MultiChainLipidAccInfo::MultiChainLipidAccInfo() :
    basic::datacache::CacheableData()
    {}
    
    /// @brief Conpy Constructor
    MultiChainLipidAccInfo::MultiChainLipidAccInfo( MultiChainLipidAccInfo const & src ) :
    basic::datacache::CacheableData()
    {
        copy_data( src, *this );
    }
    
    /// @brief Destructor
    MultiChainLipidAccInfo::~MultiChainLipidAccInfo() {}
    
    /// @brief Get Lipid Accessibility Data Map
    /// @details Maps Chain to lipid Acc info object
    std::map< int, LipidAccInfo& >
    MultiChainLipidAccInfo::get_multichain_lips_exp() { return multichain_lips_exp_; }
    
    /// @brief Get number of chains
    /// @details Get number of chains stored by the multi chain map
    core::Size
    MultiChainLipidAccInfo::get_num_chains() { return multichain_lips_exp_.size(); }
    
    /// @brief Add By Chain
    /// @details Add a lipid acc data object by chain (int)
    void
    MultiChainLipidAccInfo::add_by_chain( std::pair< int, LipidAccInfo& > pair ) {
        multichain_lips_exp_.insert( pair );
    }
    
    /// @brief Copy Data
    /// @details Copy constructor helper function
    void
    MultiChainLipidAccInfo::copy_data( MultiChainLipidAccInfo src, MultiChainLipidAccInfo copy ) {
        src.multichain_lips_exp_ = copy.multichain_lips_exp_;
    }
    
    /// @brief Add to Pose Cache (const)
    /// @details Add a const data object to the pose
    MultiChainLipidAccInfo const & MultiChainLipidAccInfo_from_pose( pose::Pose const & pose ) {
        
        assert( pose.data().has( core::pose::datacache::CacheableDataType::MULTICHAIN_LIPID_ACC_DATA ) );
        return *( static_cast< MultiChainLipidAccInfo const * >( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MULTICHAIN_LIPID_ACC_DATA )() ));
    }
    
    /// @brief Add to pose cache (nonconst)
    /// @details Add a non const object to the pose cache
    MultiChainLipidAccInfo & nonconst_MultiChainLipidAccInfo_from_pose( pose::Pose & pose ) {
        
        if ( pose.data().has( core::pose::datacache::CacheableDataType::MULTICHAIN_LIPID_ACC_DATA ) ) {
            return *( static_cast< MultiChainLipidAccInfo * >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MULTICHAIN_LIPID_ACC_DATA )() ));
        }
        // else
        MultiChainLipidAccInfoOP lips_exp = new MultiChainLipidAccInfo;
        pose.data().set( core::pose::datacache::CacheableDataType::MULTICHAIN_LIPID_ACC_DATA, lips_exp );
        return *lips_exp;
        
    }
    
} // properties
} // membrane
} // core

#endif // INCLUDED_core_membrane_properties_MultiChainLipidAccInfo_cc





