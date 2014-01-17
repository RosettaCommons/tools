// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/properties/MultiChainLipidAccInfo.hh
///
/// @brief      MultiChain Container for Lipid accessibility info
/// @details    Cacheable data object to the pose - stores lipid accessibility data by chain
///             for continuous and noncontinuous chains
///
/// @note       Last Modified: 1/2/14
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_properties_MultiChainLipidAccInfo_hh
#define INCLUDED_core_membrane_properties_MultiChainLipidAccInfo_hh

// Unit headers
#include <core/membrane/properties/MultiChainLipidAccInfo.fwd.hh>

// Project Headers
#include <core/membrane/properties/LipidAccInfo.hh>
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
    
    /// @brief      MultiChain Container for Lipid accessibility info
    /// @details    Cacheable data object to the pose - stores lipid accessibility data by chain
    ///             for continuous and noncontinuous chains
    class MultiChainLipidAccInfo : public basic::datacache::CacheableData {
        
    public: // constructors
        
        /// @brief Constructor
        MultiChainLipidAccInfo();
        
        /// @brief Conpy Constructor
        MultiChainLipidAccInfo( MultiChainLipidAccInfo const & src );
        
        /// @brief Inherited Clone Method
        /// @details Clone data as cacheable object
        basic::datacache::CacheableDataOP
        clone() const
        {
            return new MultiChainLipidAccInfo( *this );
        }
        
        /// @brief Destructor
        ~MultiChainLipidAccInfo();
        
    public: // getter/setter
        
        /// @brief Get Lipid Acc Info Multichain map
        std::map< int, LipidAccInfo& > get_multichain_lips_exp();
        
        /// @brief Get number of chains in the multichain object
        core::Size get_num_chains();
        
    public: // methods
        
        /// @brief Add By Chain
        void add_by_chain( std::pair< int, LipidAccInfo& > pair );
        
    private: // helper methods
        
        /// @brief Copy Data
        void copy_data( MultiChainLipidAccInfo src, MultiChainLipidAccInfo copy );
        
    private: // data
        
        // Map Spanning topology to chain
        std::map< int, LipidAccInfo& > multichain_lips_exp_;
        
    }; // class MultiChainLipidAccInfo
    
    /// @brief Add to Pose Cache (const and nonconst)
    MultiChainLipidAccInfo const & MultiChainLipidAccInfo_from_pose( pose::Pose const & pose );
    MultiChainLipidAccInfo & nonconst_MultiChainLipidAccInfo_from_pose( pose::Pose & pose );
    
} // properties
} // membrane
} // core

#endif // INCLUDED_core_membrane_properties_MultiChainLipidAccInfo_hh





