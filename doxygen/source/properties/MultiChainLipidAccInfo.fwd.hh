// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/properties/MultiChainLipidAccInfo.fwd.hh
///
/// @brief      MultiChain Container for Lipid accessibility info
/// @details    Cacheable data object to the pose - stores lipid accessibility data by chain
///             for continuous and noncontinuous chains
///
/// @note       Last Modified: 1/2/14
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_properties_MultiChainLipidAccInfo_fwd_hh
#define INCLUDED_core_membrane_properties_MultiChainLipidAccInfo_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace properties {
    
    /// @brief      MultiChain Container for Lipid accessibility info
    /// @details    Cacheable data object to the pose - stores lipid accessibility data by chain
    ///             for continuous and noncontinuous chains
    class MultiChainLipidAccInfo;
    typedef utility::pointer::owning_ptr< MultiChainLipidAccInfo > MultiChainLipidAccInfoOP;
    typedef utility::pointer::owning_ptr< MultiChainLipidAccInfo const > MultiChainLipidAccInfoCOP;
    
    
} // properties
} // membrane
} // core

#endif // INCLUDED_core_membrane_properties_MultiChainLipidAccInfo_fwd_hh

