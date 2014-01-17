// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/membrane/properties/MultiChainSpanningTopology.fwd.hh
///
/// @brief      MultiChain Container for Spanning Topology Data
/// @details    Cacheable data object to the pose - acts as a container for spanning topology data to be
///             stored per-chain either continuously or noncontinuously.
///
/// @note       Last Modified: 1/2/14
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_properties_MultiChainSpanningTopology_fwd_hh
#define INCLUDED_core_membrane_properties_MultiChainSpanningTopology_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace membrane {
namespace properties {

    /// @brief      MultiChain Container for Spanning Topology Data
    /// @details    Cacheable data object to the pose - acts as a container for spanning topology data to be
    ///             stored per-chain either continuously or noncontinuously.
    class MultiChainSpanningTopology;
    typedef utility::pointer::owning_ptr< MultiChainSpanningTopology > MultiChainSpanningTopologyOP;
    typedef utility::pointer::owning_ptr< MultiChainSpanningTopology const > MultiChainSpanningTopologyCOP;
    
    
} // properties
} // membrane
} // core

#endif // INCLUDED_core_membrane_properties_MultiChainSpanningTopology_fwd_hh



