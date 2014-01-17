// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/membrane/MembraneProteinOptions.hh
/// @brief  options header for membrane protein options
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinOptions_hh
#define INCLUDED_core_membrane_MembraneProteinOptions_hh

//package headers
#include <basic/resource_manager/ResourceOptions.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

namespace core {
namespace membrane {
    
    /// @brief Membrane Protein Options Class
    class MembraneProteinOptions : public basic::resource_manager::ResourceOptions
    {
    public:
        
        /// @brief Constructor
        MembraneProteinOptions();
        
        /// @brief Destructor
        virtual ~MembraneProteinOptions();
        
        /// @brief Parse XML Tag (get options)
        virtual
        void
        parse_my_tag(
                     utility::tag::TagCOP tag
                     );
        
        /// @brief Return options class type
        virtual
        std::string
        type() const;
        
        /// @brief Membrane options setters and getters
        bool membrane() const;
        void membrane( bool setting );
        
        
        
    private:
        
        bool membrane_;
    };
        
} // namespace membrane
} // namespace core


#endif // INCLUDED_core_membrane_MembraneProteinOptions_hh

