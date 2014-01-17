// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/membrane/MembraneProteinOptions.cc
/// @brief  options header for membrane protein options
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinOptions_cc
#define INCLUDED_core_membrane_MembraneProteinOptions_cc

//unit headers
#include <core/membrane/MembraneProteinOptions.hh>
#include <core/membrane/MembraneProteinOptionsCreator.hh>

//utility headers
#include <utility/tag/Tag.hh>

namespace core {
namespace membrane {
    
    /// @brief Return options type (string)
    std::string
    MembraneProteinOptionsCreator::options_type() const { return "MembraneProteinOptionsCreator"; }
    
    /// @brief Create a resource options class from parsed xml input
    basic::resource_manager::ResourceOptionsOP
    MembraneProteinOptionsCreator::create_options() const { return new MembraneProteinOptions; }
    
    /// @brief Constructor
    MembraneProteinOptions::MembraneProteinOptions() : membrane_( true ) {}
    
    /// @brief Destructor
    MembraneProteinOptions::~MembraneProteinOptions() {}
    
    /// @brief Parse XML File for options
    void
    MembraneProteinOptions::parse_my_tag(
                                   utility::tag::TagCOP tag
                                   )
    {
        membrane( tag->getOption< bool >( "membrane", true ));
    }
    
    /// @brief Return membrane protein options class type
    std::string
    MembraneProteinOptions::type() const
    {
        return "MembraneProteinOptions";
    }
    
    /// @brief Membrane Protein Options setters and getters
    bool MembraneProteinOptions::membrane() const { return membrane_; }
    void MembraneProteinOptions::membrane( bool setting ) { membrane_ = setting; }

} // namespace membrane
} // namespace core


#endif // INCLUDED_core_membrane_MembraneProteinOptions_cc

