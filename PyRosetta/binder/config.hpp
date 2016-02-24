// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   binder/config.hpp
/// @brief  Support for Binder Config file
/// @author Sergey Lyskov


#ifndef _INCLUDED_config_hpp_
#define _INCLUDED_config_hpp_

#include <string>
#include <vector>

namespace binder {

struct Config
{
	typedef std::string string;

	Config(string const &root_module_, std::vector<string> namespaces_to_bind_, std::vector<string> namespaces_to_skip_, string const &prefix_, uint maximum_file_length_) :
		root_module(root_module_), namespaces_to_bind(namespaces_to_bind_), namespaces_to_skip(namespaces_to_skip_), prefix(prefix_), maximum_file_length(maximum_file_length_) {}

	string root_module;

	std::vector<string> namespaces_to_bind, namespaces_to_skip;
	//std::vector<string> classes_to_bind, classes_to_skip;

	string prefix;

	uint maximum_file_length;


	/// check if user requested binding for given declaration
	bool is_namespace_binding_requested(string const &namespace_) const;
};


} // namespace binder

#endif // _INCLUDED_config_hpp_
