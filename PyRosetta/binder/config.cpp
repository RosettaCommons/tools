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

#include <config.hpp>

#include <util.hpp>

#include <stdexcept>

using std::string;

namespace binder {

bool Config::is_namespace_binding_requested(string const &namespace_) const
{
	bool to_bind_flag=false, to_skip_flag=false;
	string to_bind, to_skip;

	for(auto &n : namespaces_to_bind) {
		if( begins_wtih(namespace_, n) ) {
			if( n.size() > to_bind.size() ) { to_bind = n; to_bind_flag=true; }
		}
	}

	for(auto &s : namespaces_to_skip) {
		if( begins_wtih(namespace_, s) ) {
			if( s.size() > to_skip.size() ) { to_skip = s; to_skip_flag=true; }
		}
	}

	if( to_bind.size() > to_skip.size() ) return true;
	if( to_bind.size() < to_skip.size() ) return false;

	if( to_bind_flag and to_skip_flag ) throw std::runtime_error("Could not determent if namespace '" + namespace_ + "' should be binded or not... please check if options --bind and --skip conflicting!!!");

	if( to_bind_flag ) return true;

	return false;
}

} // namespace binder
