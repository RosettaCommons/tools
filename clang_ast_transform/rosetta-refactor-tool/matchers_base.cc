// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include "matchers_base.hh"
#include <iostream>

ReplaceMatchCallback::ReplaceMatchCallback(clang::tooling::Replacements *Replace, const char *tag)
	: Replace(Replace), tag(tag), debug_( false )
{
	std::cout << "Setting debug to false" << std::endl;
}

ReplaceMatchCallback::ReplaceMatchCallback(clang::tooling::Replacements *Replace, const char *tag, bool debug)
	: Replace(Replace), tag(tag), debug_( debug )
{}

