// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//===- RosettaRefactorTool.hh - Clean up for boost ptr conversion -===//
//
// This tool generates a log of source code changes to be applies to
// the original source once all changes have been collected.
// The log is writted to stdout, one change per line.
//
//===-----------------------------------------------------------===//

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_RosettaRefactorTool_HH
#define INCLUDED_RosettaRefactorTool_HH

#include "clang/Basic/SourceManager.h"
#include "clang/Tooling/CompilationDatabase.h"
#include "clang/Tooling/Refactoring.h"
#include "clang/Tooling/Tooling.h"

#include <string>

class RosettaRefactorTool {

public:
	RosettaRefactorTool(int argc, const char **argv);
	int Run();
	int runMatchers();
	int saveOutput();

	std::unique_ptr<clang::tooling::CompilationDatabase> Compilations;
	clang::tooling::RefactoringTool * Tool;
	std::string prefix_, suffix_;
};

#endif
