set(LLVM_LINK_COMPONENTS
  Support
  )

add_clang_executable(rosetta-refactor-tool
  RosettaRefactorTool.cpp
  )

target_link_libraries(rosetta-refactor-tool
  clangAST
  clangASTMatchers
  clangBasic
  clangFrontend
  clangTooling
  )
