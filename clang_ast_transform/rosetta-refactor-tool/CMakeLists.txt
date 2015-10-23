set(LLVM_LINK_COMPONENTS
  Support
  )

# add_subdirectory( matchers )

# add_library( BasicRosettaMatchers
#   matchers_base
#   utils
#   )

# target_link_libraries( BasicRosettaMatchers
#   LLVMSupport
#   clangAST
#   clangASTMatchers
#   clangBasic
#   clangFrontend
#   clangTooling
#   )

add_clang_executable(rosetta-refactor-tool
  RosettaRefactorTool.cpp
  matchers_base
  utils
  matchers/find/serialization_funcs
  )

target_link_libraries(rosetta-refactor-tool
  LLVMSupport
  clangAST
  clangASTMatchers
  clangBasic
  clangFrontend
  clangTooling
  )
