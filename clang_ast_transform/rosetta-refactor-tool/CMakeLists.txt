set(LLVM_LINK_COMPONENTS
  Support
  )

#add_subdirectory( matchers )

# add_clang_library( BasicRosettaMatchers
#   matchers_base
#   utils
# 
#   LINK_LIBS
#   clangBasic
#   clangASTMatchers
#   LLVMSupport
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
  matchers/find/calls
  matchers/find/constructor_decl
  matchers/find/field_decl
  matchers/find/record_decl
  matchers/find/self_ptr_in_ctor
  matchers/find/serialization_funcs
)

target_link_libraries(rosetta-refactor-tool
  LLVMSupport
  clangAST
  clangASTMatchers
  clangBasic
  clangFrontend
  clangTooling
  #BasicRosettaMatchers
  #RosettaFinders
  )
