set(LLVM_LINK_COMPONENTS
  Support
  )

set(LLVM_USED_LIBS clangTooling clangBasic clangAST clangASTMatchers )

add_clang_library( RosettaFinders
  serialization_funcs

  LINK_LIBS
  BasicRosettaMatchers
  clangASTMatchers
  LLVMSupport
  )

# target_link_libraries( RosettaFinders
#   LLVMSupport
#   clangASTMatchers
# )