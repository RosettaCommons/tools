set(LLVM_LINK_COMPONENTS
  Support
  )

set(LLVM_USED_LIBS clangTooling clangBasic clangAST clangASTMatchers )

add_library( RosettaFinders
  serialization_funcs
  )

target_link_libraries( RosettaFinders
  LLVMSupport
  clangASTMatchers
)