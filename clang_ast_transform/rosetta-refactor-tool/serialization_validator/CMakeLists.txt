set(LLVM_LINK_COMPONENTS
  Support
  )


add_clang_executable( serialization_validator
  all_members_serialized_validator.cc
  ../matchers/find/serialization_funcs
  ../matchers/find/field_decl
)

target_link_libraries( serialization_validator
  LLVMSupport
  clangAST
  clangASTMatchers
  clangBasic
  clangFrontend
  clangTooling
  #BasicRosettaMatchers
  #RosettaFinders
  )