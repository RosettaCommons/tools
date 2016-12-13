set(LLVM_LINK_COMPONENTS
  Support
  )

# add_subdirectory( matchers )

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
  ../matchers_base
  ../utils
  ../matchers/code_quality/bad_pointer_casts
  ../matchers/code_quality/naked_ptr_op_casts
  ../matchers/code_quality/obj_on_stack
  ../matchers/find/calls
  ../matchers/find/constructor_decl
  ../matchers/find/field_decl
  ../matchers/find/record_decl
  ../matchers/find/self_ptr_in_ctor
  ../matchers/find/serialization_funcs
  ../matchers/rewrite/add_serialization_code
  ../matchers/rewrite/call_operator
  ../matchers/rewrite/cast_from_new
  ../matchers/rewrite/cast_from_new_expr
  ../matchers/rewrite/cast_from_new_vardecl
  ../matchers/rewrite/cast_in_assignment
  ../matchers/rewrite/ctor_initializer
  ../matchers/rewrite/datamap_get
  ../matchers/rewrite/dynamic_cast
  ../matchers/rewrite/member_calls
  ../matchers/rewrite/pointer_name
  ../matchers/rewrite/pose_dynamic_cast
  ../matchers/rewrite/real_comparison
  ../matchers/rewrite/typedef
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

