set(LLVM_LINK_COMPONENTS support)
set(LLVM_USED_LIBS clangTooling clangBasic clangAST)

add_clang_executable(binder
  binder.cpp
  )
target_link_libraries(binder
  clangTooling
  clangBasic
  clangASTMatchers
  )
