Module Documentation
====================

.. automodule:: willclang
	:members:

Util
----
.. automodule:: willclang.util
	:members:
	
LibClangWrap
------------
.. automodule:: willclang.LibclangWrap
.. autofunction:: get_doctest_file
.. autofunction:: get_doctest_ast
.. autofunction:: print_raw_ast
.. autofunction:: hashcursor

SourceFile
----------
SourceFile class represents a single C++ sourcefile and should hold all data and logic
that doesn'd require parsing at AST
WORK IN PROGRESS!!!

.. autoclass:: SourceFile
	:members:

RosettaSourceFile
~~~~~~~~~~~~~~~~~
.. autoclass:: RosettaSourceFile
	:members:

AST
---
.. autoclass:: AST
	:members:
	
Node
~~~~
.. autoclass:: Node
	:members:

Exceptions
~~~~~~~~~~
.. autoexception:: ASTException
