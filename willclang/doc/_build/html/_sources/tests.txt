Tests
=====

doctest
-------

http://docs.python.org/library/doctest.html

doctest is really sweet... just put example code in your docstrings, and doctest will run it. all the 83 tests I've written today are done this way. It's easy and they double as nice documentation! 

using doctest
~~~~~~~~~~~~~

If you put the following at the bottom of *any* python file (even some deep API only thing) it'll run anything that looks like code in your docstrings in that file and compare the output to what

to doctest::

	if __name__ == '__main__':
		import doctest
		tr = doctest.testmod()
		print "tests passed:",tr.attempted-tr.failed
		print "tests failed:",tr.failed


unittests
---------

options
~~~~~~~

* nose_? 
* py.test_? 
* unittest2_? 
* or is doctest_ sufficient?

.. _nose: http://code.google.com/p/python-nose/
.. _py.test: http://pytest.org/latest/
.. _unittest2: http://docs.python.org/library/unittest.html
.. _doctest: http://docs.python.org/library/doctest.html
