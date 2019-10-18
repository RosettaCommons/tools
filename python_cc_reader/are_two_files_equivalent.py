from python_cc_reader.external.blargs import blargs
from python_cc_reader.beauty.are_two_files_equivalent import compare_whether_two_files_are_equivalent

with blargs.Parser(locals()) as p:
    p.str("file1")
    p.str("file2")

compare_whether_two_files_are_equivalent(file1, file2)
