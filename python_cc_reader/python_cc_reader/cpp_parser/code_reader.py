
# This class is intended to read a c++ file and
# determine which lines will be read by the compiler.
# That is, it keeps track of which lines end up inside
# a long-comment block, and which lines are inside
# ifdef blocks for variables that are not used
# The class should be given the opportunity to read
# each line in a file.  It maintains a stack of macro
# variables have been defined -- it also maintains
# a scope depth.  A scope depth of 0 is at global scope.

import re, sys
from .preprocessor_parser import PreProcParser

def strip_toendofline_comment(line):
    start = line.find("//")
    if start == -1:
        return line
    else:
        return line.partition("//")[0] + "\n"


class CodeReader:
    def __init__(self):
        self.re_pound_if = re.compile(r"\s*#if ")
        self.re_pound_ifdef = re.compile(r"\s*#ifdef")
        self.re_pound_ifndef = re.compile(r"\s*#ifndef")
        self.re_pound_else = re.compile(r"\s*#else")
        self.re_pound_endif = re.compile(r"\s*#endif")
        self.re_pound_define = re.compile(r"\s*#define")
        self.if_directive = "if_directive"
        self.re_if_directive = re.compile(self.if_directive)
        self.re_splitter = re.compile(r"\S+")
        self.long_comment = False
        self.nested_ifdefs = [("None", True)]
        self.defined_macros = []
        self.scope_level = 0
        self.file_name_stack = []
        self.line_num_stack = []
        self.full_line = ""
        self.commentless_line = ""
        self.preprocessor_parser = PreProcParser()

    def strip_line_comment(self, line):
        start = line.find("//")
        if start == -1:
            return line
        else:
            return line.partition("//")[0] + "\n"

    def line_is_visible(self):
        return self.commentless_line != "\n" and self.nested_ifdefs[-1][1]

    def curr_file(self):
        return self.file_name_stack[-1]

    def curr_line(self):
        return self.line_num_stack[-1]

    # tell the CodeReader that you're about to start processing a new file
    def push_new_file(self, filename):
        self.file_name_stack.append(filename)
        self.line_num_stack.append(0)

    def pop_last_file(self):
        self.file_name_stack.pop()
        self.line_num_stack.pop()

    # def recursive_commentless_line( self, line ) :
    #   if self.long_comment :
    #      if line.find( "*/" ) == -1 :
    #         return "\n"
    #      else :
    #         self.long_comment = False
    #         split_string = line.partition( "*/" )
    #         return self.recursive_commentless_line( split_string[ 2 ] )
    #   else :
    #      long_comment_start = line.find( "/*" )
    #      if long_comment_start == -1 :
    #         return line
    #      else :
    #         self.long_comment = True
    #         split_string = line.partition( "/*" )
    #         return split_string[ 0 ] + self.recursive_commentless_line( split_string[ 2 ] )

    def decomment_line(self, line):
        if len(line) == 0 or len(line) == 1:
            return line
        in_string = False
        if self.long_comment:
            newline = ""
        else:
            newline = line[0]
        last_char = line[0]
        for i in range(1, len(line)):
            if self.long_comment:
                if line[i] == "/" and last_char == "*":
                    self.long_comment = False
                    last_char = " "  # consider the "/" in "*/" consumed
                    continue
                # // in long comments is meaningless!
                # elif line[ i ] == "/" and last_char == "/" :
                #   return newline[:-1] + "\n" # remove the last /
            else:
                if in_string:
                    if (
                        line[i] == '"' and last_char != "\\"
                    ):  # escaped quotes don't end strings.
                        in_string = False
                    newline += line[i]
                else:
                    if line[i] == "/" and last_char == "/":
                        return newline[:-1] + "\n"  # remove the last /
                    elif line[i] == "*" and last_char == "/":
                        newline = newline[
                            :-1
                        ]  # delete the / at the beginning of the block comment
                        self.long_comment = True
                        last_char = (
                            " "
                        )  # don't let /*/ be read as a start and end of a block comment!
                        continue
                    else:
                        if line[i] == '"':
                            in_string = True
                        newline += line[i]
            last_char = line[i]
        return newline

    def determine_visibility(self):
        # print(self.long_comment, self.nested_ifdefs[-1][0], self.nested_ifdefs[-1][ 1], self.commentless_line,)
        if self.commentless_line == "\n":
            return
        if self.re_pound_if.match(self.commentless_line):
            if_expression = self.commentless_line.partition("#if")[2]
            ast = self.preprocessor_parser.tree_from_line(if_expression.strip())
            # print("if expression", if_expression, "ast", ast)
            if ast:
                for defnode in self.preprocessor_parser.defined_nodes:
                    defnode.defined_set = self.defined_macros
                next_visibility = ast.eval() and self.nested_ifdefs[-1][1]
                # print("visibility", ast.eval())
                self.nested_ifdefs.append(("parsed " + if_expression, next_visibility))

            else:
                # cannot process this "if" directive -- conservatively say that
                # both this and the "else" (if present) are visible.
                self.nested_ifdefs.append(
                    (
                        self.if_directive
                        + " "
                        + self.re_splitter.findall(self.commentless_line)[1],
                        True,
                    )
                )
        elif self.re_pound_ifdef.match(
            self.commentless_line
        ) or self.re_pound_ifndef.match(self.commentless_line):
            toks = self.re_splitter.findall(self.commentless_line)
            if len(toks) >= 2:
                name_list = toks[1].split("(")
                if len(name_list) >= 1:
                    name = name_list[0]
                    positive = self.re_pound_ifdef.match(self.commentless_line) != None
                    if (name in self.defined_macros) == positive:
                        # Logical AND: was the parent visible?
                        self.nested_ifdefs.append((name, self.nested_ifdefs[-1][1]))
                    else:
                        self.nested_ifdefs.append((name, False))
        elif self.re_pound_else.match(self.commentless_line):
            if self.re_if_directive.match(self.nested_ifdefs[-1][0]):
                # cannot process "if" directives -- leave the visibility set to True.
                return
            if len(self.nested_ifdefs) >= 2:
                last_ifdef = self.nested_ifdefs.pop()
                if last_ifdef[1]:
                    self.nested_ifdefs.append((last_ifdef[0], False))
                else:
                    # Logical AND: was the parent visible?
                    self.nested_ifdefs.append(
                        (last_ifdef[0], self.nested_ifdefs[-1][1])
                    )
            else:
                print("Error: in", self.file_name_stack[-1], " on line ")
                print(self.line_num_stack[-1], ":", self.full_line)
                print(self.line_num_stack[-1], ":", self.commentless_line)
                print(
                    "Encountered a #else directive when the size of the nested-ifdef stack is less than 2"
                )
                sys.exit(1)
        elif self.re_pound_endif.match(self.commentless_line):
            if len(self.nested_ifdefs) < 2:
                print("Error: in", self.file_name_stack[-1], " on line ")
                print(self.line_num_stack[-1], ":", self.full_line)
                print(self.line_num_stack[-1], ":", self.commentless_line)
                print(
                    "Encountered a #endif directive when the size of the nested-ifdef stack is less than 2"
                )
                sys.exit(1)
            self.nested_ifdefs.pop()
        elif self.nested_ifdefs[-1][1] and self.re_pound_define.match(self.commentless_line):
            toks = self.re_splitter.findall(self.commentless_line)
            if len(toks) < 2:
                print("Error: in", self.file_name_stack[-1], " on line ")
                print(self.line_num_stack[-1], ":", self.full_line)
                print(self.line_num_stack[-1], ":", self.commentless_line)
                print("#define directive should be followed by a macro variable name")
                sys.exit(1)
            self.defined_macros.append(toks[1])

    def stringless_line(self):
        # take the commentless line and remove all string and character literals
        newline = ""
        escape = False
        in_single = False
        in_double = False
        for i in range(len(self.commentless_line)):
            curr = self.commentless_line[i]
            if in_single or in_double:
                if curr == "\\":
                    if escape:
                        escape = False
                    else:
                        # print("ESCAPE", in_single, in_double)
                        escape = True
                    continue
            if escape:
                escape = False
                continue

            if in_single:
                if curr == "'":
                    in_single = False
            elif in_double:
                if curr == '"':
                    in_double = False
            else:
                if curr == "'":
                    in_single = True
                elif curr == '"':
                    in_double = True
                else:
                    newline += curr
        return newline

    def count_nesting(self):
        stringless = self.stringless_line()
        self.scope_level += stringless.count("{") - stringless.count("}")
        if self.scope_level < 0:
            print("Error: in", self.file_name_stack[-1], " on line ")
            print(self.line_num_stack[-1], ":", self.full_line)
            print(self.line_num_stack[-1], ":", self.commentless_line)
            print("Scoping level goes negative")
            sys.exit(1)

    # read the contents of the line, removing any comments
    # and then determine the visibility of the line
    def examine_line(self, line):
        self.line_num_stack[-1] += 1
        self.commentless_line = self.decomment_line(line)
        # print(self.long_comment, self.scope_level, self.commentless_line,)
        self.determine_visibility()
        if self.line_is_visible():
            self.count_nesting()


# This class is designed to keep track of class and function declarations as a header
# file is being read.
class HeavyCodeReader(CodeReader):
    def __init__(self):
        CodeReader.__init__(self)
        self.namespace_stack = []
        self.class_stack = []
        self.privacy_stack = []
        self.parent_class_name = None
        self.function_name = ""
        self.function_is_explicitly_virtual = False
        self.function_args = []
        self.function_arg_names = []
        self.function_return_type = ""
        self.function_body_present = False
        self.vartype = ""
        self.varname = ""
        self.statement_string = ""
        self.last_scope_level = 0
        self.found_lparen = False
        self.paren_depth = 0
        self.processing_function_declaration = False
        self.processing_class_declaration = False
        self.re_class_dec = re.compile(r"\s*class\s+")
        self.re_class_fwd_dec = re.compile(r"\s*class\s+\w+;")
        self.re_class_dataop_dec = re.compile(r"\s*[\w:]+OP\s+\w+;")
        self.re_optype = re.compile(r"[\w:]+OP$")
        self.re_namespace_dec = re.compile(r"\s*namespace\s+\w")
        self.reserve_words = set(
            [
                "using",
                "namespace",
                "class",
                "for",
                "while",
                "do",
                "repeat",
                "public",
                "private",
                "protected",
                "template",
                "typedef",
                "typename",
                "operator",
                "inline",
                "explicit",
                "static",
                "mutable",
                "virtual",
                "friend",
                "unsigned",
            ]
        )
        self.ignorable_statements = set(["typedef", "class", "using", "template"])
        self.skipable_statements = set(["static", "friend"])
        self.basic_types = set(["void", "char", "int", "short", "float", "double"])
        self.re_pound_directive = re.compile("#")
        # self.privacy_settings = [ "private", "protected", "public" ]
        self.re_privacy = [
            ("private", re.compile(r"\s*private\s*:\s*$")),
            ("protected", re.compile(r"\s*protected\s*:\s*$")),
            ("public", re.compile(r"\s*public\s*:\s*$")),
        ]
        self.re_assert_only = re.compile(r"(.*)ASSERT_ONLY\( (\w+) \)(.*)")
        self.re_whitespace = re.compile(r"\s+")
        self.re_whitespace_line = re.compile(r"\s*$")
        self.re_ends_in_operator_parens = re.compile(r"operator\s*\(\s*$")
        self.re_regular_varname = re.compile(r"\w:$")
        self.func_patterns = [
            ("varname", re.compile(r"[\w!:+=/$~.]+")),
            ("&", re.compile(r"&")),
            ("(", re.compile(r"\(")),
            (")", re.compile(r"\)")),
            ("*", re.compile(r"\*")),
            (",", re.compile(r",")),
            ("-", re.compile(r"-")),
            ("=", re.compile(r"=")),
            ('"', re.compile(r'"')),
            ("'", re.compile(r"'")),
            ("<", re.compile(r"<")),
            (">", re.compile(r">")),
            ("]", re.compile(r"\]")),
            ("[", re.compile(r"\[")),
        ]
        self.enum_pattern = re.compile(r"\s*enum\s")

    def examine_line(self, line):
        self.prep_for_newline()
        CodeReader.examine_line(self, line)
        if self.line_is_visible():
            # print(self.found_lparen, self.commentless_line,)
            self.examine_visible_line()
        # print(self.curr_line(), self.scope_name(), self.scope_level, self.vartype,)
        # print(self.processing_function_declaration, self.function_name, len( self.function_args ),)
        # print(self.function_body_present, self.commentless_line,)

    def at_class_scope(self):
        if len(self.class_stack) == 0:
            return False
        if (
            self.class_stack[-1][1]
            and self.last_scope_level == self.class_stack[-1][2] + 1
        ):
            return True
        return False

    def at_global_scope(self):
        if self.last_scope_level == 0:
            return True
        return False

    def just_parsed_a_class_function_declaration(self):
        if self.function_name == "" or self.processing_function_declaration:
            return False
        return True

    def just_parsed_a_class_variable_declaration(self):
        if self.vartype != "":
            return True

    def at_namespace_scope(self):
        if self.at_class_scope():
            return False
        if len(self.namespace_stack) == 0:
            return False
        if self.namespace_stack[-1][1] + 1 == self.last_scope_level:
            return True
        return False

    def current_class_name(self):
        if len(self.class_stack) == 0:
            return "none"
        if self.class_stack[-1][1] == False:
            return "none"
        return self.class_stack[-1][0]

    def current_class_parent_name(self):
        if len(self.class_stack) == 0:
            return None
        if self.class_stack[-1][1] == False:
            return None
        return self.class_stack[-1][3]

    def current_privacy_level(self):
        if len(self.privacy_stack) == 0:
            return "none"
        return self.privacy_stack[-1]

    def class_hasnt_begun(self):
        if len(self.class_stack) == 0:
            return False
        return not self.class_stack[-1][1]

    def remove_assert_only(self, line):
        retline = ""
        rematch = self.re_assert_only.search(line)
        if rematch:
            g1 = rematch.group(1)[:]
            g3 = rematch.group(3)[:]  # excise the ASSERT_ONLY( blah ) block
            retline = g1 + g3
        else:
            retline = line
        return retline

    def full_scope_name(self):
        if self.scope_level == 0:
            return ""
        scope = ""
        for ns in self.namespace_stack:
            if scope != "":
                scope += "::"
            scope += ns[0]
        for cname in self.class_stack:
            if scope != "":
                scope += "::"
            scope += cname[0]
        return scope

    def scope_name(self):
        if self.scope_level == 0:
            return "Global"
        elif self.at_class_scope():
            return self.class_stack[-1][0]
        elif self.at_namespace_scope():
            return self.namespace_stack[-1][0]
        else:
            return "other"

    def prep_for_newline(self):
        self.last_scope_level = self.scope_level
        if not self.processing_function_declaration:
            self.function_name = ""
            self.function_is_explicitly_virtual = False
            self.function_args = []
            self.function_arg_names = []
            self.function_return_type = ""
            self.function_body_present = False
            self.vartype = ""

    def examine_visible_line(self):
        if self.re_whitespace_line.match(self.commentless_line):
            return
        if self.re_pound_directive.match(self.commentless_line):
            return
        if self.at_class_scope():
            if self.scope_level == self.class_stack[-1][2]:
                self.class_stack.pop()
                self.privacy_stack.pop()
        if self.at_namespace_scope() or self.at_global_scope():
            if self.at_namespace_scope():
                if self.scope_level == self.namespace_stack[-1][1]:
                    self.namespace_stack.pop()
            if self.re_namespace_dec.match(self.commentless_line):
                ns = (
                    self.commentless_line.split("namespace")[1]
                    .strip()
                    .split(" ")[0]
                    .strip()
                )
                self.namespace_stack.append((ns, self.last_scope_level))
                return

        if self.re_class_dec.match(self.commentless_line):
            if self.class_hasnt_begun():
                # print("CLASS HASNT BEGUN")
                print("Discarding empty class", self.class_stack[-1][0])
                self.class_stack.pop()
                self.privacy_stack.pop()
                self.parent_class_name = None
            classname = (
                self.commentless_line.split("class")[1]
                .split(":")[0]
                .split("{")[0]
                .strip()
            )
            if self.re_class_fwd_dec.match(self.commentless_line):
                self.fwd_dec_class = classname
            else:
                class_has_begun = False
                pname = None
                if self.commentless_line.find("{") != -1:
                    class_has_begun = True

                if len(self.commentless_line.split(":")) > 1:
                    self.parent_class_name = (
                        self.commentless_line.partition(":")[2].split("{")[0].strip()
                    )
                if class_has_begun:
                    if self.parent_class_name:
                        self.strip_parent_class_privacy()
                        pname = self.parent_class_name
                        # print("   ",classname, "with base class", self.parent_class_name)
                    else:
                        pname = None
                        # print("   ",classname, "without a base class")
                    self.parent_class_name = None

                self.class_stack.append(
                    (classname, class_has_begun, self.last_scope_level, pname)
                )
                self.privacy_stack.append("private")
        elif self.class_hasnt_begun():
            parent_str = self.commentless_line.split("{")[0].strip()
            if self.parent_class_name:
                self.parent_class_name += parent_str
            else:
                self.parent_class_name = parent_str
            if (
                self.scope_level > self.class_stack[-1][2]
            ):  # we've found a "{" on this line
                last_class = self.class_stack.pop()

                if self.parent_class_name:  # figure out who my parent class is
                    self.strip_parent_class_privacy()
                    # print("   ", last_class[0], "with base class", self.parent_class_name)
                else:
                    pass
                    # print("   ", last_class[0], "without a base class")
                pname = self.parent_class_name

                new_class = (last_class[0], True, last_class[2], pname)
                self.class_stack.append(new_class)
                self.parent_class_name = None
        elif self.processing_function_declaration:
            # print("proc func dec")
            lcurly_ind = self.commentless_line.find("{")
            semicolon_ind = self.commentless_line.find(";")
            if lcurly_ind == -1 and semicolon_ind == -1:
                return
            self.processing_function_declaration = False
            if lcurly_ind == -1:
                self.function_body_present = False
            elif semicolon_ind == -1:
                self.function_body_present = True
            else:
                if semicolon_ind < lcurly_ind:
                    self.function_body_present = False
                else:
                    self.function_body_present = True

        elif self.at_class_scope():
            priv_statement = False
            for re_priv in self.re_privacy:
                if re_priv[1].match(self.commentless_line):
                    self.privacy_stack[-1] = re_priv[0]
                    if self.statement_string != "":
                        if self.statement_string.find("template") == -1:
                            print(
                                "ERROR: statement_string is not empty!",
                                len(self.statement_string),
                            )
                            print(self.statement_string)
                            assert self.statement_string == ""
                        else:
                            self.statement_string = ""
                    priv_statement = True
                    break
            if priv_statement:
                return
            assertonlyless_line = self.remove_assert_only(self.commentless_line)
            uptosemicolon_line = assertonlyless_line.split(";")[0]
            start_paren_depth = self.paren_depth
            if self.found_lparen == False:
                if uptosemicolon_line.find("(") != -1:
                    uptosemicolon_line = uptosemicolon_line.partition("(")[2]
                    self.paren_depth = 1
                    self.found_lparen = True
            self.paren_depth += uptosemicolon_line.count("(")
            self.paren_depth -= uptosemicolon_line.count(")")
            declaration_done = assertonlyless_line.find(";") != -1 or (
                self.found_lparen and self.paren_depth == 0
            )  # assertonlyless_line.find( ")" ) != -1
            # print(declaration_done)
            if declaration_done:
                # print("Statement string: ", self.statement_string)
                if not self.found_lparen:
                    # print("dec done, substr = ", assertonlyless_line.split(";")[ 0 ])
                    substr = assertonlyless_line.split(";")[0]
                else:
                    # print("Assertonly-less line:", assertonlyless_line,)
                    count_paren_depth = start_paren_depth
                    # print(count_paren_depth)
                    remaining = assertonlyless_line[:]
                    substr = ""
                    if count_paren_depth == 0:
                        # everything starts and ends on this line
                        assert remaining.find("(") != -1
                        substr += remaining[: remaining.find("(") + 1]
                        count_paren_depth = 1
                        remaining = remaining[remaining.find("(") + 1 :]
                    while count_paren_depth != 0:
                        # print("remaining:", remaining)
                        next_rparen_ind = remaining.find(")")
                        next_lparen_ind = remaining.find("(")
                        if next_rparen_ind == -1:
                            assert False  # there should be a ")" on this line
                        elif next_lparen_ind == -1:
                            substr += remaining[:next_rparen_ind]
                            remaining = remaining[next_rparen_ind:]
                            count_paren_depth -= 1
                        elif next_lparen_ind < next_rparen_ind:
                            substr += remaining[: next_lparen_ind + 1]
                            remaining = remaining[next_lparen_ind + 1 :]
                            count_paren_depth += 1
                        else:
                            if count_paren_depth == 1:
                                substr += remaining[
                                    :next_rparen_ind
                                ]  # clip the last ")"
                            else:
                                substr += remaining[: next_rparen_ind + 1]
                                remaining = remaining[next_rparen_ind + 1 :]
                            count_paren_depth -= 1
                    # print("substr:", substr)
            else:
                substr = assertonlyless_line
            is_operator_parens = False
            if self.re_ends_in_operator_parens.search(substr):
                is_operator_parens = True
                # print("OPERATOR ()", substr)
                remnant = assertonlyless_line.partition(")")[2]
                # print("remnant", remnant,)
                declaration_done = remnant.find(";") != -1 or (
                    self.found_lparen and self.paren_depth == 0
                )  # remnant.find( ")" ) != -1
                substr = substr + ")" + remnant.split(";")[0].split(")")[0]
                # print(declaration_done)
            last_tok = ""
            if declaration_done:
                semi_ind = assertonlyless_line.find(";")
                rparen_ind = assertonlyless_line.find(")")
                # print(semi_ind, rparen_ind)
                if semi_ind == -1:
                    last_tok = ")"
                elif rparen_ind == -1:
                    last_tok = ";"
                elif semi_ind < rparen_ind:
                    last_tok = ";"
                else:
                    last_tok = ")"
            self.statement_string += substr

            if declaration_done:
                if (
                    self.statement_string.strip() == "}"
                ):  # At the end of a class declaration, you see "};"
                    self.statement_string = ""
                    return
                # print(self.statement_string)
                self.found_lparen = False
                if last_tok == ")":
                    # function declaration

                    rest_of_line = assertonlyless_line.split(")")[1]
                    if rest_of_line.find("{") != -1 or rest_of_line.find(";") != -1:
                        lcurly_ind = self.commentless_line.find("{")
                        semicolon_ind = self.commentless_line.find(";")
                        if lcurly_ind == -1:
                            self.function_body_present = False
                        elif semicolon_ind == -1:
                            self.function_body_present = True
                        else:
                            if semicolon_ind < lcurly_ind:
                                self.function_body_present = False
                            else:
                                self.function_body_present = True

                        self.processing_function_declaration = False
                    else:
                        self.processing_function_declaration = True
                    toks = self.tokenize_func_dec_string(self.statement_string)
                    # print("TOKS:")
                    # for tok in toks :
                    #   print(tok)
                    count = -1
                    ind_of_lparen = -1
                    for tok in toks:
                        count += 1
                        if tok == "(":
                            ind_of_lparen = count
                            if ind_of_lparen > 0:
                                if (
                                    toks[ind_of_lparen - 1] != "operator"
                                ):  # stop now if we're not parsing operator()
                                    # print("Found lparen: ", ind_of_lparen)
                                    break
                    if ind_of_lparen == -1:
                        print(
                            "Error: failed to find left paren while parsing function in class"
                        )
                        print(self.class_stack[-1][0])
                        self.statement_string = ""
                        return
                    if ind_of_lparen == 0:
                        print(
                            "Error: failed to find function name before left parenthesis while parsing function in class"
                        )
                        print(self.class_stack[-1][0])
                        self.statement_string = ""
                        return
                    # print("FUNC TOKS:")
                    # for tok in toks :
                    #   print("   ", tok)
                    self.function_name = toks[ind_of_lparen - 1]
                    # print("FUNCTION", self.function_name)
                    function_name_begin_ind = ind_of_lparen - 1
                    if not self.re_regular_varname.match(self.function_name):
                        if ind_of_lparen > 1:
                            if toks[ind_of_lparen - 2] == "operator":
                                function_name_begin_ind = ind_of_lparen - 2
                                self.function_name = (
                                    "operator" + toks[ind_of_lparen - 1]
                                )
                            elif ind_of_lparen > 2:
                                if toks[ind_of_lparen - 3] == "operator":
                                    self.function_name = (
                                        "operator"
                                        + toks[ind_of_lparen - 2]
                                        + toks[ind_of_lparen - 1]
                                    )
                                    function_name_begin_ind = ind_of_lparen - 3
                    if self.function_name == ">":
                        # ok, we're looking at a case like "template < class T > T funcname<> ( T const & one )"
                        # where the function name is actually ind_of_lparen-3
                        self.function_name = toks[ind_of_lparen - 3]
                        function_name_begin_ind = ind_of_lparen - 3
                        # print("template function function_name_begin_ind reset: ", self.function_name)
                    dstor_string = "~" + self.class_stack[-1][0]
                    # print(dstor_string)
                    if self.function_name.find(dstor_string) != -1:
                        self.function_name = "Destructor"
                    elif self.function_name.find(self.class_stack[-1][0]) != -1:
                        self.function_name = "Constructor"
                    return_type_toks = toks[:function_name_begin_ind]
                    if len(return_type_toks) != 0:
                        self.function_return_type, dummy = self.interpret_toks_as_funcparam(
                            return_type_toks
                        )
                        # print(self.function_return_type, dummy)
                        if dummy != "":
                            print("ERROR: dummy return type is not empty!")
                            print("FUNCTION NAME: ", self.function_name)
                            print("RETURN TYPE: ", self.function_return_type)
                            print("DUMMY:", dummy)
                            print("LINE:", self.curr_line())
                            assert dummy == ""
                    else:
                        self.function_return_type = ""
                    count = ind_of_lparen
                    start = ind_of_lparen + 1
                    self.function_args = []
                    self.function_arg_names = []
                    while count < len(toks) - 1:
                        count += 1
                        if toks[count] == ",":
                            paramtype, paramname = self.interpret_toks_as_funcparam(
                                toks[start:count]
                            )
                            if paramtype:
                                self.function_args.append(paramtype)
                                self.function_arg_names.append(paramname)
                            start = count + 1
                    paramtype, paramname = self.interpret_toks_as_funcparam(
                        toks[start:count]
                    )
                    if paramtype:
                        self.function_args.append(paramtype)
                        self.function_arg_names.append(paramname)
                    self.statement_string = ""
                elif self.enum_pattern.match(self.statement_string):
                    # don't read enums listed in classes -- clear out the statement string.  I don't know if this will work for multi-line enums
                    print("enum encountered; skipping", self.statement_string)
                    self.statement_string = ""
                    pass
                else:  # variable declaration
                    toks = self.tokenize_func_dec_string(self.statement_string)
                    # print("VAR TOKS:")
                    # for tok in toks :
                    #   print("   ",tok)
                    tok_to_read = 0
                    while tok_to_read != len(toks):
                        if toks[tok_to_read] in self.skipable_statements:
                            tok_to_read += 1
                        elif toks[tok_to_read] in self.ignorable_statements:
                            self.statement_string = ""
                            return
                        else:
                            break
                    vartype, varname = self.interpret_toks_as_vardec(toks)
                    self.vartype = vartype
                    self.varname = varname
                    self.statement_string = ""
        else:
            self.statement_string = ""

    def tokenize_func_dec_string(self, instring):
        strcopy = instring[:]
        retlist = []
        re_whitespace = re.compile(r"\s+")
        while strcopy:
            goon = False
            # print(strcopy)
            ws_match = None
            ws_match = self.re_whitespace.match(strcopy)
            if ws_match:
                strcopy = strcopy[ws_match.end() :]
            if strcopy == "":
                break
            for tokentype, pattern in self.func_patterns:
                m = pattern.match(strcopy)
                if m:
                    if tokentype == "varname":
                        # print("varname", strcopy[ : m.end() ])
                        retlist.append(strcopy[: m.end()])
                    else:
                        retlist.append(tokentype)
                    strcopy = strcopy[m.end() :]
                    goon = True
                    break
            if goon:
                continue
            print("ERROR: while loop failed to terminate while processing string")
            print("instring:", instring)
            print("strcopy:", strcopy)
            break

        return retlist

    def interpret_toks_as_funcparam(self, toks):
        # print(toks)
        paramtype = ""
        paramname = ""
        if len(toks) == 0:
            return None, None
        count = -1
        template_depth = 0
        depth_0_type = ""
        n_valid_toks = len(toks)
        for tok in toks:
            count += 1
            # print("TOK: ", tok, "TYPE:", paramtype, "NAME:", paramname, template_depth, depth_0_type)
            if tok in self.reserve_words:
                n_valid_toks -= 1
                if tok == "virtual":
                    self.function_is_explicitly_virtual = True
                continue
            paramtype += " "
            if tok == "<":
                template_depth += 1
                paramtype += tok
            elif tok == ">":
                template_depth -= 1
                paramtype += tok
            elif tok == "const" or tok == "*" or tok == "&":
                paramtype += tok
            elif tok == "=":
                break
            elif template_depth != 0:
                paramtype += tok
            elif depth_0_type == "":
                paramtype += tok
                depth_0_type = tok
            else:
                assigned_to_paramtype = False
                if count > 0:
                    if toks[count - 1] == ">" and tok[0] == ":":
                        paramtype += tok
                        assigned_to_paramtype = True
                if not assigned_to_paramtype:
                    paramname = tok
        return paramtype.strip(), paramname

    def interpret_toks_as_vardec(self, toks):
        vartype, varname = self.interpret_toks_as_funcparam(toks)
        if varname == "":
            print("ERROR: interpret_toks_as_vardec failed to find a variable name.")
            for tok in toks:
                print(tok)
            print()
            varname = "DUMMY"
        return vartype, varname

    def strip_parent_class_privacy(self):
        parts = self.parent_class_name.partition(" ")
        if parts[0] == "virtual":
            self.parent_class_name = parts[2]
            self.strip_parent_class_privacy()
        if parts[0] == "public" or parts[0] == "private" or parts[0] == "protected":
            self.parent_class_name = parts[2]


class FunctionDeclaration:
    def __init__(self):
        self.name = ""
        self.is_explicitly_virtual = False
        self.return_type = ""
        self.parameters = ""
        self.definition_present = True
        self.definition_begin_line = -1
        self.definition_end_line = -1
        self.privacy_level = "unk"


class ClassDeclaration:
    def __init__(self):
        self.name = ""
        self.base = None
        self.scope = ""
        self.is_templated = False
        self.functions = []
        self.data_members = []  # a list of ordered pairs, where [0] is the type, and [1] is the name of the variable
        self.file_declared_within = ""


class HCROptions:
    def __init__(self):
        self.defined_preprocessor_variables = []


def read_classes_from_header(fname, flines, options=None):
    classes = []
    cr = HeavyCodeReader()
    if options:
        cr.defined_macros = list(options.defined_preprocessor_variables)
    cr.push_new_file(fname)
    curr_class = None
    func = None
    scope_defining_function_end = -1
    for line in flines:
        cr.examine_line(line)

        if cr.at_class_scope():
            if not curr_class or (
                curr_class and curr_class.name != cr.current_class_name()
            ):
                if curr_class:
                    if func:
                        curr_class.functions.append(func)
                        func = None
                    classes.append(curr_class)
                curr_class = ClassDeclaration()
                curr_class.name = cr.current_class_name()
                curr_class.base = cr.current_class_parent_name()
                curr_class.scope = cr.full_scope_name()
                curr_class.file_declared_within = fname
        if cr.just_parsed_a_class_function_declaration():
            # print("Parsed function", cr.function_name)
            if func:
                curr_class.functions.append(func)
            func = FunctionDeclaration()
            # deep cpy
            func.name = str(cr.function_name)
            func.is_explicitly_virtual = cr.function_is_explicitly_virtual
            func.return_type = str(cr.function_return_type)
            func.parameters = list(cr.function_args)
            func.definition_present = cr.function_body_present
            # print("def present?", func.definition_present)
            func.privacy_level = cr.current_privacy_level()
            if not cr.function_body_present:
                curr_class.functions.append(func)
                func = None
            else:
                func.definition_begin_line = cr.curr_line()
                scope_defining_function_end = cr.last_scope_level

        if cr.just_parsed_a_class_variable_declaration():
            curr_class.data_members.append((cr.vartype, cr.varname))

        if scope_defining_function_end >= cr.scope_level:
            assert func
            func.definition_end_line = cr.curr_line()
            curr_class.functions.append(func)
            func = None
            scope_defining_function_end = -1
    if curr_class:
        if func:
            curr_class.functions.append(func)
        classes.append(curr_class)
    return classes


def find_data_members_assigned_in_assignment_operator(fname, flines, class_data):
    asgnop_regex = re.compile(class_data.name + r"::operator[\s]*=[\s]*\(")
    inline_asgnop_regex = re.compile(r"operator[\s]*=[\s]*\(")
    assignment_regex = re.compile("[^!]=[^=]")
    cr = HeavyCodeReader()
    cr.push_new_file(fname)
    reading_asgnop = False
    asgnop_begun = False
    data_members_assigned = set([])
    for line in flines:
        cr.examine_line(line)
        if cr.at_class_scope() and cr.current_class_name() == class_data.name:
            if inline_asgnop_regex.search(cr.commentless_line):
                reading_asgnop = True
                ns_scope = cr.last_scope_level
        elif asgnop_regex.search(cr.commentless_line) and cr.at_namespace_scope():
            reading_asgnop = True
            ns_scope = cr.namespace_stack[-1][1] + 1
            # print(ns_scope, cr.last_scope_level, cr.scope_level, cr.at_namespace_scope(), line)
        if reading_asgnop:
            # print(cr.scope_level, cr.commentless_line,)
            if cr.scope_level > ns_scope:
                asgnop_begun = True
                if assignment_regex.search(cr.commentless_line):
                    lhs = (
                        cr.commentless_line.rpartition("=")[0]
                        + cr.commentless_line.rpartition("=")[1]
                    )
                    # print(cr.scope_level, cr.commentless_line, lhs)
                    for dm in class_data.data_members:
                        # print("looking for datamember", dm[1], "on line", lhs)
                        restring = r"[\s.>]*" + dm[1] + r"[\s=[(]"
                        # print(restring)
                        if re.search(restring, lhs):
                            # print("  found data member ", dm[ 1 ], "on the left hand side of the assignment statement: ", lhs)
                            data_members_assigned.add(dm[1])
            else:
                if asgnop_begun:
                    reading_asgnop = False
    if not asgnop_begun:
        for func in class_data.functions:
            if func.name == "operator=" and not func.privacy_level == "private":
                print(
                    "Did not find assgnment operator for",
                    class_data.name,
                    " in ",
                    fname,
                )
        return None
    return data_members_assigned


def find_data_members_assigned_in_copy_constructor(fname, flines, class_data):
    cname = class_data.name
    restring = cname + r"[\s]*::[\s]*" + cname + r"[\s]*\([^,]*" + cname + r"[^,]*\)"
    # print(restring)
    copy_ctor_regex = re.compile(restring)
    cr = HeavyCodeReader()
    cr.push_new_file(fname)
    on_own_line = re.compile(r"[;\)}]")
    this_equals_re = re.compile(r"\*this[\s]*=[^=]")
    this_equals_re2 = re.compile(r"\([\s]*\*this[\s]*\)[\s]*=[^=]")
    opequals_re = re.compile(r"operator[\s]*=[\s]*\(")

    last_line = ""
    processing_copy_ctor = False
    cc_scope_begun = False

    dmass = None

    for line in flines:
        cr.examine_line(line)
        if processing_copy_ctor:
            # print(processing_copy_ctor, cc_scope_begun, line,)
            if cr.at_namespace_scope():
                if cc_scope_begun:
                    processing_copy_ctor = False
                    continue
                else:
                    for dm in class_data.data_members:
                        if re.search(r"[\s]*" + dm[1] + r"\(", cr.commentless_line):
                            # print("   found data member", dm[1], "initialized in the intialize list for", cname)
                            dmass.add(dm[1])
                    if cr.commentless_line.find("{") != -1:
                        cc_scope_begun = True
            else:
                if (
                    this_equals_re.search(cr.commentless_line)
                    or opequals_re.search(cr.commentless_line)
                    or this_equals_re2.search(cr.commentless_line)
                ):
                    # print("Found an invocation of the assignment operator!", cr.commentless_line,)
                    # for dm in dmass :
                    #   print("before", dm)
                    dmassop = find_data_members_assigned_in_assignment_operator(
                        fname, flines, class_data
                    )
                    if dmassop:
                        dmass |= dmassop
                    # for dm in dmass :
                    #   print("after", dm)
                for dm in class_data.data_members:
                    # print("looking for datamember", dm[1], "on line", cr.commentless_line,)
                    restring = r"[\s.>]*" + dm[1] + r"[\s=[(]"
                    # print(restring)
                    if re.search(restring, cr.commentless_line):
                        # print("  found data member ", dm[ 1 ], "on the left hand side of the assignment statement: ", cr.commentless_line,)
                        dmass.add(dm[1])

        else:

            if (
                cr.paren_depth != 0 or cr.at_namespace_scope()
            ) and not on_own_line.search(cr.commentless_line):
                last_line += " " + cr.commentless_line.strip()
            else:
                last_line += " " + cr.commentless_line.strip()
                # print("one line", last_line)
                if copy_ctor_regex.search(last_line):
                    # print("DMASS CREATION", last_line)
                    dmass = set([])
                    processing_copy_ctor = True
                    if cr.commentless_line.find("{") != -1:
                        cc_scope_begun = True
                last_line = ""
    return dmass


def find_data_members_unassigned_in_assignment_operator(fname, flines, class_data):
    assigned = find_data_members_assigned_in_assignment_operator(
        fname, flines, class_data
    )
    unassigned = []
    if assigned == None:
        return unassigned
    for dm in class_data.data_members:
        if not dm[1] in assigned:
            unassigned.append(dm[1])
    return unassigned


def find_data_members_unassigned_in_copy_constructor(fname, flines, class_data):
    assigned = find_data_members_assigned_in_copy_constructor(fname, flines, class_data)
    # for ass in assigned :
    #   print(ass, "is assigned in copy constructor")
    unassigned = []
    if assigned == None:
        return unassigned
    for dm in class_data.data_members:
        if not dm[1] in assigned:
            unassigned.append(dm[1])
    return unassigned
