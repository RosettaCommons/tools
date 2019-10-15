from enum import Enum


# partial support for preprocessor conditional logic:
# allows for arbitrary combinations of
# 1. !
# 2. ()
# 3. defined
# 4. variable-names
# 5. &&
# 6. ||
# 7. integer literals (0 and 1 being the most interesting)
#
# it does not support looking at the value that is stored in
# any preprocessor variable or doing any math. If any of the
# following characters are present in an expression, then
# the parser/tokenizer will bail:
#    +, -, *, /, ==, !=, <=, >=, <, >, ,,
# the parser will also bail if macros are encountered --
# the only time a variable name can appear is if it is following
# a "defined" (or inside a parenthesis-pair following a defined).
# The entirety of the preprocessor logic must also live on a
# single line


class PreProcTokenTypes(Enum):
     NOT = "not"
     LPAREN = "lparen"
     RPAREN = "rparen"
     DEFINED = "defined"
     AND = "and"
     OR = "or"
     VAR = "var"
     LITERAL = "literal"
     UNPROCESSABLE = "unprocessable"


class PreProcToken:
    def __init__(self):
        self.spelling = ""
        self.start = -1
        self.one_past_end = -1
        self.type = None


class PreProcScanner:
    def __init__(self):
        self.whitespace = set([" ", "\t"])
        self.dividers = set(["(", ")", "[", "]", "=", "&", "|", "\\", '"', "'", ",", "%", ">", "<"])
        self.two_character_dividers = set(["&&", "||", "<=", ">=", "==", "!="])
        self.unprocessable_tokens = set(["==", "!=", "<=", "<", ">=", ">", "*", "+", "/", ",", "%", "[", "]", "'", '"', "=", "&", "|", "-"])
        self.literal_starts = set(["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "+", "-"])
                                         
    def tokenize_line(self, line):
        """The #if must have already been removed from the line"""

        tokens = []
        start = -1
        i = 0
        while i < len(line):
            if line[i] in self.whitespace:
                if start >= 0:
                    tokens.append(self.take_token_for_range(line, start, i))
                start = -1
                i += 1
            elif line[i] in self.dividers:
                if start >= 0:
                    tokens.append(self.take_token_for_range(line, start, i))
                start = -1
                if i+1 < len(line) and line[i:i+2] in self.two_character_dividers:
                    tokens.append(self.take_token_for_range(line, i,i+2))
                    i += 2
                else:
                    tokens.append(self.take_token_for_range(line, i,i+1))
                    i += 1
            else:
                if start < 0:
                    start = i
                i += 1
        if start >= 0:
            tokens.append(self.take_token_for_range(line, start, len(line)))

        return tokens
    def take_token_for_range(self, line, start, one_past_end):
        tok = PreProcToken()
        tok.spelling = line[start:one_past_end]
        tok.start = start
        tok.one_past_end = one_past_end
        #print("token", tok.spelling, tok.start, tok.one_past_end)

        if tok.spelling == "defined":
            tok.type = PreProcTokenTypes.DEFINED
        elif tok.spelling == "!":
            tok.type = PreProcTokenTypes.NOT
        elif tok.spelling == "(":
            tok.type = PreProcTokenTypes.LPAREN
        elif tok.spelling == ")":
            tok.type = PreProcTokenTypes.RPAREN
        elif tok.spelling == "&&":
            tok.type = PreProcTokenTypes.AND
        elif tok.spelling == "||":
            tok.type = PreProcTokenTypes.OR
        elif tok.spelling in self.unprocessable_tokens:
            tok.type = PreProcTokenTypes.UNPROCESSABLE
        else:
            # we either have a variable, or a literal
            tok.type = PreProcTokenTypes.VAR if tok.spelling[0] not in self.literal_starts else PreProcTokenTypes.LITERAL
        return tok

class PreProcASTAtom:
    # variable or literal
    def __init__(self):
        self.type = None
        self.value = None
        self.varname = None

    def eval(self):
        return bool(self.value)

class PreProcASTUnary:
    # defined, not, or parens
    def __init__(self):
        self.type = None
        self.sub_expr = None
        self.defined_set = None

    def eval(self):
        if self.type == PreProcTokenTypes.DEFINED:
            return self.sub_expr.varname in self.defined_set
        else:
            return not self.sub_expr.eval()
            

class PreProcASTFactor:
    # "&&"
    def __init__(self):
        self.type = None
        self.left = None
        self.right = None
    def eval(self):
        assert self.type == PreProcTokenTypes.AND
        assert self.left
        assert self.right
        return self.left.eval() and self.right.eval()

class PreProcASTTerm:
    # "||"
    def __init__(self):
        self.type = None
        self.left = None
        self.right = None
    def eval(self):
        assert self.type == PreProcTokenTypes.OR
        assert self.left
        assert self.right
        return self.left.eval() or self.right.eval()


class PreProcParser:
    def __init__(self):
        self.scanner = PreProcScanner()
        self.tokens = []
        self.defined_nodes = []
    def tree_from_line(self, line):
        self.defined_nodes = []
        self.tokens = self.scanner.tokenize_line(line)
        # for i,tok in enumerate(self.tokens):
        #      print(i, tok.spelling)
        if self.tokens_contain_unprocessable():
            return None
        root, ind = self.parse_expression(0)
        # print("tree from lines:", root, ind, len(self.tokens))
        if ind < len(self.tokens):
            return None
        else:
            return root

    def tokens_contain_unprocessable(self):
        for tok in self.tokens:
            if tok.type == PreProcTokenTypes.UNPROCESSABLE:
                return True
        return False
        

    def parse_expression(self, tok_ind):
        #print("parse expression", tok_ind)
        term, ind = self.parse_term(tok_ind)
        if term is None:
            return None, ind
        if ind < len(self.tokens):
            #print("next token type?", ind, self.tokens[ind].type)
            if self.tokens[ind].type == PreProcTokenTypes.OR:
                term2, ind2 = self.parse_expression(ind+1)
                if term2 is None:
                    return None, ind2
                root = PreProcASTTerm()
                root.type = PreProcTokenTypes.OR
                root.left = term
                root.right = term2
                return root, ind2
            
        return term, ind

    def parse_term(self, tok_ind):
        #print("parse term", tok_ind)
        fact, ind = self.parse_factor(tok_ind)
        if fact is None:
            return None, ind
        if ind < len(self.tokens):
            if self.tokens[ind].type == PreProcTokenTypes.AND:
                fact2, ind2 = self.parse_term(ind+1)
                if fact2 is None:
                    return None, ind2
                root = PreProcASTFactor()
                root.type = PreProcTokenTypes.AND
                root.left = fact
                root.right = fact2
                return root, ind2
        return fact, ind

    def parse_factor(self, tok_ind):
        #print("parse factor", tok_ind, self.tokens[tok_ind].type)
        if self.tokens[tok_ind].type == PreProcTokenTypes.LITERAL:
            root = PreProcASTAtom()
            root.type = PreProcTokenTypes.LITERAL
            try:
                root.value = float(self.tokens[tok_ind].spelling)
            except:
                return None, tok_ind
            return root, tok_ind+1
        elif self.tokens[tok_ind].type == PreProcTokenTypes.VAR:
            root = PreProcASTAtom()
            root.type = PreProcTokenTypes.VAR
            root.varname = self.tokens[tok_ind].spelling
            return root, tok_ind+1
        elif self.tokens[tok_ind].type == PreProcTokenTypes.LPAREN:
            root, ind = self.parse_expression(tok_ind+1)
            #print("rparen next?", ind, self.tokens[ind].type)
            if self.tokens[ind].type == PreProcTokenTypes.RPAREN:
                return root, ind+1
            else:
                return None, ind
        elif self.tokens[tok_ind].type == PreProcTokenTypes.NOT or self.tokens[tok_ind].type == PreProcTokenTypes.DEFINED:
            factor, ind = self.parse_factor(tok_ind+1)
            if factor is None:
                return factor, ind
            if self.tokens[tok_ind].type == PreProcTokenTypes.DEFINED:
                if factor.type != PreProcTokenTypes.VAR:
                    return None, ind
            root = PreProcASTUnary()
            root.type = self.tokens[tok_ind].type
            root.sub_expr = factor
            if self.tokens[tok_ind].type == PreProcTokenTypes.DEFINED:
                self.defined_nodes.append(root)
            return root, ind
        else:
            return None, tok_ind
    

            
