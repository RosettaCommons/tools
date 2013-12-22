'''

    blargs command line parser
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Easy argument parsing with powerful dependency support.
    :copyright: (c) 2012 by Karl Gyllstrom
    :license: BSD (see LICENSE.txt)

'''


from __future__ import print_function

import os
import operator
from functools import partial, wraps
from itertools import starmap, permutations
import sys


if sys.version_info[0] == 3:
    iterkeys = lambda x: x.keys()
    iteritems = lambda x: iter(x.items())
    isstring = lambda x: isinstance(x, str)
    from urllib.parse import urlparse
    xrange = range
    import configparser as cpars
else:
    iterkeys = lambda x: x.iterkeys()
    iteritems = lambda x: x.iteritems()
    isstring = lambda x: isinstance(x, basestring)
    from urlparse import urlparse
    import ConfigParser as cpars


class Multidict(object):
    def __init__(self, dictionary=None):
        self._values = {}
        if dictionary is not None:
            # XXX should we make sure dictionary conforms to our schema?
            self._values.update(dictionary)

    def __delitem__(self, key):
        return self._values.__delitem__(key)

    def __contains__(self, key):
        return self._values.__contains__(key)

    def copy(self):
        return Multidict(self._values.copy())

    def get(self, key):
        return self._values.get(key)

    def __str__(self):
        return self._values.__str__()

    def overwrite(self, key, value):
        self._values[key] = value

    def __setitem__(self, key, value):
        if key in self._values:
            v = self._values[key]
            if not isinstance(v, list):
                v = [v]

            value = v + [value]

        self._values[key] = value

    def __getitem__(self, key):
        return self._values.__getitem__(key)

    def __iter__(self):
        return iteritems(self._values)


class _ConfigCaster(object):
    def __init__(self, parent):
        self._parent = parent

    def __call__(self, filename):
        cfp = cpars.ConfigParser()
        cfp.readfp(open(filename))
        for sec in cfp.sections():
            for key, value in cfp.items(sec):
                yield key, value


class _RangeCaster(object):
    def __call__(self, value):
        def raise_error():
            raise FormatError(('%s is not range format: N:N+i, N-N+i, or N'
                    + ' N+i') % value)

        splitter = None
        for char in (' :-'):
            if char in value:
                splitter = char
                break

        if splitter:
            toks = value.split(splitter)
        else:
            toks = [value]

        if not (1 <= len(toks) <= 3):
            raise_error()
        try:
            return xrange(*[int(y) for y in toks])
        except ValueError:
            raise_error()


class _DirectoryOpenerCaster(object):
    def __init__(self, create):
        self._create = create

    def __call__(self, name):
        if not os.path.exists(name):
            if self._create:
                os.makedirs(name)
                return name
            raise IOError('%s does not exist' % name)

        if not os.path.isdir(name):
            raise IOError('%s is not directory' % name)

        return name


class _FileOpenerCaster(object):
    def __init__(self, mode=None, buffering=None):
        self._kw = {}
        if mode is not None:
            self._kw['mode'] = mode
        if buffering is not None:
            self._kw['buffering'] = buffering

    def __call__(self, *args):
        return open(*args, **self._kw)


# ---------- decorators ---------- #


def _names_to_options(f):
    @wraps(f)
    def inner(*args, **kwargs):
        new_args = []
        for arg in args[1:]:
            if isstring(arg):
                arg = args[0]._options[arg]

            new_args.append(arg)

        return f(*[args[0]] + new_args, **kwargs)

    return inner


def _options_to_names(f):
    ''' Convert any :class:`Option`s to names. '''

    @wraps(f)
    def inner(*args, **kwargs):
        new_args = []
        for arg in args[1:]:
            if isinstance(arg, Option) and not isinstance(arg, Group):
                arg = arg.argname

            new_args.append(arg)

        return f(*[args[0]] + new_args, **kwargs)

    return inner


def _verify_args_exist(f):
    @wraps(f)
    def inner(*args, **kwargs):
        def raise_error(name):
            raise ValueError('%s not known' % name)

        self = args[0]
        for arg in args[1:]:
            if isstring(arg) and arg not in self._readers:
                raise_error(arg)

        return f(*args, **kwargs)

    return inner


def _localize_all(f):
    @wraps(f)
    def inner(*args, **kwargs):
        args = list(args)
        self = args[0]
        new_args = []
        for arg in args[1:]:
            if isstring(arg):
                arg = self._localize(arg)
            new_args.append(arg)
        args[1:] = new_args
        return f(*args, **kwargs)

    return inner


def localize(f):
    @wraps(f)
    def inner(*args, **kwargs):
        args = list(args)
        self = args[0]
#        args[1:] = [self._localize(x) for x in args[1]]
        args[1] = self._localize(args[1])

        return f(*args, **kwargs)

    return inner

# ---------- end decorators ---------- #


# ---------- exceptions ---------- #


class ArgumentError(ValueError):
    ''' Root class of all arguments that are thrown due to user input which
    violates a rule set by the parser. In other words, errors of this type
    should be caught and communicated to the user some how. The default
    behavior is to signal the particular error and show the `usage`.'''

    pass


class FormatError(ArgumentError):
    ''' Argument not formatted correctly. '''
    pass


class ConditionError(ArgumentError):
    ''' Condition not met. '''

    def __init__(self, argname, condition):
        super(ConditionError, self).__init__()
        self.argname = argname
        self.condition = condition

    def __str__(self):
        return '%s required unless %s' % (self.argname, self.condition)


class MissingRequiredArgumentError(ArgumentError):
    ''' Required argument not specified. '''

    def __init__(self, arg):
        if isinstance(arg, Option):
            arg = arg.argname

        super(MissingRequiredArgumentError, self).__init__(
                'No value passed for %s' % arg)


class ManyAllowedNoneSpecifiedArgumentError(ArgumentError):
    ''' An argument out of a list of required is missing.
    e.g., one of -a and -b is required and neither is specified. '''

    def __init__(self, allowed):
        super(ManyAllowedNoneSpecifiedArgumentError, self).__init__(('[%s] not'
                + ' specified') % ', '.join(sorted(map(str, allowed))))


class UnspecifiedArgumentError(ArgumentError):
    ''' User supplies argument that isn't specified. '''

    message = 'illegal option {0}'

    def __init__(self, arg):
        super(UnspecifiedArgumentError,
                self).__init__(UnspecifiedArgumentError.message.format(arg))


class MultipleSpecifiedArgumentError(ArgumentError):
    ''' Multiple of the same argument specified. '''
    pass


class DependencyError(ArgumentError):
    ''' User specified an argument that requires another, unspecified argument.
    '''

    def __init__(self, arg1, arg2):
        super(DependencyError, self).__init__('%s requires %s' % (arg1, arg2))


class ConflictError(ArgumentError):
    ''' User specified an argument that conflicts another specified argument.
    '''

    def __init__(self, offender1, offender2):
        super(ConflictError, self).__init__('%s conflicts with %s' %
                (offender1, offender2))


class FailedConditionError(ArgumentError):
    ''' Condition failed. '''
    pass


class MissingValueError(ArgumentError):
    ''' Argument flag specified without value. '''
    pass


class InvalidEnumValueError(ArgumentError):
    ''' Enum value provided not allowed. '''
    pass

# ---------- end exceptions ---------- #


class Condition(object):
    def __init__(self):
        self._other_conditions = []
        self._and = True
        self._neg = False

    def _copy(self):
        c = self.__new__(self.__class__)
        c._other_conditions = self._other_conditions[:]
        c._and = self._and
        c._neg = self._neg
        return c

    def __neg__(self):
        c = self._copy()
        c._neg = True
        return c

    def _inner_satisfied(self, parsed):
        raise NotImplementedError

    def and_(self, condition):
        if not self._and:
            raise ValueError('and/or both specified')

        c = self._copy()
        c._other_conditions.append(condition)
        return c

    def or_(self, condition):
        c = self._copy()
        c._other_conditions.append(condition)
        c._and = False
        return c

    def _is_satisfied(self, parsed):
        def _inner():
            for n in self._other_conditions:
                if not n._is_satisfied(parsed):
                    if self._and:
                        return False
                elif not self._and:
                    return True

            return self._inner_satisfied(parsed)

        result = _inner()
        if self._neg:
            result = not result
        return result


class _CallableCondition(Condition):
    def __init__(self, call, main, other):
        super(_CallableCondition, self).__init__()
        self._call = call
        self._main = main
        self._other = other

    def _copy(self):
        c = super(_CallableCondition, self)._copy()
        c._call = self._call
        c._main = self._main
        c._other = self._other
        return c

    def _inner_satisfied(self, parsed):
        operands = []
        for item in (self._main, self._other):
            if isinstance(item, Option):
                v = parsed.get(item.argname)
                if isinstance(v, list):
                    v = [vi.getvalue() for vi in v]
                else:
                    v = [v.getvalue()]
            else:
                v = [item]

            operands.append(v)

        for x1 in operands[0]:
            for x2 in operands[1]:
                if not self._call(x1, x2):
                    return False

        return True

    def __repr__(self):
        x = self._main.argname

        names = {operator.le: '<=',
                 operator.lt: '<',
                 operator.gt: '>',
                 operator.ge: '>=',
                 operator.eq: '==',
                 operator.ne: '!='}

        x += ' ' + names[self._call] + ' ' + str(self._other)
        return x


class Option(Condition):
    ''' Do not construct directly, as it will not be tethered to a
    :class:`Parser` object and so will not be handled in argument parsing. '''

    def __init__(self, argname, parser):
        super(Option, self).__init__()

        self.argname = argname
        self._parser = parser
        self._conditions = []
        self._allows_multiple = False
        self._description = None

    def _copy(self):
        c = super(Option, self)._copy()
        c.argname = self.argname
        c._parser = self._parser
        c._conditions = self._conditions
        return c

    def described_as(self, description):
        '''

        Add a description to this argument, to be displayed when users trigger the help message.

        :param description: description of argument

        '''

        self._description = description
        return self

    def requires(self, *conditions):
        '''
        
        Specify other options/conditions which this argument requires.

        :param conditions: required conditions
        :type others: sequence of either :class:`Option` or :class:`Condition`

        '''

        [self._parser._set_requires(self.argname, x) for x in conditions]
        return self

    def conflicts(self, *conditions):
        ''' Specify other conditions which this argument conflicts with.

        :param conditions: conflicting options/conditions
        :type conditions: sequence of either :class:`Option` or
            :class:`Condition`
        '''

        [self._parser._set_conflicts(self.argname, x) for x in conditions]
        return self

    def shorthand(self, alias):
        '''
        
        Set shorthand for this argument. Shorthand arguments are 1
        character in length and are prefixed by a single '-'. For example:

        ::

            parser.str('option').shorthand('o')

        would cause '--option' and '-o' to be alias argument labels when
        invoked on the command line.

        :param alias: alias of argument

        '''

        self._parser._add_shorthand(self.argname, alias)
        return self

    def default(self, value):
        ''' Provide a default value for this argument. '''
        self._parser._set_default(self.argname, value)
        return self

    def environment(self):
        ''' Pull argument value from OS environment if unspecified. The case of
        the argument name, all lower, and all upper are all tried. For example,
        if the argument name is `Port`, the following names will be used for
        environment lookups: `Port`, `port`, `PORT`.

        ::

            with Parser(locals()) as p:
                p.int('port').environment()

        Both command lines work:

        ::

            python test.py --port 5000
            export PORT=5000; python test.py

        '''

        default = os.environ.get(self.argname)
        if default is None:
            default = os.environ.get(self.argname.lower())

            if default is None:
                default = os.environ.get(self.argname.upper())

        if default is not None:
            # XXX not tested
            return self.default(default)

        return self

    def cast(self, cast):
        ''' Provide a casting value for this argument. '''

        self._parser._readers[self.argname] = Caster(
                self._parser._readers[self.argname], cast)

        return self

    def required(self):
        ''' Indicate that this argument is required. '''

        self._parser._set_required(self.argname)
        return self

    def if_(self, condition):
        ''' Argument is required if ``conditions``. '''
        self._parser._set_required(self.argname, [-condition])
        return self

    def unless(self, condition):
        ''' Argument is required unless ``conditions``. '''

        self._parser._set_required(self.argname, [condition])
        return self

    def unspecified_default(self):
        ''' Indicate that values passed without argument labels will be
        attributed to this argument. '''

        self._parser._set_unspecified_default(self.argname)
        return self

    def multiple(self):
        ''' Indicate that the argument can be specified multiple times. '''

        self._allows_multiple = True
        return self

    # --- conditions

    def _inner_satisfied(self, parsed):
        v = parsed.get(self.argname)

        if not isinstance(v, list):
            v = [v]

        return all(x.is_resolvable() for x in v)

    def _make_condition(self, func, other):
        return _CallableCondition(func, self, other)

    def __le__(self, other):
        return self._make_condition(operator.__le__, other)

    def __lt__(self, other):
        return self._make_condition(operator.__lt__, other)

    def __gt__(self, other):
        return self._make_condition(operator.__gt__, other)

    def __ge__(self, other):
        return self._make_condition(operator.__ge__, other)

    def __eq__(self, other):
        return self._make_condition(operator.__eq__, other)

    def __ne__(self, other):
        return self._make_condition(operator.__ne__, other)

    def __hash__(self):
        return hash(str(self.argname))

    # -- private access methods

    def _isrequired(self):
        return self in self._parser._required

    def _getconflicts(self):
        return self._parser._conflicts.get(self, [])

    def _getreqs(self):
       return self._parser._requires.get(self, [])

    def _alias(self):
        return self._parser._source_to_alias.get(self.argname)

    def __str__(self):
        v = '%s%s' % (self._parser._double_prefix, self.argname)
        alias = self._alias()
        if alias:
            v += '/%s%s' % (self._parser._single_prefix, alias)
        return v


# ---------- Argument readers ---------- #

# Argument readers help parse the command line. The flow is as follows:
#
# 1) We create a reader for each argument based on the type specified by the
#    programmer (e.g., a bool will get a _FlagArgumentReader
#
# 2) The reader is used when we hit the argument label


class _ArgumentReader(object):
    class _UNSPECIFIED(object):
        def __repr__(self):
            return 'UNSPECIFIED'

    UNSPECIFIED = _UNSPECIFIED()

    def __init__(self, parent):
        self.value = _ArgumentReader.UNSPECIFIED
        self.parent = parent
        self._default = _ArgumentReader.UNSPECIFIED
        self._init()

    def fresh_copy(self):
        return self.__class__(self.parent)

    def _init(self):
        pass

    def activate(self):
        pass

    def _set_default(self, default):
        self._default = default

    def consume_or_skip(self, arg):
        raise NotImplementedError

    def is_specified(self):
        return self.value is not _ArgumentReader.UNSPECIFIED

    def is_resolvable(self):
        return (self.is_specified() or self._default is not
                _ArgumentReader.UNSPECIFIED)

    def default(self):
        if self._default is not _ArgumentReader.UNSPECIFIED:
            return self._default
        return self.__class__.class_default()

    @classmethod
    def class_default(cls):
        return _ArgumentReader.UNSPECIFIED

    def _get(self):
        raise NotImplementedError

    def getvalue(self):
        if self.is_specified():
            return self._get()

        return self.default()


class _MultiWordArgumentReader(_ArgumentReader):
    def consume_or_skip(self, arg):
        if self.parent._is_argument_label(arg):
            return False

        if not self.is_specified():
            self.value = []

        self.value.append(arg)
        return True

    def _get(self):
        if not self.is_specified() or len(self.value) == 0:
            # XXX
            raise MissingValueError

        return ' '.join(self.value)


class _FlagArgumentReader(_ArgumentReader):
    def _init(self):
        self.value = False

    def activate(self):
        self.value = True

    def is_resolvable(self):
        return self.value

    def consume_or_skip(self, arg):
        return False

    def _get(self):
        return self.value

    @classmethod
    def class_default(cls):
        return False


class _SingleWordReader(_ArgumentReader):
    def consume_or_skip(self, arg):
        if self.is_specified():
            return False

        self.value = arg
        return True

    def _get(self):
        if not self.is_specified():
            raise MissingValueError
        return self.value


class Caster(object):
    def __init__(self, reader, cast):
        self._reader = reader
        self._cast = cast

    def fresh_copy(self):
        # XXX does _cast need to be copied too?
        return self.__class__(self._reader.fresh_copy(), self._cast)

    def activate(self):
        self._reader.activate()

    def getvalue(self):
        try:
            v = self._reader.getvalue()
            if v is _ArgumentReader.UNSPECIFIED:
                return None

            return self._cast(v)
        except ValueError:
            raise FormatError

    def is_resolvable(self):
        return self._reader.is_resolvable()

    def consume_or_skip(self, arg):
        return self._reader.consume_or_skip(arg)

    def is_specified(self):
        return self._reader.is_specified()

    def _set_default(self, default):
        self._reader._set_default(default)


# ---------- Argument readers ---------- #


class Group(Option):
    def __init__(self, parser, *names):
        self._parser = parser
#        super(Group, self).__init__('group', parser)
        self._names = names
        self.argname = self

    def default(self, name):
        if name not in self._names:
            raise ValueError('%s not in group' % name)

        self._default = name
        return self

    def _is_satisfied(self, parsed):
        for item in self._names:
            if item._is_satisfied(parsed):
                return True
        return False

    def __str__(self):
        return ', '.join([str(name) for name in self._names])
            

class Parser(object):
    ''' Command line parser. '''

    def __init__(self, store=None):
        self._options = {}
        self._readers = {}
        self._extras = []
        self._unspecified_default = None

        # dict of A (required) -> args that could replace A if A is not
        # specified
        self._required = {}

        # dict of A -> args that A depends on
        self._requires = {}

        # dict of A -> args that A conflicts with
        self._conflicts = {}

        self._alias = {}
        self._source_to_alias = {}

        # convert arg labels to understore in value dict
        self._to_underscore = False

        # prefix to shorthand args
        self._single_prefix = '-'

        # prefix to fullname args
        self._double_prefix = '--'

        # help message
        self._help_prefix = None

        self._sys_exit_error = SystemExit
        self.out = sys.stdout # XXX not documented

        # set by user
        self._init_user_set(store)

        self.flag('help').shorthand('h').described_as('Print help message.')

    def _init_user_set(self, store=None):
        if store is None:
            store = {}
        self._store = store
        self._namemaps = {}
        self._rnamemaps = {}

    def _argument_exists(self, name_or_alias):
        return name_or_alias in self._readers or name_or_alias in self._alias

    def set_help_prefix(self, message):
        ''' Indicate text to appear before argument list when the ``help``
        function is triggered. '''

        self._help_prefix = message
        return self

    def underscore(self):
        ''' Convert '-' to '_' in argument names. This is enabled if
        ``with_locals`` is used, as variable naming rules are applied. '''

        self._to_underscore = True
        return self

    def set_single_prefix(self, flag):
        ''' Set the single flag prefix. This appears before short arguments
        (e.g., -a). '''

        if self._double_prefix in flag:
            raise ValueError(
                    'single_prefix cannot be superset of double_prefix')

        self._single_prefix = flag
        return self

    def set_double_prefix(self, flag):
        ''' Set the double flag prefix. This appears before long arguments
        (e.g., --arg). '''

        if flag in self._single_prefix:
            raise ValueError('single_flag cannot be superset of double_flag')

        self._double_prefix = flag
        return self

    def use_aliases(self):
        raise NotImplementedError

    @classmethod
    def with_locals(cls):
        ''' Create :class:`Parser` using locals() dict. '''

        import inspect
        vals = inspect.currentframe().f_back.f_locals
        return Parser(vals).underscore()

# --- types --- #

    def config(self, name):
        '''
        
        Add configuration file, whose key/value pairs will provide/replace any
        arguments created for this parser. Configuration format should be INI.
        For example:

            ::

                with Parser() as p:
                    p.int('a')
                    p.str('b')
                    p.config('conf')

            Now, arg ``a`` can be specfied on the command line, or in the
            configuration file passed to ``conf``. For example:

            ::

                python test.py --a 3
                python test.py --conf myconfig.cfg

            Where myconfig.cfg:

            ::

                [myconfigfile]
                a = 5
                b = 9
                x = 'hello'

            Note that any parameters in the config that aren't created as
            arguments via this parser are ignored. In the example above, the
            values of variables ``a`` and ``b`` would be assigned, while ``x``
            would be ignored (as the developer did not create an ``x``
            argument).

            If an argument is specified both on the command line and within the
            config file, it must be indicated as allowing multiple (via
            :py:meth:`.multiple`).

            '''

        return self.str(name).cast(_ConfigCaster(self))

    def enum(self, name, values):
        arg = self.str(name)
        cond = arg == values[0]
        for v in values[1:]:
            cond = cond.or_(arg == v)

        return arg.requires(cond)

    def int(self, name):
        ''' Add integer argument. '''
        return self._add_option(name).cast(int)

    def float(self, name):
        ''' Add float argument. '''

        return self._add_option(name).cast(float)

#    def enum(self, name, values):
#        ''' Add enum type. '''
#
#        def inner(value):
#            if not value in values:
#                raise InvalidEnumValueError()
#            return value
#
#        return self._add_option(name).cast(inner)

    def str(self, name):
        ''' Add :py:class:`str` argument. '''
        return self._add_option(name)

    def range(self, name):
        ''' Range type. Accepts similar values to that of python's py:`range`
            and py:`xrange`. Accepted delimiters are space, -, and :.

            ::

                with Parser() as p:
                    p.range('values')

            Now accepts:

            ::

              python test.py --values 10  # -> xrange(10)
              python test.py --values 0-1  # -> xrange(0, 1)
              python test.py --values 0:10:2  # -> xrange(0, 10, 2)
              python test.py --values 0 10 3  # -> xrange(0, 10, 3)
        '''

        return self.multiword(name).cast(_RangeCaster())

    def multiword(self, name):
        ''' Accepts multiple terms as an argument. For example:

            ::

                with Parser() as p:
                    p.multiword('multi')

            Now accepts:

            ::

               python test.py --multi path to something

        '''

        result = self._add_option(name)
        self._set_reader(name, _MultiWordArgumentReader(self))
        return result

    def bool(self, name):
        ''' Alias of :func:`flag`. '''

        return self.flag(name)

    def flag(self, name):
        ''' Boolean value. The presence of this flag indicates a true value,
        while an absence indicates false. No arguments. '''

        result = self._add_option(name)
        self._set_reader(name, _FlagArgumentReader(self))
        return result

    def file(self, name, mode=None, buffering=None):
        ''' Opens the file indicated by the name passed by the user. ``mode``
        and ``buffering`` are arguments passed to ``open``.

        The example below implements a file copy operation:

        ::

            with Parser(locals()) as p:
                p.file('input_file')
                p.file('output_file', mode='w')
                
            output_file.write(input_file.read())

        '''

        return self.multiword(name).cast(_FileOpenerCaster(mode, buffering))

    def directory(self, name, create=False):
        ''' File directory value. Checks to ensure that the user passed file
        name exists and is a directory (i.e., not some other file object). If
        ``create`` is specified, creates the directory using ``os.makedirs``;
        any intermediate directories are also created. '''

        return self.multiword(name).cast(_DirectoryOpenerCaster(create))

    def url(self, name):
        ''' URL value; verifies that argument has a scheme (e.g., http, ftp,
        file). '''

        def parse(value):
            if urlparse(value).scheme == '':
                raise FormatError('%s not valid URL' % value)
            return value

        return self.str(name).cast(parse)

# --- aggregate calls --- #

    def at_least_one(self, *args):
        ''' Require at least one of ``args``. '''

        return self._require_at_least_one(*args)

    def require_one(self, *args):
        ''' Require only and only one of ``args``. '''

        self._set_one_required(*args)
        return Group(self, *args)

    def all_if_any(self, *args):
        ''' If *any* of ``args`` is specified, then all of ``args`` must be
        specified. '''

        list(starmap(self._set_requires, permutations(args, 2)))
        return Group(self, *args)

    def only_one_if_any(self, *args):
        ''' If *any* of ``args`` is specified, then none of the remaining
        ``args`` may be specified.'''

        list(starmap(self._set_conflicts, permutations(args, 2)))
        return Group(self, *args)

    def __getitem__(self, name):
        return Option(name, self)

# --- private --- #

    def _require_at_least_one(self, *names):
        ''' At least one of the arguments is required. '''

        s = set(names)
        for v in s:
            self._set_required(v, s - set([v]))

        return Group(self, *names)

    @localize
    def _set_unspecified_default(self, name):
        if self._unspecified_default is not None:
            raise ValueError('Trying to specify multiple unspecified defaults')

        self._unspecified_default = name

    @localize
    @_verify_args_exist
    @_names_to_options
    def _set_required(self, arg, replacements=None):
        newreplacements = []
        for item in replacements or []:
            if not isinstance(item, (Option, Condition)):
                item = self._options[item]
            newreplacements.append(item)

        # XXX [] allows redundancy?
        self._required.setdefault(arg, []).extend(newreplacements)

    @localize
    def _set_reader(self, name, option):
        self._readers[name] = option

    def _add_shorthand(self, source, alias):
        if source not in self._readers:
            raise ValueError('%s not an option' % source)
        if alias in self._alias:
            raise ValueError('{0} already shorthand for {1}'.format(alias,
                self._alias[alias]))

        self._alias[alias] = source
        self._source_to_alias[source] = alias

    def _add_option(self, name):
        name = self._localize(name)

        if name in self._readers:
            raise ValueError('multiple types specified for %s' % name)

        self._set_reader(name, _SingleWordReader(self))

        o = Option(name, self)
        self._options[name] = o
        return o

    def _getoption(self, option):
        o = self._readers.get(option, None)
        if o is not None:
            return option, o

        source = self._alias.get(option, None)
        if source is None:
            return option, None

        return source, self._readers.get(source, None)

    def _localize(self, key):
        if self._to_underscore:
            modified = key.replace('-', '_')
            self._rnamemaps.setdefault(modified, key)
            return self._namemaps.setdefault(key, modified)
        return key

    def _unlocalize(self, key):
        v = self._rnamemaps.get(key, None)
        if v is None:
            # already unlocalized(?)
            return key
        return v

    def _tokenize(self, args):
        new_args = []
        for arg in args:
            if '=' in arg:
                new_args += arg.split('=')
            else:
                new_args.append(arg)

        return new_args

    def _is_argument_label(self, arg):
        return (arg.startswith(self._single_prefix) or
                arg.startswith(self._double_prefix))

    def _parse(self, tokenized):
        current_reader = None
        parsed = Multidict()

        for arg in tokenized:
            if current_reader is not None:
                if current_reader.consume_or_skip(arg):
                    continue
                current_reader = None

            argument_name = None

            if self._is_argument_label(arg):
                is_full = arg.startswith(self._double_prefix)

                if is_full:
                    prefix = self._double_prefix
                else:
                    prefix = self._single_prefix

                arg = arg[len(prefix):]

                argument_name = self._localize(arg)

                if is_full:
                    argument_name = self._options.get(argument_name)
                else:
                    argument_name = self._options.get(
                            self._alias.get(argument_name))

                if argument_name is not None:
                    argument_name = argument_name.argname

                current_reader = self._readers.get(argument_name)

                if current_reader is None:
                    if argument_name is None:
                        argument_name = arg
                    raise UnspecifiedArgumentError(argument_name)

                current_reader = current_reader.fresh_copy()
                current_reader.activate()

            elif self._unspecified_default is not None:
                argument_name = self._unspecified_default

                # push value onto _SingleWordReader
                current_reader = _SingleWordReader(self)
                current_reader.consume_or_skip(arg)

            if argument_name:
                parsed[argument_name] = current_reader
            else:
                self._extras.append(arg)

        for k, v in parsed:
            if not isinstance(v, list):
                v = [v]
            for item in v:
                if not item.is_specified():
                    raise MissingValueError

        return parsed

    def _help_if_necessary(self, processed):
        if 'help' in processed:
            self.print_help()
            raise self._sys_exit_error(0)

    def _config_values(self, parsed):
        pc = parsed.copy()
        for key, value in parsed:
            if isinstance(value, Caster) and isinstance(value._cast,
                    _ConfigCaster):  # XXX ugly

                for k, v in value.getvalue():
                    current_reader = pc.get(k)

                    if current_reader is None:
                        # developer didn't specify this argument
                        continue

                    is_specified = current_reader.is_specified()
                    if is_specified:
                        current_reader = current_reader.fresh_copy()

                    current_reader.consume_or_skip(v)

                    if is_specified:
                        pc[k] = current_reader

                del pc[key]

        return pc

    def _assign(self, combined):
        assigned = {}
        for key, values in combined:
            try:
                if not self._options[key]._allows_multiple:
                    value = values.getvalue()
                else:
                    if not isinstance(values, list):
                        values = [values]

                    value = [v.getvalue() for v in values]

                if value is _ArgumentReader.UNSPECIFIED:
                    value = None

                assigned[key] = value

            except MissingValueError:
                raise MissingValueError('%s specified but missing given value'
                        % key)

        return assigned

    def _check_multiple(self, assigned):
        for key, values in assigned:
            if isinstance(values, list) and not self._options[key]._allows_multiple:
                raise MultipleSpecifiedArgumentError(('%s specified multiple' +
                    ' times') % self._options[key])

    def _check_required(self, assigned):
        for arg, replacements in iteritems(self._required):
            missing = []
            if not arg._is_satisfied(assigned):
                for v in replacements:
                    if isinstance(v, Group):
                        missing += v._names
                    elif not isinstance(v, Condition):
                        missing.append(v)

                    if v._is_satisfied(assigned):
                        break
                else:
                    if missing:
                        raise ManyAllowedNoneSpecifiedArgumentError([arg] +
                                missing)
                    else:
                        raise MissingRequiredArgumentError(arg)

    def _check_dependencies(self, assigned):
        for arg, deps in iteritems(self._requires):
            if arg._is_satisfied(assigned):
                for v in deps:
                    if not v._is_satisfied(assigned):
                        if isinstance(v, _CallableCondition):
                            raise ConditionError(arg.argname, v)
                        raise DependencyError(arg, v)

    def _check_conflicts(self, assigned):
        for arg, conflicts in iteritems(self._conflicts):
            if arg._is_satisfied(assigned):
                for conflict in conflicts:
                    if conflict._is_satisfied(assigned):
                        raise ConflictError(arg.argname, conflict.argname)

    def _verify(self, assigned):
        self._check_required(assigned)
        self._check_dependencies(assigned)
        self._check_conflicts(assigned)

    def _assign_to_store(self, assigned):
        for key, value in iteritems(assigned):
            self._store[key] = value

    @_options_to_names
    def _set_one_required(self, *names):
        self.only_one_if_any(*names)
        self._require_at_least_one(*names)

    def _get_args(self, args):
        if args is None:
            args = sys.argv[1:]

        if not isinstance(args, list):
            raise TypeError('%s not list' % args)

        return args

    def _combine_with_defaults(self, user_args):
        copy = Multidict(self._readers.copy())

        for k, v in user_args:
            copy.overwrite(k, v)

        return copy

    def _process_command_line(self, args=None):
        try:
            args = self._get_args(args)
            tokenized = self._tokenize(args)
            user_args = self._parse(tokenized)
            self._help_if_necessary(user_args)
            user_args = self._combine_with_defaults(user_args)
            user_args = self._config_values(user_args)
            self._check_multiple(user_args)
            self._verify(user_args)

            assigned = self._assign(user_args)
            self._assign_to_store(assigned)
        except ArgumentError as e:
            raise e

#        self._init_user_set()  # reset

        return self._store

    def _emit(self, *args):
        print(*args, file=self.out)

    def process_command_line(self, args=None):
        '''

        Parses the command line specified by the user. Consider this a manual alternative to the ``with`` idiom. The following examples are equivalent:

        ::

            with Parser(locals()) as p:
                p.int('var')

            print var

        ::

            p = Parser(locals())
            p.int('var')
            p.process_command_line()
            print var

        '''

        try:
            return self._process_command_line(args)
        except ArgumentError as e:
            self.bail(e)

    def _usage(self):
        return ('Usage: %s ' % sys.argv[0]) + ' '.join('[%s]' %
                self._label(value) for value in self._options.values())

    def _print_table(self, t):
        # XXX what about empty list?
        column_max_lengths = [len(x) for x in t[0]]

        for row in t[1:]:
            column_max_lengths = [max(column_max_lengths[i], len(column)) for i, column
                    in enumerate(row)]

        for row in t:
            line = ''
            for i, column in enumerate(row):
                fmt = '   %-' + str(column_max_lengths[i]) + 's'
                line += fmt

            self._emit(line % row)

    def bail(self, e):
        msg = []
        msg.append('Error: ' + str(e))
        msg.append(self._usage())

        self._emit('\n'.join(msg))

        raise self._sys_exit_error(1)

    @localize
    @_options_to_names
    def _set_default(self, name, value):
        self._readers[name]._set_default(value)
    #    self._defaults[name] = value

    @_localize_all
    @_verify_args_exist
    @_names_to_options
    def _set_requires(self, a, b):
        self._requires.setdefault(a, set()).add(b)

    @_localize_all
    @_verify_args_exist
    @_names_to_options
    def _set_conflicts(self, a, b):
        self._conflicts.setdefault(a, set()).add(b)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is None or isinstance(exc_type, ArgumentError):
            self.process_command_line()
            return True
        return False

    def _option_label(self, reader):
        if not isinstance(reader, Caster):
            return 'option'

        if reader._cast is int:
            return 'int'

        if reader._cast is float:
            return 'float'

        if reader._cast is float:
            return 'float'

        if isinstance(reader._cast, _RangeCaster):
            return 'range'

        return 'option'

    def _label(self, opt):
        pkey = str(opt)

        reader = self._readers[opt.argname]
        if not isinstance(reader, _FlagArgumentReader):
            pkey = '%s <%s>' % (pkey, self._option_label(reader))

        return pkey

    def print_help(self):
        self._emit(self._usage())

        if self._help_prefix:
            self._emit(self._help_prefix)

        self._emit('Options: (! denotes required argument)')
        column_width = -1

        labels = []
        for key, opt in iteritems(self._options):
            name = self._label(opt)
            desc = ''
            if opt._description is not None:
                desc = opt._description

            if opt._isrequired():
                name = '!' + name

            conflict_str = ''
            conflicts = opt._getconflicts()
            if conflicts:
                conflict_str = 'Conflicts with %s' % ', '.join(str(item) for item in
                        conflicts)

            requirement_str = ''
            reqs = opt._getreqs()
            if reqs:
                requirement_str = 'Requires %s' % ', '.join(str(item) for item in reqs)

            labels.append((name, desc, conflict_str, requirement_str))

        self._print_table(labels)


__all__ = ['Parser']
__version__ = '0.2.29b'
