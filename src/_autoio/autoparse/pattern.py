"""
autoparse.pattern
*****************

Generate re patterns for text parsing.
"""
# pattern generators
from autoparse._pattern import escape
from autoparse._pattern import maybe
from autoparse._pattern import preceded_by
from autoparse._pattern import not_preceded_by
from autoparse._pattern import followed_by
from autoparse._pattern import not_followed_by
from autoparse._pattern import zero_or_more
from autoparse._pattern import one_or_more
from autoparse._pattern import one_of_these
from autoparse._pattern import capturing
from autoparse._pattern import named_capturing
from autoparse._pattern import series
from autoparse._more_patterns import block_pattern
from autoparse._more_patterns import lpadded
from autoparse._more_patterns import rpadded
from autoparse._more_patterns import padded
# pattern constants #
from autoparse._lib import STRING_START
from autoparse._lib import STRING_END
from autoparse._lib import LINE_START
from autoparse._lib import LINE_END
from autoparse._lib import WILDCARD
from autoparse._lib import WILDCARD2
from autoparse._lib import NEWLINE
from autoparse._lib import NONNEWLINE
from autoparse._lib import LINE_FILL
from autoparse._lib import LINE
from autoparse._lib import SPACE
from autoparse._lib import SPACES
from autoparse._lib import ZSPACES
from autoparse._lib import LINESPACE
from autoparse._lib import LINESPACES
from autoparse._lib import PADDING
from autoparse._lib import NONSPACE
from autoparse._lib import UPPERCASE_LETTER
from autoparse._lib import LOWERCASE_LETTER
from autoparse._lib import PLUS
from autoparse._lib import MINUS
from autoparse._lib import PERIOD
from autoparse._lib import UNDERSCORE
from autoparse._lib import LETTER
from autoparse._lib import DIGIT
from autoparse._lib import URLSAFE_CHAR
from autoparse._lib import CKINSAFE_CHAR
from autoparse._lib import SIGN
from autoparse._lib import UNSIGNED_INTEGER
from autoparse._lib import UNSIGNED_FLOAT
from autoparse._lib import INTEGER
from autoparse._lib import FLOAT
from autoparse._lib import EXPONENTIAL_INTEGER
from autoparse._lib import EXPONENTIAL_FLOAT
from autoparse._lib import NUMBER
from autoparse._lib import EXPONENTIAL_INTEGER_D
from autoparse._lib import EXPONENTIAL_FLOAT_D
from autoparse._lib import VARIABLE_STRING
from autoparse._lib import VARIABLE_NAME

__all__ = [
    # pattern generators
    'escape',
    'maybe',
    'preceded_by',
    'not_preceded_by',
    'followed_by',
    'not_followed_by',
    'zero_or_more',
    'one_or_more',
    'one_of_these',
    'capturing',
    'named_capturing',
    'series',
    'block_pattern',
    'lpadded',
    'rpadded',
    'padded',
    # pattern constants
    'STRING_START',
    'STRING_END',
    'LINE_START',
    'LINE_END',
    'WILDCARD',
    'WILDCARD2',
    'NEWLINE',
    'NONNEWLINE',
    'LINE_FILL',
    'LINE',
    'SPACE',
    'SPACES',
    'ZSPACES',
    'LINESPACE',
    'LINESPACES',
    'PADDING',
    'NONSPACE',
    'UPPERCASE_LETTER',
    'LOWERCASE_LETTER',
    'PLUS',
    'MINUS',
    'PERIOD',
    'UNDERSCORE',
    'LETTER',
    'DIGIT',
    'URLSAFE_CHAR',
    'CKINSAFE_CHAR',
    'SIGN',
    'UNSIGNED_INTEGER',
    'UNSIGNED_FLOAT',
    'INTEGER',
    'FLOAT',
    'EXPONENTIAL_INTEGER',
    'EXPONENTIAL_FLOAT',
    'NUMBER',
    'EXPONENTIAL_INTEGER_D',
    'EXPONENTIAL_FLOAT_D',
    'VARIABLE_STRING',
    'VARIABLE_NAME',
]
