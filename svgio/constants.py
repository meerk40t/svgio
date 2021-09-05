import re

SVGIO_VERSION = "0.0.1"

MIN_DEPTH = 5
ERROR = 1e-12

max_depth = 0

# SVG STATIC VALUES
DEFAULT_PPI = 96.0
SVG_NAME_TAG = "svg"
SVG_ATTR_VERSION = "version"
SVG_VALUE_VERSION = "1.1"
SVG_ATTR_XMLNS = "xmlns"
SVG_VALUE_XMLNS = "http://www.w3.org/2000/svg"
SVG_ATTR_XMLNS_LINK = "xmlns:xlink"
SVG_VALUE_XLINK = "http://www.w3.org/1999/xlink"
SVG_ATTR_XMLNS_EV = "xmlns:ev"
SVG_VALUE_XMLNS_EV = "http://www.w3.org/2001/xml-events"

XLINK_HREF = "{http://www.w3.org/1999/xlink}href"
SVG_HREF = "href"
SVG_ATTR_WIDTH = "width"
SVG_ATTR_HEIGHT = "height"
SVG_ATTR_VIEWBOX = "viewBox"
SVG_VIEWBOX_TRANSFORM = "viewbox_transform"
SVG_TAG_PATH = "path"
SVG_TAG_GROUP = "g"
SVG_TAG_RECT = "rect"
SVG_TAG_CIRCLE = "circle"
SVG_TAG_ELLIPSE = "ellipse"
SVG_TAG_LINE = "line"
SVG_TAG_POLYLINE = "polyline"
SVG_TAG_POLYGON = "polygon"
SVG_TAG_TEXT = "text"
SVG_TAG_TSPAN = "tspan"
SVG_TAG_IMAGE = "image"
SVG_TAG_DESC = "desc"
SVG_TAG_TITLE = "title"
SVG_TAG_METADATA = "metadata"
SVG_TAG_STYLE = "style"
SVG_TAG_DEFS = "defs"
SVG_TAG_USE = "use"
SVG_TAG_CLIPPATH = "clipPath"
SVG_TAG_PATTERN = "pattern"

SVG_STRUCT_ATTRIB = "attributes"
SVG_ATTR_ID = "id"
SVG_ATTR_DATA = "d"
SVG_ATTR_DISPLAY = "display"
SVG_ATTR_COLOR = "color"
SVG_ATTR_FILL = "fill"
SVG_ATTR_FILL_OPACITY = "fill-opacity"
SVG_ATTR_STROKE = "stroke"
SVG_ATTR_STROKE_OPACITY = "stroke-opacity"
SVG_ATTR_STROKE_WIDTH = "stroke-width"
SVG_ATTR_TRANSFORM = "transform"
SVG_ATTR_STYLE = "style"
SVG_ATTR_CLASS = "class"
SVG_ATTR_CLIP_PATH = "clip-path"
SVG_ATTR_CLIP_RULE = "clip-rule"
SVG_ATTR_CLIP_UNIT_TYPE = "clipPathUnits"
SVG_ATTR_CENTER_X = "cx"
SVG_ATTR_CENTER_Y = "cy"
SVG_ATTR_RADIUS_X = "rx"
SVG_ATTR_RADIUS_Y = "ry"
SVG_ATTR_RADIUS = "r"
SVG_ATTR_POINTS = "points"
SVG_ATTR_PRESERVEASPECTRATIO = "preserveAspectRatio"
SVG_ATTR_X = "x"
SVG_ATTR_Y = "y"
SVG_ATTR_X0 = "x0"
SVG_ATTR_Y0 = "y0"
SVG_ATTR_X1 = "x1"
SVG_ATTR_Y1 = "y1"
SVG_ATTR_X2 = "x2"
SVG_ATTR_Y2 = "y2"
SVG_ATTR_DX = "dx"
SVG_ATTR_DY = "dy"
SVG_ATTR_TAG = "tag"
SVG_ATTR_FONT = "font"
SVG_ATTR_FONT_FAMILY = "font-family"  # Serif, sans-serif, cursive, fantasy, monospace
SVG_ATTR_FONT_FACE = "font-face"
SVG_ATTR_FONT_SIZE = "font-size"
SVG_ATTR_FONT_WEIGHT = "font-weight"  # normal, bold, bolder, lighter, 100-900
SVG_ATTR_TEXT_ANCHOR = "text-anchor"
SVG_ATTR_PATTERN_CONTENT_UNITS = "patternContentUnits"
SVG_ATTR_PATTERN_TRANSFORM = "patternTransform"
SVG_ATTR_PATTERN_UNITS = "patternUnits"

SVG_ATTR_VECTOR_EFFECT = "vector-effect"

SVG_UNIT_TYPE_USERSPACEONUSE = "userSpaceOnUse"
SVG_UNIT_TYPE_OBJECTBOUNDINGBOX = "objectBoundingBox"

SVG_RULE_NONZERO = "nonzero"
SVG_RULE_EVENODD = "evenodd"

SVG_TRANSFORM_MATRIX = "matrix"
SVG_TRANSFORM_TRANSLATE = "translate"
SVG_TRANSFORM_SCALE = "scale"
SVG_TRANSFORM_ROTATE = "rotate"
SVG_TRANSFORM_SKEW_X = "skewx"
SVG_TRANSFORM_SKEW_Y = "skewy"
SVG_TRANSFORM_SKEW = "skew"
SVG_TRANSFORM_TRANSLATE_X = "translatex"
SVG_TRANSFORM_TRANSLATE_Y = "translatey"
SVG_TRANSFORM_SCALE_X = "scalex"
SVG_TRANSFORM_SCALE_Y = "scaley"

SVG_VALUE_NONE = "none"
SVG_VALUE_CURRENT_COLOR = "currentColor"

SVG_VALUE_NON_SCALING_STROKE = "non-scaling-stroke"

PATTERN_WS = r"[\s\t\n]*"
PATTERN_COMMA = r"(?:\s*,\s*|\s+|(?=-))"
PATTERN_COMMAWSP = r"[ ,\t\n\x09\x0A\x0C\x0D]+"
PATTERN_FLOAT = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
PATTERN_LENGTH_UNITS = "cm|mm|Q|in|pt|pc|px|em|cx|ch|rem|vw|vh|vmin|vmax"
PATTERN_ANGLE_UNITS = "deg|grad|rad|turn"
PATTERN_TIME_UNITS = "s|ms"
PATTERN_FREQUENCY_UNITS = "Hz|kHz"
PATTERN_RESOLUTION_UNITS = "dpi|dpcm|dppx"
PATTERN_PERCENT = "%"
PATTERN_TRANSFORM = (
    SVG_TRANSFORM_MATRIX
    + "|"
    + SVG_TRANSFORM_TRANSLATE
    + "|"
    + SVG_TRANSFORM_TRANSLATE_X
    + "|"
    + SVG_TRANSFORM_TRANSLATE_Y
    + "|"
    + SVG_TRANSFORM_SCALE
    + "|"
    + SVG_TRANSFORM_SCALE_X
    + "|"
    + SVG_TRANSFORM_SCALE_Y
    + "|"
    + SVG_TRANSFORM_ROTATE
    + "|"
    + SVG_TRANSFORM_SKEW
    + "|"
    + SVG_TRANSFORM_SKEW_X
    + "|"
    + SVG_TRANSFORM_SKEW_Y
)
PATTERN_TRANSFORM_UNITS = (
    PATTERN_LENGTH_UNITS + "|" + PATTERN_ANGLE_UNITS + "|" + PATTERN_PERCENT
)

REGEX_IRI = re.compile(r"url\(#?(.*)\)")
REGEX_FLOAT = re.compile(PATTERN_FLOAT)
REGEX_COORD_PAIR = re.compile(
    "(%s)%s(%s)" % (PATTERN_FLOAT, PATTERN_COMMA, PATTERN_FLOAT)
)
REGEX_TRANSFORM_TEMPLATE = re.compile(
    r"(?u)(%s)%s\(([^)]+)\)" % (PATTERN_TRANSFORM, PATTERN_WS)
)
REGEX_TRANSFORM_PARAMETER = re.compile(
    "(%s)%s(%s)?" % (PATTERN_FLOAT, PATTERN_WS, PATTERN_TRANSFORM_UNITS)
)
REGEX_COLOR_HEX = re.compile(r"^#?([0-9A-Fa-f]{3,8})$")
REGEX_COLOR_RGB = re.compile(
    r"rgba?\(\s*(%s)\s*,\s*(%s)\s*,\s*(%s)\s*(?:,\s*(%s)\s*)?\)"
    % (PATTERN_FLOAT, PATTERN_FLOAT, PATTERN_FLOAT, PATTERN_FLOAT)
)
REGEX_COLOR_RGB_PERCENT = re.compile(
    r"rgba?\(\s*(%s)%%\s*,\s*(%s)%%\s*,\s*(%s)%%\s*(?:,\s*(%s)\s*)?\)"
    % (PATTERN_FLOAT, PATTERN_FLOAT, PATTERN_FLOAT, PATTERN_FLOAT)
)
REGEX_COLOR_HSL = re.compile(
    r"hsla?\(\s*(%s)\s*,\s*(%s)%%\s*,\s*(%s)%%\s*(?:,\s*(%s)\s*)?\)"
    % (PATTERN_FLOAT, PATTERN_FLOAT, PATTERN_FLOAT, PATTERN_FLOAT)
)
REGEX_LENGTH = re.compile(r"(%s)([A-Za-z%%]*)" % PATTERN_FLOAT)
REGEX_CSS_STYLE = re.compile(r"([^{]+)\s*\{\s*([^}]+)\s*\}")
REGEX_CSS_FONT = re.compile(
    r"(?:(normal|italic|oblique)\s|(normal|small-caps)\s|(normal|bold|bolder|lighter|\d{3})\s|(normal|ultra-condensed|extra-condensed|condensed|semi-condensed|semi-expanded|expanded|extra-expanded|ultra-expanded)\s)*\s*(xx-small|x-small|small|medium|large|x-large|xx-large|larger|smaller|\d+(?:em|pt|pc|px|%))(?:/(xx-small|x-small|small|medium|large|x-large|xx-large|larger|smaller|\d+(?:em|pt|pc|px|%)))?\s*(.*),?\s+(serif|sans-serif|cursive|fantasy|monospace);?"
)

svg_parse = [("COMMAND", r"[MmZzLlHhVvCcSsQqTtAa]"), ("SKIP", PATTERN_COMMAWSP)]
svg_re = re.compile("|".join("(?P<%s>%s)" % pair for pair in svg_parse))
num_parse = [("FLOAT", PATTERN_FLOAT), ("CLOSE", r"[Zz]"), ("SKIP", PATTERN_COMMAWSP)]
num_re = re.compile("|".join("(?P<%s>%s)" % pair for pair in num_parse))
flag_parse = [("FLAG", r"[01]"), ("SKIP", PATTERN_COMMAWSP)]
flag_re = re.compile("|".join("(?P<%s>%s)" % pair for pair in flag_parse))
