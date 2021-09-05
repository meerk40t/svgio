from copy import copy
from math import (
    acos,
    atan,
    atan2,
    ceil,
    cos,
    degrees,
    hypot,
    log,
    radians,
    sin,
    sqrt,
    tan,
)
try:
    from math import tau
except ImportError:
    from math import pi

    tau = pi * 2
from .constants import *

########################
#  SVG/CSS Types
########################


class Length(object):
    """
    Length as used in SVG/CSS

    Length are lazy solving values. Several conversion values are typically unknown by default and length simply
    stores that ambiguity. So we can have a length of 50% and without calling .value(relative_length=3000) it will
    simply store as 50%. Likewise you can have discrete values like 30cm or 20in which have knowable discrete values
    but are not knowable in pixels unless a PPI value is supplied. We can say .value(relative_length=30cm, PPI=96) and
    solve this for a value like 12%. We can also convert values between knowable lengths. So 30cm in 300mm regardless
    whether we know how to convert this to pixels. 0% is 0 in any units or relative values. We can convert pixels to
    pc and pt without issue. We can convert vh, vw, vmax, vmin values if we know viewbox values. We can convert em
    values if we know the font_size. We can add values together if they are convertible units. Length("20in") + "3cm".

    If .value() cannot solve for the value with the given information then it will return a Length value. If it can
    be solved it will return a float.
    """

    def __init__(self, *args, **kwargs):
        if len(args) == 1:
            value = args[0]
            if value is None:
                self.amount = None
                self.units = None
                return
            s = str(value)
            for m in REGEX_LENGTH.findall(s):
                self.amount = float(m[0])
                self.units = m[1]
                return
        elif len(args) == 2:
            self.amount = args[0]
            self.units = args[1]
            return
        self.amount = 0.0
        self.units = ""

    def __float__(self):
        if self.amount is None:
            return None
        if self.units == "pt":
            return self.amount * 1.3333
        elif self.units == "pc":
            return self.amount * 16.0
        return self.amount

    def __imul__(self, other):
        if isinstance(other, (int, float)):
            self.amount *= other
            return self
        if self.amount == 0.0:
            return 0.0
        if isinstance(other, str):
            other = Length(other)
        if isinstance(other, Length):
            if other.amount == 0.0:
                self.amount = 0.0
                return self
            if self.units == other.units:
                self.amount *= other.amount
                return self
            if self.units == "%":
                self.units = other.units
                self.amount = self.amount * other.amount / 100.0
                return self
            elif other.units == "%":
                self.amount = self.amount * other.amount / 100.0
                return self
        raise ValueError

    def __iadd__(self, other):
        if not isinstance(other, Length):
            other = Length(other)
        if self.units == other.units:
            self.amount += other.amount
            return self
        if self.amount == 0:
            self.amount = other.amount
            self.units = other.units
            return self
        if other.amount == 0:
            return self
        if self.units == "px" or self.units == "":
            if other.units == "px" or other.units == "":
                self.amount += other.amount
            elif other.units == "pt":
                self.amount += other.amount * 1.3333
            elif other.units == "pc":
                self.amount += other.amount * 16.0
            else:
                raise ValueError
            return self
        if self.units == "pt":
            if other.units == "px" or other.units == "":
                self.amount += other.amount / 1.3333
            elif other.units == "pc":
                self.amount += other.amount * 12.0
            else:
                raise ValueError
            return self
        elif self.units == "pc":
            if other.units == "px" or other.units == "":
                self.amount += other.amount / 16.0
            elif other.units == "pt":
                self.amount += other.amount / 12.0
            else:
                raise ValueError
            return self
        elif self.units == "cm":
            if other.units == "mm":
                self.amount += other.amount / 10.0
            elif other.units == "in":
                self.amount += other.amount / 0.393701
            else:
                raise ValueError
            return self
        elif self.units == "mm":
            if other.units == "cm":
                self.amount += other.amount * 10.0
            elif other.units == "in":
                self.amount += other.amount / 0.0393701
            else:
                raise ValueError
            return self
        elif self.units == "in":
            if other.units == "cm":
                self.amount += other.amount * 0.393701
            elif other.units == "mm":
                self.amount += other.amount * 0.0393701
            else:
                raise ValueError
            return self
        raise ValueError("%s units were not determined." % self.units)

    def __abs__(self):
        c = self.__copy__()
        c.amount = abs(c.amount)
        return c

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            c = self.__copy__()
            c.amount /= other
            return c
        if self.amount == 0.0:
            return 0.0
        if isinstance(other, str):
            other = Length(other)
        if isinstance(other, Length):
            if self.units == other.units:
                q = self.amount / other.amount
                return q  # no units
        if self.units == "px" or self.units == "":
            if other.units == "px" or other.units == "":
                return self.amount / other.amount
            elif other.units == "pt":
                return self.amount / (other.amount * 1.3333)
            elif other.units == "pc":
                return self.amount / (other.amount * 16.0)
            else:
                raise ValueError
        if self.units == "pt":
            if other.units == "px" or other.units == "":
                return self.amount / (other.amount / 1.3333)
            elif other.units == "pc":
                return self.amount / (other.amount * 12.0)
            else:
                raise ValueError
        if self.units == "pc":
            if other.units == "px" or other.units == "":
                return self.amount / (other.amount / 16.0)
            elif other.units == "pt":
                return self.amount / (other.amount / 12.0)
            else:
                raise ValueError
        if self.units == "cm":
            if other.units == "mm":
                return self.amount / (other.amount / 10.0)
            elif other.units == "in":
                return self.amount / (other.amount / 0.393701)
            else:
                raise ValueError
        if self.units == "mm":
            if other.units == "cm":
                return self.amount / (other.amount * 10.0)
            elif other.units == "in":
                return self.amount / (other.amount / 0.0393701)
            else:
                raise ValueError
        if self.units == "in":
            if other.units == "cm":
                return self.amount / (other.amount * 0.393701)
            elif other.units == "mm":
                return self.amount / (other.amount * 0.0393701)
            else:
                raise ValueError
        raise ValueError

    __floordiv__ = __truediv__
    __div__ = __truediv__

    def __lt__(self, other):
        return (self - other).amount < 0.0

    def __le__(self, other):
        return (self - other).amount <= 0.0

    def __gt__(self, other):
        return (self - other).amount > 0.0

    def __ge__(self, other):
        return (self - other).amount >= 0.0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __add__(self, other):
        if isinstance(other, (str, float, int)):
            other = Length(other)
        c = self.__copy__()
        c += other
        return c

    __radd__ = __add__

    def __mul__(self, other):
        c = copy(self)
        c *= other
        return c

    def __rdiv__(self, other):
        c = copy(self)
        c *= 1.0 / other.amount
        return c

    def __neg__(self):
        s = self.__copy__()
        s.amount = -s.amount
        return s

    def __isub__(self, other):
        if isinstance(other, (str, float, int)):
            other = Length(other)
        self += -other
        return self

    def __sub__(self, other):
        s = self.__copy__()
        s -= other
        return s

    def __rsub__(self, other):
        if isinstance(other, (str, float, int)):
            other = Length(other)
        return (-self) + other

    def __copy__(self):
        return Length(self.amount, self.units)

    __rmul__ = __mul__

    def __repr__(self):
        return "Length('%s')" % (str(self))

    def __str__(self):
        if self.amount is None:
            return SVG_VALUE_NONE
        return "%s%s" % (Length.str(self.amount), self.units)

    def __eq__(self, other):
        if other is None:
            return False
        s = self.in_pixels()
        if isinstance(other, (float, int)):
            if s is not None:
                return abs(s - other) <= ERROR
            else:
                return other == 0 and self.amount == 0
        if isinstance(other, str):
            other = Length(other)
        if self.amount == other.amount and self.units == other.units:
            return True
        if s is not None:
            o = self.in_pixels()
            if abs(s - o) <= ERROR:
                return True
        s = self.in_inches()
        if s is not None:
            o = self.in_inches()
            if abs(s - o) <= ERROR:
                return True
        return False

    @property
    def value_in_units(self):
        return self.amount

    def in_pixels(self):
        if self.units == "px" or self.units == "":
            return self.amount
        if self.units == "pt":
            return self.amount / 1.3333
        if self.units == "pc":
            return self.amount / 16.0
        return None

    def in_inches(self):
        if self.units == "mm":
            return self.amount * 0.0393701
        if self.units == "cm":
            return self.amount * 0.393701
        if self.units == "in":
            return self.amount
        return None

    def to_mm(
        self,
        ppi=DEFAULT_PPI,
        relative_length=None,
        font_size=None,
        font_height=None,
        viewbox=None,
    ):
        value = self.value(
            ppi=ppi,
            relative_length=relative_length,
            font_size=font_size,
            font_height=font_height,
            viewbox=viewbox,
        )
        v = value / (ppi * 0.0393701)
        return Length("%smm" % (Length.str(v)))

    def to_cm(
        self,
        ppi=DEFAULT_PPI,
        relative_length=None,
        font_size=None,
        font_height=None,
        viewbox=None,
    ):
        value = self.value(
            ppi=ppi,
            relative_length=relative_length,
            font_size=font_size,
            font_height=font_height,
            viewbox=viewbox,
        )
        v = value / (ppi * 0.393701)
        return Length("%scm" % (Length.str(v)))

    def to_inch(
        self,
        ppi=DEFAULT_PPI,
        relative_length=None,
        font_size=None,
        font_height=None,
        viewbox=None,
    ):
        value = self.value(
            ppi=ppi,
            relative_length=relative_length,
            font_size=font_size,
            font_height=font_height,
            viewbox=viewbox,
        )
        v = value / ppi
        return Length("%sin" % (Length.str(v)))

    def value(
        self,
        ppi=None,
        relative_length=None,
        font_size=None,
        font_height=None,
        viewbox=None,
        **kwargs
    ):
        if self.amount is None:
            return None
        if self.units == "%":
            if relative_length is None:
                return self
            fraction = self.amount / 100.0
            if isinstance(relative_length, (float, int)):
                return fraction * relative_length
            elif isinstance(relative_length, (str, Length)):
                length = relative_length * self
                if isinstance(length, Length):
                    return length.value(
                        ppi=ppi,
                        font_size=font_size,
                        font_height=font_height,
                        viewbox=viewbox,
                    )
                return length
            return self
        if self.units == "mm":
            if ppi is None:
                return self
            return self.amount * ppi * 0.0393701
        if self.units == "cm":
            if ppi is None:
                return self
            return self.amount * ppi * 0.393701
        if self.units == "in":
            if ppi is None:
                return self
            return self.amount * ppi
        if self.units == "px" or self.units == "":
            return self.amount
        if self.units == "pt":
            return self.amount * 1.3333
        if self.units == "pc":
            return self.amount * 16.0
        if self.units == "em":
            if font_size is None:
                return self
            return self.amount * float(font_size)
        if self.units == "ex":
            if font_height is None:
                return self
            return self.amount * float(font_height)
        if self.units == "vw":
            if viewbox is None:
                return self
            v = Viewbox(viewbox)
            return self.amount * v.width / 100.0
        if self.units == "vh":
            if viewbox is None:
                return self
            v = Viewbox(viewbox)
            return self.amount * v.height / 100.0
        if self.units == "vmin":
            if viewbox is None:
                return self
            v = Viewbox(viewbox)
            m = min(v.height, v.height)
            return self.amount * m / 100.0
        if self.units == "vmax":
            if viewbox is None:
                return self
            v = Viewbox(viewbox)
            m = max(v.height, v.height)
            return self.amount * m / 100.0
        try:
            return float(self)
        except ValueError:
            return self

    @staticmethod
    def str(s):
        if s is None:
            return "n/a"
        if isinstance(s, Length):
            if s.units == "":
                s = s.amount
            else:
                a = "%.12f" % s.amount
                if "." in a:
                    a = a.rstrip("0").rstrip(".")
                return "'%s%s'" % (a, s.units)
        try:
            s = "%.12f" % s
        except TypeError:
            return str(s)
        if "." in s:
            s = s.rstrip("0").rstrip(".")
        return s


class Color(object):
    """
    SVG Color Parsing
    Parses different forms of defining colors.

    Including keyword: https://www.w3.org/TR/SVG11/types.html#ColorKeywords
    """

    def __init__(self, *args, **kwargs):
        self.value = 0
        arglen = len(args)
        if arglen == 1:
            v = args[0]
            if isinstance(v, Color):
                self.value = v.value
            elif isinstance(v, int):
                self.rgb = v
            else:
                self.value = Color.parse(v)
        elif arglen == 2:
            v = args[0]
            if isinstance(v, Color):
                self.value = v.value
            elif isinstance(v, int):
                self.rgb = v
            else:
                self.value = Color.parse(v)
            self.opacity = float(args[1])
        elif arglen == 3:
            r = args[0]
            g = args[1]
            b = args[2]
            self.value = Color.rgb_to_int(r, g, b)
        elif arglen == 4:
            r = args[0]
            g = args[1]
            b = args[2]
            opacity = args[3] / 255.0
            self.value = Color.rgb_to_int(r, g, b, opacity)
        if "red" in kwargs:
            self.red = kwargs["red"]
        if "green" in kwargs:
            self.green = kwargs["green"]
        if "blue" in kwargs:
            self.blue = kwargs["blue"]
        if "alpha" in kwargs:
            self.alpha = kwargs["alpha"]
        if "opacity" in kwargs:
            self.opacity = kwargs["opacity"]
        if "r" in kwargs:
            self.red = kwargs["r"]
        if "g" in kwargs:
            self.green = kwargs["g"]
        if "b" in kwargs:
            self.blue = kwargs["b"]
        if "rgb" in kwargs:
            self.rgb = kwargs["rgb"]
        if "argb" in kwargs:
            self.argb = kwargs["argb"]
        if "rgba" in kwargs:
            self.rgba = kwargs["rgba"]
        if "h" in kwargs:
            self.hue = kwargs["h"]
        if "s" in kwargs:
            self.saturation = kwargs["s"]
        if "l" in kwargs:
            self.lightness = kwargs["l"]
        if "hue" in kwargs:
            self.hue = kwargs["hue"]
        if "saturation" in kwargs:
            self.saturation = kwargs["saturation"]
        if "lightness" in kwargs:
            self.lightness = kwargs["lightness"]

    def __int__(self):
        return self.value

    def __str__(self):
        if self.value is None:
            return str(self.value)
        return self.hex

    def __repr__(self):
        if self.value is None:
            return "Color('%s')" % self.value
        return "Color('%s')" % self.hex

    def __eq__(self, other):
        if self is other:
            return True
        first = self.value
        second = other
        if isinstance(second, str):
            second = Color(second)
        if isinstance(second, Color):
            second = second.value
        if first is None:
            return second is None
        if second is None:
            return first is None
        return first & 0xFFFFFFFF == second & 0xFFFFFFFF

    def __ne__(self, other):
        return not self == other

    def __abs__(self):
        # Return opaque color.
        if self.value is None:
            return Color(self.value)
        return Color(self.red, self.green, self.blue)

    @staticmethod
    def rgb_to_int(r, g, b, opacity=1.0):
        if opacity > 1:
            opacity = 1.0
        if opacity < 0:
            opacity = 0
        r = Color.crimp(r)
        g = Color.crimp(g)
        b = Color.crimp(b)
        a = int(round(opacity * 255.0))
        a = Color.crimp(a)
        r <<= 24
        g <<= 16
        b <<= 8
        c = r | g | b | a
        return c

    @staticmethod
    def hsl_to_int(h, s, l, opacity=1.0):
        def hue_2_rgb(v1, v2, vh):
            if vh < 0:
                vh += 1
            if vh > 1:
                vh -= 1
            if 6.0 * vh < 1.0:
                return v1 + (v2 - v1) * 6.0 * vh
            if 2.0 * vh < 1:
                return v2
            if 3 * vh < 2.0:
                return v1 + (v2 - v1) * ((2.0 / 3.0) - vh) * 6.0
            return v1

        if s == 0.0:
            r = 255.0 * l
            g = 255.0 * l
            b = 255.0 * l
        else:
            if l < 0.5:
                v2 = l * (1.0 + s)
            else:
                v2 = (l + s) - (s * l)
            v1 = 2 * l - v2
            r = 255.0 * hue_2_rgb(v1, v2, h + (1.0 / 3.0))
            g = 255.0 * hue_2_rgb(v1, v2, h)
            b = 255.0 * hue_2_rgb(v1, v2, h - (1.0 / 3.0))
        value = Color.rgb_to_int(r, g, b, opacity=opacity)
        return value

    @staticmethod
    def parse(color_string):
        """Parse SVG color, will return a set value."""
        if color_string is None or color_string == SVG_VALUE_NONE:
            return None
        match = REGEX_COLOR_HEX.match(color_string)
        if match:
            return Color.parse_color_hex(color_string)
        match = REGEX_COLOR_RGB.match(color_string)
        if match:
            return Color.parse_color_rgb(match.groups())
        match = REGEX_COLOR_RGB_PERCENT.match(color_string)
        if match:
            return Color.parse_color_rgbp(match.groups())
        match = REGEX_COLOR_HSL.match(color_string)
        if match:
            return Color.parse_color_hsl(match.groups())
        return Color.parse_color_lookup(color_string)

    @staticmethod
    def parse_color_lookup(v):
        """Parse SVG Color by Keyword on dictionary lookup"""
        if not isinstance(v, str):
            return Color.rgb_to_int(0, 0, 0)
        else:
            v = v.replace(" ", "").lower()
        if v == "transparent":
            return Color.rgb_to_int(0, 0, 0, 0.0)
        if v == "aliceblue":
            return Color.rgb_to_int(250, 248, 255)
        if v == "aliceblue":
            return Color.rgb_to_int(240, 248, 255)
        if v == "antiquewhite":
            return Color.rgb_to_int(250, 235, 215)
        if v == "aqua":
            return Color.rgb_to_int(0, 255, 255)
        if v == "aquamarine":
            return Color.rgb_to_int(127, 255, 212)
        if v == "azure":
            return Color.rgb_to_int(240, 255, 255)
        if v == "beige":
            return Color.rgb_to_int(245, 245, 220)
        if v == "bisque":
            return Color.rgb_to_int(255, 228, 196)
        if v == "black":
            return Color.rgb_to_int(0, 0, 0)
        if v == "blanchedalmond":
            return Color.rgb_to_int(255, 235, 205)
        if v == "blue":
            return Color.rgb_to_int(0, 0, 255)
        if v == "blueviolet":
            return Color.rgb_to_int(138, 43, 226)
        if v == "brown":
            return Color.rgb_to_int(165, 42, 42)
        if v == "burlywood":
            return Color.rgb_to_int(222, 184, 135)
        if v == "cadetblue":
            return Color.rgb_to_int(95, 158, 160)
        if v == "chartreuse":
            return Color.rgb_to_int(127, 255, 0)
        if v == "chocolate":
            return Color.rgb_to_int(210, 105, 30)
        if v == "coral":
            return Color.rgb_to_int(255, 127, 80)
        if v == "cornflowerblue":
            return Color.rgb_to_int(100, 149, 237)
        if v == "cornsilk":
            return Color.rgb_to_int(255, 248, 220)
        if v == "crimson":
            return Color.rgb_to_int(220, 20, 60)
        if v == "cyan":
            return Color.rgb_to_int(0, 255, 255)
        if v == "darkblue":
            return Color.rgb_to_int(0, 0, 139)
        if v == "darkcyan":
            return Color.rgb_to_int(0, 139, 139)
        if v == "darkgoldenrod":
            return Color.rgb_to_int(184, 134, 11)
        if v == "darkgray":
            return Color.rgb_to_int(169, 169, 169)
        if v == "darkgreen":
            return Color.rgb_to_int(0, 100, 0)
        if v == "darkgrey":
            return Color.rgb_to_int(169, 169, 169)
        if v == "darkkhaki":
            return Color.rgb_to_int(189, 183, 107)
        if v == "darkmagenta":
            return Color.rgb_to_int(139, 0, 139)
        if v == "darkolivegreen":
            return Color.rgb_to_int(85, 107, 47)
        if v == "darkorange":
            return Color.rgb_to_int(255, 140, 0)
        if v == "darkorchid":
            return Color.rgb_to_int(153, 50, 204)
        if v == "darkred":
            return Color.rgb_to_int(139, 0, 0)
        if v == "darksalmon":
            return Color.rgb_to_int(233, 150, 122)
        if v == "darkseagreen":
            return Color.rgb_to_int(143, 188, 143)
        if v == "darkslateblue":
            return Color.rgb_to_int(72, 61, 139)
        if v == "darkslategray":
            return Color.rgb_to_int(47, 79, 79)
        if v == "darkslategrey":
            return Color.rgb_to_int(47, 79, 79)
        if v == "darkturquoise":
            return Color.rgb_to_int(0, 206, 209)
        if v == "darkviolet":
            return Color.rgb_to_int(148, 0, 211)
        if v == "deeppink":
            return Color.rgb_to_int(255, 20, 147)
        if v == "deepskyblue":
            return Color.rgb_to_int(0, 191, 255)
        if v == "dimgray":
            return Color.rgb_to_int(105, 105, 105)
        if v == "dimgrey":
            return Color.rgb_to_int(105, 105, 105)
        if v == "dodgerblue":
            return Color.rgb_to_int(30, 144, 255)
        if v == "firebrick":
            return Color.rgb_to_int(178, 34, 34)
        if v == "floralwhite":
            return Color.rgb_to_int(255, 250, 240)
        if v == "forestgreen":
            return Color.rgb_to_int(34, 139, 34)
        if v == "fuchsia":
            return Color.rgb_to_int(255, 0, 255)
        if v == "gainsboro":
            return Color.rgb_to_int(220, 220, 220)
        if v == "ghostwhite":
            return Color.rgb_to_int(248, 248, 255)
        if v == "gold":
            return Color.rgb_to_int(255, 215, 0)
        if v == "goldenrod":
            return Color.rgb_to_int(218, 165, 32)
        if v == "gray":
            return Color.rgb_to_int(128, 128, 128)
        if v == "grey":
            return Color.rgb_to_int(128, 128, 128)
        if v == "green":
            return Color.rgb_to_int(0, 128, 0)
        if v == "greenyellow":
            return Color.rgb_to_int(173, 255, 47)
        if v == "honeydew":
            return Color.rgb_to_int(240, 255, 240)
        if v == "hotpink":
            return Color.rgb_to_int(255, 105, 180)
        if v == "indianred":
            return Color.rgb_to_int(205, 92, 92)
        if v == "indigo":
            return Color.rgb_to_int(75, 0, 130)
        if v == "ivory":
            return Color.rgb_to_int(255, 255, 240)
        if v == "khaki":
            return Color.rgb_to_int(240, 230, 140)
        if v == "lavender":
            return Color.rgb_to_int(230, 230, 250)
        if v == "lavenderblush":
            return Color.rgb_to_int(255, 240, 245)
        if v == "lawngreen":
            return Color.rgb_to_int(124, 252, 0)
        if v == "lemonchiffon":
            return Color.rgb_to_int(255, 250, 205)
        if v == "lightblue":
            return Color.rgb_to_int(173, 216, 230)
        if v == "lightcoral":
            return Color.rgb_to_int(240, 128, 128)
        if v == "lightcyan":
            return Color.rgb_to_int(224, 255, 255)
        if v == "lightgoldenrodyellow":
            return Color.rgb_to_int(250, 250, 210)
        if v == "lightgray":
            return Color.rgb_to_int(211, 211, 211)
        if v == "lightgreen":
            return Color.rgb_to_int(144, 238, 144)
        if v == "lightgrey":
            return Color.rgb_to_int(211, 211, 211)
        if v == "lightpink":
            return Color.rgb_to_int(255, 182, 193)
        if v == "lightsalmon":
            return Color.rgb_to_int(255, 160, 122)
        if v == "lightseagreen":
            return Color.rgb_to_int(32, 178, 170)
        if v == "lightskyblue":
            return Color.rgb_to_int(135, 206, 250)
        if v == "lightslategray":
            return Color.rgb_to_int(119, 136, 153)
        if v == "lightslategrey":
            return Color.rgb_to_int(119, 136, 153)
        if v == "lightsteelblue":
            return Color.rgb_to_int(176, 196, 222)
        if v == "lightyellow":
            return Color.rgb_to_int(255, 255, 224)
        if v == "lime":
            return Color.rgb_to_int(0, 255, 0)
        if v == "limegreen":
            return Color.rgb_to_int(50, 205, 50)
        if v == "linen":
            return Color.rgb_to_int(250, 240, 230)
        if v == "magenta":
            return Color.rgb_to_int(255, 0, 255)
        if v == "maroon":
            return Color.rgb_to_int(128, 0, 0)
        if v == "mediumaquamarine":
            return Color.rgb_to_int(102, 205, 170)
        if v == "mediumblue":
            return Color.rgb_to_int(0, 0, 205)
        if v == "mediumorchid":
            return Color.rgb_to_int(186, 85, 211)
        if v == "mediumpurple":
            return Color.rgb_to_int(147, 112, 219)
        if v == "mediumseagreen":
            return Color.rgb_to_int(60, 179, 113)
        if v == "mediumslateblue":
            return Color.rgb_to_int(123, 104, 238)
        if v == "mediumspringgreen":
            return Color.rgb_to_int(0, 250, 154)
        if v == "mediumturquoise":
            return Color.rgb_to_int(72, 209, 204)
        if v == "mediumvioletred":
            return Color.rgb_to_int(199, 21, 133)
        if v == "midnightblue":
            return Color.rgb_to_int(25, 25, 112)
        if v == "mintcream":
            return Color.rgb_to_int(245, 255, 250)
        if v == "mistyrose":
            return Color.rgb_to_int(255, 228, 225)
        if v == "moccasin":
            return Color.rgb_to_int(255, 228, 181)
        if v == "navajowhite":
            return Color.rgb_to_int(255, 222, 173)
        if v == "navy":
            return Color.rgb_to_int(0, 0, 128)
        if v == "oldlace":
            return Color.rgb_to_int(253, 245, 230)
        if v == "olive":
            return Color.rgb_to_int(128, 128, 0)
        if v == "olivedrab":
            return Color.rgb_to_int(107, 142, 35)
        if v == "orange":
            return Color.rgb_to_int(255, 165, 0)
        if v == "orangered":
            return Color.rgb_to_int(255, 69, 0)
        if v == "orchid":
            return Color.rgb_to_int(218, 112, 214)
        if v == "palegoldenrod":
            return Color.rgb_to_int(238, 232, 170)
        if v == "palegreen":
            return Color.rgb_to_int(152, 251, 152)
        if v == "paleturquoise":
            return Color.rgb_to_int(175, 238, 238)
        if v == "palevioletred":
            return Color.rgb_to_int(219, 112, 147)
        if v == "papayawhip":
            return Color.rgb_to_int(255, 239, 213)
        if v == "peachpuff":
            return Color.rgb_to_int(255, 218, 185)
        if v == "peru":
            return Color.rgb_to_int(205, 133, 63)
        if v == "pink":
            return Color.rgb_to_int(255, 192, 203)
        if v == "plum":
            return Color.rgb_to_int(221, 160, 221)
        if v == "powderblue":
            return Color.rgb_to_int(176, 224, 230)
        if v == "purple":
            return Color.rgb_to_int(128, 0, 128)
        if v == "red":
            return Color.rgb_to_int(255, 0, 0)
        if v == "rosybrown":
            return Color.rgb_to_int(188, 143, 143)
        if v == "royalblue":
            return Color.rgb_to_int(65, 105, 225)
        if v == "saddlebrown":
            return Color.rgb_to_int(139, 69, 19)
        if v == "salmon":
            return Color.rgb_to_int(250, 128, 114)
        if v == "sandybrown":
            return Color.rgb_to_int(244, 164, 96)
        if v == "seagreen":
            return Color.rgb_to_int(46, 139, 87)
        if v == "seashell":
            return Color.rgb_to_int(255, 245, 238)
        if v == "sienna":
            return Color.rgb_to_int(160, 82, 45)
        if v == "silver":
            return Color.rgb_to_int(192, 192, 192)
        if v == "skyblue":
            return Color.rgb_to_int(135, 206, 235)
        if v == "slateblue":
            return Color.rgb_to_int(106, 90, 205)
        if v == "slategray":
            return Color.rgb_to_int(112, 128, 144)
        if v == "slategrey":
            return Color.rgb_to_int(112, 128, 144)
        if v == "snow":
            return Color.rgb_to_int(255, 250, 250)
        if v == "springgreen":
            return Color.rgb_to_int(0, 255, 127)
        if v == "steelblue":
            return Color.rgb_to_int(70, 130, 180)
        if v == "tan":
            return Color.rgb_to_int(210, 180, 140)
        if v == "teal":
            return Color.rgb_to_int(0, 128, 128)
        if v == "thistle":
            return Color.rgb_to_int(216, 191, 216)
        if v == "tomato":
            return Color.rgb_to_int(255, 99, 71)
        if v == "turquoise":
            return Color.rgb_to_int(64, 224, 208)
        if v == "violet":
            return Color.rgb_to_int(238, 130, 238)
        if v == "wheat":
            return Color.rgb_to_int(245, 222, 179)
        if v == "white":
            return Color.rgb_to_int(255, 255, 255)
        if v == "whitesmoke":
            return Color.rgb_to_int(245, 245, 245)
        if v == "yellow":
            return Color.rgb_to_int(255, 255, 0)
        if v == "yellowgreen":
            return Color.rgb_to_int(154, 205, 50)
        try:
            return int(v)
        except ValueError:
            return Color.rgb_to_int(0, 0, 0)

    @staticmethod
    def parse_color_hex(hex_string):
        """Parse SVG Color by Hex String"""
        h = hex_string.lstrip("#")
        size = len(h)
        if size == 8:
            return int(h[:8], 16)
        elif size == 6:
            s = "{0}FF".format(h[:6])
            v = int(s, 16)
            return v
        elif size == 4:
            s = h[0] + h[0] + h[1] + h[1] + h[2] + h[2] + h[3] + h[3]
            return int(s, 16)
        elif size == 3:
            s = "{0}{0}{1}{1}{2}{2}FF".format(h[0], h[1], h[2])
            v = int(s, 16)
            return v
        return Color.rgb_to_int(0, 0, 0)

    @staticmethod
    def parse_color_rgb(values):
        """Parse SVG Color, RGB value declarations """
        r = int(values[0])
        g = int(values[1])
        b = int(values[2])
        if values[3] is not None:
            opacity = float(values[3])
        else:
            opacity = 1
        return Color.rgb_to_int(r, g, b, opacity)

    @staticmethod
    def parse_color_rgbp(values):
        """Parse SVG color, RGB percent value declarations"""
        ratio = 255.0 / 100.0
        r = round(float(values[0]) * ratio)
        g = round(float(values[1]) * ratio)
        b = round(float(values[2]) * ratio)
        if values[3] is not None:
            opacity = float(values[3])
        else:
            opacity = 1
        return Color.rgb_to_int(r, g, b, opacity)

    @staticmethod
    def parse_color_hsl(values):
        """Parse SVG color, HSL value declarations"""
        h = Angle.parse(values[0])
        h = h.as_turns
        s = float(values[1]) / 100.0
        if s > 1:
            s = 1.0
        if s < 0:
            s = 0.0
        l = float(values[2]) / 100.0
        if l > 1:
            l = 1.0
        if l < 0:
            l = 0.0
        if values[3] is not None:
            opacity = float(values[3])
        else:
            opacity = 1
        return Color.hsl_to_int(h, s, l, opacity)

    @classmethod
    def distinct(cls, index):
        """
        Produces a deterministic distinct color for the given index.
        """

        def _pattern(pattern: int):
            n = int(pattern ** (1.0 / 3.0))
            pattern -= n * n * n
            p = [n] * 3
            if pattern == 0:
                return p
            pattern -= 1

            v = int(pattern % 3)
            pattern = int(pattern // 3)
            if pattern < n:
                p[v] = pattern % n
                return p
            pattern -= n

            p[v] = pattern // n
            v += 1
            p[v % 3] = pattern % n
            return p

        def _8bit_reverse(r: int):
            value = r - 1
            v = 0
            for i in range(0, 8):
                v = v | (value & 1)
                v <<= 1
                value >>= 1
            v >>= 1
            return v & 0xFF

        p = _pattern(index)
        return Color(
            _8bit_reverse(p[0]),
            _8bit_reverse(p[1]),
            _8bit_reverse(p[2]),
        )

    @property
    def rgb(self):
        if self.value is None:
            return None
        return self.value >> 8

    @rgb.setter
    def rgb(self, rgb):
        rgb <<= 8
        rgb |= 0xFF
        self.value = rgb

    @property
    def rgba(self):
        return self.value

    @rgba.setter
    def rgba(self, rgba):
        self.value = rgba

    @property
    def argb(self):
        if self.value is None:
            return None
        return ((self.value >> 8) & 0xFFFFFF) | (self.alpha << 24)

    @argb.setter
    def argb(self, argb):
        self.value = ((argb << 8) & 0xFFFFFF00) | (argb >> 24 & 0xFF)

    @property
    def opacity(self):
        return self.alpha / 255.0 if self.value is not None else None

    @opacity.setter
    def opacity(self, opacity):
        if self.value is None:
            raise ValueError
        a = int(round(opacity * 255.0))
        a = Color.crimp(a)
        self.alpha = a

    @property
    def alpha(self):
        return self.value & 0xFF if self.value is not None else None

    @alpha.setter
    def alpha(self, a):
        if self.value is None:
            raise ValueError
        a = Color.crimp(a)
        self.value &= ~0xFF
        self.value |= a

    @property
    def red(self):
        return (self.value >> 24) & 0xFF if self.value is not None else None

    @red.setter
    def red(self, r):
        if self.value is None:
            raise ValueError
        r = Color.crimp(r)
        self.value &= ~0xFF000000
        self.value |= r << 24

    @property
    def green(self):
        return (self.value >> 16) & 0xFF if self.value is not None else None

    @green.setter
    def green(self, g):
        if self.value is None:
            raise ValueError
        g = Color.crimp(g)
        self.value &= ~0xFF0000
        self.value |= g << 16

    @property
    def blue(self):
        return (self.value >> 8) & 0xFF if self.value is not None else None

    @blue.setter
    def blue(self, b):
        if self.value is None:
            raise ValueError
        b = Color.crimp(b)
        self.value &= ~0xFF00
        self.value |= b << 8

    @property
    def hexa(self):
        return (
            "#%02x%02x%02x%02x" % (self.red, self.green, self.blue, self.alpha)
            if self.value is not None
            else None
        )

    @property
    def hex(self):
        if self.alpha == 0xFF:
            return (
                "#%02x%02x%02x" % (self.red, self.green, self.blue)
                if self.value is not None
                else None
            )
        else:
            return self.hexa

    @property
    def hue(self):
        if self.value is None:
            return None
        r = self.red / 255.0
        g = self.green / 255.0
        b = self.blue / 255.0
        var_min = min(r, g, b)
        var_max = max(r, g, b)
        delta_max = var_max - var_min
        if delta_max == 0:
            return 0
        dr = (((var_max - r) / 6.0) + delta_max / 2.0) / delta_max
        dg = (((var_max - g) / 6.0) + delta_max / 2.0) / delta_max
        db = (((var_max - b) / 6.0) + delta_max / 2.0) / delta_max
        if r == var_max:
            h = db - dg
        elif g == var_max:
            h = (1.0 / 3.0) + dr - db
        else:  # db == max_v
            h = (2.0 / 3.0) + dg - dr
        if h < 0:
            h += 1
        if h > 1:
            h -= 1
        return Angle.turns(h).as_degrees

    @hue.setter
    def hue(self, v):
        if self.value is None:
            raise ValueError
        h, s, l = self.hsl
        self.hsl = v, s, l

    @property
    def saturation(self):
        if self.value is None:
            return None
        r = self.red / 255.0
        g = self.green / 255.0
        b = self.blue / 255.0
        min_v = min(r, g, b)
        max_v = max(r, g, b)
        delta = max_v - min_v
        if max_v == min_v:
            return 0.0
        if (max_v + min_v) < 1:
            return delta / (max_v + min_v)
        else:
            return delta / (2.0 - max_v - min_v)

    @saturation.setter
    def saturation(self, v):
        if self.value is None:
            raise ValueError
        h, s, l = self.hsl
        self.hsl = h, v, l

    @property
    def lightness(self):
        if self.value is None:
            return None
        r = self.red / 255.0
        g = self.green / 255.0
        b = self.blue / 255.0
        min_v = min(r, g, b)
        max_v = max(r, g, b)
        return (max_v + min_v) / 2.0

    @lightness.setter
    def lightness(self, v):
        if self.value is None:
            raise ValueError
        h, s, l = self.hsl
        self.hsl = h, s, v

    @property
    def intensity(self):
        if self.value is None:
            return None
        r = self.red
        g = self.green
        b = self.blue
        return (r + b + g) / 768.0

    @property
    def brightness(self):
        if self.value is None:
            return None
        r = self.red
        g = self.green
        b = self.blue
        cmax = max(r, g, b)
        return cmax / 255.0

    @property
    def blackness(self):
        if self.value is None:
            return None
        return 1.0 - self.brightness

    @property
    def luminance(self):
        if self.value is None:
            return None
        r = self.red / 255.0
        g = self.green / 255.0
        b = self.blue / 255.0
        return r * 0.3 + g * 0.59 + b * 0.11

    @property
    def luma(self):
        if self.value is None:
            return None
        r = self.red / 255.0
        g = self.green / 255.0
        b = self.blue / 255.0
        return r * 0.2126 + g * 0.7152 + b * 0.0722

    @staticmethod
    def over(c1, c2):
        """
        Porter Duff Alpha compositing operation over.
        Returns c1 over c2. This is the standard painter algorithm.
        """
        if isinstance(c1, str):
            c1 = Color.parse(c1)
        elif isinstance(c1, int):
            c1 = Color(c1)
        if isinstance(c2, str):
            c2 = Color.parse(c2)
        elif isinstance(c2, int):
            c2 = Color(c2)
        r1 = c1.red
        g1 = c1.green
        b1 = c1.blue
        a1 = c1.alpha
        if a1 == 255:
            return c1.value
        if a1 == 0:
            return c2.value
        r2 = c2.red
        g2 = c2.green
        b2 = c2.blue
        a2 = c2.alpha

        q = 255.0 - a1

        sr = r1 * a1 * 255.0 + r2 * a2 * q
        sg = g1 * a1 * 255.0 + g2 * a2 * q
        sb = b1 * a1 * 255.0 + b2 * a2 * q
        sa = a1 * 255.0 + a2 * q
        sr /= sa
        sg /= sa
        sb /= sa
        sa /= 255.0 * 255.0
        return Color.rgb_to_int(sr, sg, sb, sa)

    @staticmethod
    def distance(c1, c2):
        return sqrt(Color.distance_sq(c1, c2))

    @staticmethod
    def distance_sq(c1, c2):
        """
        Function returns the square of colordistance. The square of the color distance will always be closer than the
        square of another color distance.

        Rather than naive Euclidean distance we use Compuphase's Redmean color distance.
        https://www.compuphase.com/cmetric.htm

        It's computationally simple, and empirical tests finds it to be on par with LabDE2000.

        :param c1: first color
        :param c2: second color
        :return: square of color distance
        """
        if isinstance(c1, str):
            c1 = Color(c1)
        elif isinstance(c1, int):
            c1 = Color(c1)
        if isinstance(c2, str):
            c2 = Color(c2)
        elif isinstance(c2, int):
            c2 = Color(c2)
        red_mean = int((c1.red + c2.red) / 2.0)
        r = c1.red - c2.red
        g = c1.green - c2.green
        b = c1.blue - c2.blue
        return (((512 + red_mean) * r * r) >> 8) + 4 * g * g + (
            (767 - red_mean) * b * b
        ) >> 8

    @staticmethod
    def crimp(v):
        if v > 255:
            return 255
        if v < 0:
            return 0
        return int(v)

    @property
    def hsl(self):
        if self.value is None:
            return None
        return self.hue, self.saturation, self.lightness

    @hsl.setter
    def hsl(self, value):
        if not isinstance(value, (tuple, list)):
            return
        h = value[0]
        s = value[1]
        l = value[2]
        self.value = Color.hsl_to_int(h, s, l, 1.0)

    def distance_to(self, other):
        return Color.distance(self, other)

    def blend(self, other, opacity=None):
        """
        Blends the given color with the current color.
        """
        if opacity is None:
            self.value = Color.over(other, self)
        else:
            color = Color(other)
            color.opacity = opacity
            self.value = Color.over(color, self)


class Point:
    """Point is a general subscriptable point class with .x and .y as well as [0] and [1]

    For compatibility with regebro svg.path we accept complex numbers as points x + yj,
    and provide .real and .imag as properties. As well as float and integer values as (v,0) elements.

    With regard to SVG 7.15.1 defining SVGPoint this class provides for matrix transformations.

    Points are only positions in real Euclidean space. This class is not intended to interact with
    the Length class.
    """

    def __init__(self, x, y=None):
        if x is not None and y is None:
            if isinstance(x, str):
                string_x, string_y = REGEX_COORD_PAIR.findall(x)[0]
                self.x = float(string_x)
                self.y = float(string_y)
                return
            try:  # Try .x .y
                self.y = x.y
                self.x = x.x
                return
            except AttributeError:
                pass
            try:  # try subscription.
                self.y = x[1]
                self.x = x[0]
                return
            except TypeError:
                pass
            try:  # try .imag .real complex values.
                self.y = x.imag
                self.x = x.real
                return
            except AttributeError:
                # Unknown.
                raise TypeError
        self.x = x
        self.y = y

    def __key(self):
        return (self.x, self.y)

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        if other is None:
            return False
        try:
            if not isinstance(other, Point):
                other = Point(other)
        except Exception:
            return NotImplemented

        return abs(self.x - other.x) <= ERROR and abs(self.y - other.y) <= ERROR

    def __ne__(self, other):
        return not self == other

    def __len__(self):
        return 2

    def __getitem__(self, item):
        if item == 0:
            return self.x
        elif item == 1:
            return self.y
        else:
            raise IndexError

    def __setitem__(self, key, value):
        if key == 0:
            self.x = value
        elif key == 1:
            self.y = value
        else:
            raise IndexError

    def __repr__(self):
        x_str = Length.str(self.x)
        y_str = Length.str(self.y)
        return "Point(%s,%s)" % (x_str, y_str)

    def __copy__(self):
        return Point(self.x, self.y)

    def __str__(self):
        try:
            x_str = "%.12G" % self.x
        except TypeError:
            return self.__repr__()
        if "." in x_str:
            x_str = x_str.rstrip("0").rstrip(".")
        y_str = "%.12G" % self.y
        if "." in y_str:
            y_str = y_str.rstrip("0").rstrip(".")
        return "%s,%s" % (x_str, y_str)

    def __imul__(self, other):
        if isinstance(other, str):
            other = Matrix(other)
        if isinstance(other, Matrix):
            v = other.point_in_matrix_space(self)
            self.x = v.x
            self.y = v.y
            return self
        try:
            c = complex(self) * complex(other.x, other.y)
            self.x = c.real
            self.y = c.imag
            return self
        except AttributeError:
            pass
        try:
            c = complex(self) * complex(other[0], other[1])
            self.x = c.real
            self.y = c.imag
            return self
        except (TypeError, IndexError):
            pass
        try:
            c = complex(self) * complex(other.real, other.imag)
            self.x = c.real
            self.y = c.imag
            return self
        except AttributeError:
            pass
        try:
            self.x *= other
            self.y *= other
            return self
        except Exception:
            return NotImplemented

    def __mul__(self, other):
        if isinstance(other, str):
            other = Matrix(other)
        if isinstance(other, Matrix):
            return other.point_in_matrix_space(self)
        try:
            return Point(complex(self) * complex(other.x, other.y))
        except AttributeError:
            pass
        try:
            return Point(complex(self) * complex(other[0], other[1]))
        except (TypeError, IndexError):
            pass
        try:
            return Point(complex(self) * complex(other.real, other.imag))
        except AttributeError:
            pass
        try:
            return Point(self.x * other, self.y * other)
        except Exception:
            return NotImplemented

    __rmul__ = __mul__

    def __iadd__(self, other):
        try:
            self.x += other.x
            self.y += other.y
            return self
        except AttributeError:
            pass
        try:
            self.y += other[1]
            self.x += other[0]
            return self
        except (TypeError, IndexError):
            pass
        try:
            self.x += other.real
            self.y += other.imag
            return self
        except AttributeError:
            pass
        try:
            self.x += other
            return self
        except Exception:
            return NotImplemented

    def __add__(self, other):
        try:
            x = self.x + other.x
            y = self.y + other.y
            return Point(x, y)
        except AttributeError:
            pass
        try:
            y = self.y + other[1]
            x = self.x + other[0]
            return Point(x, y)
        except (TypeError, IndexError):
            pass
        try:
            x = self.x + other.real
            y = self.y + other.imag
            return Point(x, y)
        except AttributeError:
            pass
        if isinstance(other, (float, int)):
            x = self.x + other
            return Point(x, self.y)
        return NotImplemented

    __radd__ = __add__

    def __isub__(self, other):
        try:
            self.x -= other.x
            self.y -= other.y
            return self
        except AttributeError:
            pass
        try:
            self.y -= other[1]
            self.x -= other[0]
            return self
        except (TypeError, IndexError):
            pass
        try:
            self.x -= other.real
            self.y -= other.imag
            return self
        except AttributeError:
            pass
        try:
            self.x -= other
            return self
        except Exception:
            return NotImplemented

    def __sub__(self, other):
        try:
            x = self.x - other.x
            y = self.y - other.y
            return Point(x, y)
        except AttributeError:
            pass
        try:
            y = self.y - other[1]
            x = self.x - other[0]
            return Point(x, y)
        except (TypeError, IndexError):
            pass
        try:
            x = self.x - other.real
            y = self.y - other.imag
            return Point(x, y)
        except AttributeError:
            pass
        if isinstance(other, (float, int)):
            x = self.x - other
            return Point(x, self.y)
        return NotImplemented

    def __rsub__(self, other):
        try:
            x = other.x - self.x
            y = other.y - self.y
            return Point(x, y)
        except AttributeError:
            pass
        try:
            y = other[1] - self.y
            x = other[0] - self.x
            return Point(x, y)
        except (TypeError, IndexError):
            pass
        try:
            x = other.real - self.x
            y = other.imag - self.y
            return Point(x, y)
        except AttributeError:
            pass
        if isinstance(other, (float, int)):
            x = other - self.x
            return Point(x, self.y)
        return NotImplemented

    def __complex__(self):
        return self.x + self.y * 1j

    def __abs__(self):
        return hypot(self.x, self.y)

    def __pow__(self, other):
        r_raised = abs(self) ** other
        argz_multiplied = self.argz() * other

        real_part = round(r_raised * cos(argz_multiplied))
        imag_part = round(r_raised * sin(argz_multiplied))
        return self.__class__(real_part, imag_part)

    def conjugate(self):
        return self.__class__(self.real, -self.imag)

    def argz(self):
        return atan(self.imag / self.real)

    @property
    def real(self):
        """Emulate svg.path use of complex numbers"""
        return self.x

    @property
    def imag(self):
        """Emulate svg.path use of complex numbers"""
        return self.y

    def matrix_transform(self, matrix):
        self *= matrix
        return self

    def move_towards(self, p2, amount=1):
        if not isinstance(p2, Point):
            p2 = Point(p2)
        self += amount * (p2 - self)

    def distance_to(self, p2):
        return abs(self - p2)

    def angle_to(self, p2):
        p = p2 - self
        return Angle.radians(atan2(p.y, p.x))

    def polar_to(self, angle, distance):
        q = Point.polar(self, angle, distance)
        self.x = q.x
        self.y = q.y
        return self

    def reflected_across(self, p):
        return p + (p - self)

    @staticmethod
    def orientation(p, q, r):
        """Determine the clockwise, linear, or counterclockwise orientation of the given points"""
        val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
        if val == 0:
            return 0
        elif val > 0:
            return 1
        else:
            return 2

    @staticmethod
    def convex_hull(pts):
        if len(pts) == 0:
            return
        points = sorted(set(pts), key=lambda p: p[0])
        first_point_on_hull = points[0]
        point_on_hull = first_point_on_hull
        while True:
            yield point_on_hull
            endpoint = point_on_hull
            for t in points:
                if (
                    point_on_hull is endpoint
                    or Point.orientation(point_on_hull, t, endpoint) == 2
                ):
                    endpoint = t
            point_on_hull = endpoint
            if first_point_on_hull is point_on_hull:
                break

    @staticmethod
    def distance(p1, p2):
        dx = p1[0] - p2[0]
        dy = p1[1] - p2[1]
        dx *= dx
        dy *= dy
        return sqrt(dx + dy)

    @staticmethod
    def polar(p1, angle, r):
        dx = cos(angle) * r
        dy = sin(angle) * r
        return Point(p1[0] + dx, p1[1] + dy)

    @staticmethod
    def angle(p1, p2):
        return Angle.radians(atan2(p2[1] - p1[1], p2[0] - p1[0]))

    @staticmethod
    def towards(p1, p2, amount):
        tx = amount * (p2[0] - p1[0]) + p1[0]
        ty = amount * (p2[1] - p1[1]) + p1[1]
        return Point(tx, ty)


class Angle(float):
    """CSS Angle defines as used in SVG/CSS"""

    def __repr__(self):
        return "Angle(%.12f)" % self

    def __copy__(self):
        return Angle(self)

    def __eq__(self, other):
        # Python 2
        c1 = abs((self % tau) - (other % tau)) <= 1e-11
        return c1

    def normalized(self):
        return Angle(self % tau)

    @classmethod
    def parse(cls, angle_string):
        if not isinstance(angle_string, str):
            return
        angle_string = angle_string.lower()
        if angle_string.endswith("deg"):
            return Angle.degrees(float(angle_string[:-3]))
        if angle_string.endswith("grad"):
            return Angle.gradians(float(angle_string[:-4]))
        if angle_string.endswith(
            "rad"
        ):  # Must be after 'grad' since 'grad' ends with 'rad' too.
            return Angle.radians(float(angle_string[:-3]))
        if angle_string.endswith("turn"):
            return Angle.turns(float(angle_string[:-4]))
        if angle_string.endswith("%"):
            return Angle.turns(float(angle_string[:-1]) / 100.0)
        return Angle.degrees(float(angle_string))

    @classmethod
    def radians(cls, radians):
        return cls(radians)

    @classmethod
    def degrees(cls, degrees):
        return cls(tau * degrees / 360.0)

    @classmethod
    def gradians(cls, gradians):
        return cls(tau * gradians / 400.0)

    @classmethod
    def turns(cls, turns):
        return cls(tau * turns)

    @property
    def as_radians(self):
        return self

    @property
    def as_degrees(self):
        return self * 360.0 / tau

    @property
    def as_positive_degrees(self):
        v = self.as_degrees
        while v < 0:
            v += 360.0
        return v

    @property
    def as_gradians(self):
        return self * 400.0 / tau

    @property
    def as_turns(self):
        return self / tau

    def is_orthogonal(self):
        return (self % (tau / 4.0)) == 0


class Matrix:
    """ "
    Provides svg matrix interfacing.

    SVG 7.15.3 defines the matrix form as:
    [a c  e]
    [b d  f]

    While e and f are defined as floats, they can be for limited periods defined as a Length.
    With regard to CSS, it's reasonable to perform operations like 'transform(20cm, 20cm)' and
    expect these to be treated consistently. Performing other matrix operations in a consistent
    way. However, render must be called to change these parameters into float locations prior to
    any operation which might be used to transform a point or polyline or path object.
    """

    def __init__(self, *components, **kwargs):
        self.a = 1.0
        self.b = 0.0
        self.c = 0.0
        self.d = 1.0
        self.e = 0.0
        self.f = 0.0
        len_args = len(components)
        if len_args == 0:
            pass
        elif len_args == 1:
            m = components[0]
            if isinstance(m, str):
                self.parse(m)
                self.render(**kwargs)
            else:
                self.a = m[0]
                self.b = m[1]
                self.c = m[2]
                self.d = m[3]
                self.e = m[4]
                self.f = m[5]
        else:
            self.a = components[0]
            self.b = components[1]
            self.c = components[2]
            self.d = components[3]
            self.e = components[4]
            self.f = components[5]
            self.render(**kwargs)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        if other is None:
            return False
        if isinstance(other, str):
            other = Matrix(other)
        if not isinstance(other, Matrix):
            return False
        if abs(self.a - other.a) > 1e-12:
            return False
        if abs(self.b - other.b) > 1e-12:
            return False
        if abs(self.c - other.c) > 1e-12:
            return False
        if abs(self.d - other.d) > 1e-12:
            return False
        if self.e != other.e and abs(self.e - other.e) > 1e-12:
            return False
        if self.f != other.f and abs(self.f - other.f) > 1e-12:
            return False
        return True

    def __len__(self):
        return 6

    def __invert__(self):
        m = self.__copy__()
        return m.inverse()

    def __matmul__(self, other):
        m = copy(self)
        m.__imatmul__(other)
        return m

    def __rmatmul__(self, other):
        m = copy(other)
        m.__imatmul__(self)
        return m

    def __imatmul__(self, other):
        if isinstance(other, str):
            other = Matrix(other)
        self.a, self.b, self.c, self.d, self.e, self.f = Matrix.matrix_multiply(
            self, other
        )
        return self

    __mul__ = __matmul__
    __rmul__ = __rmatmul__
    __imul__ = __imatmul__

    def __getitem__(self, item):
        if item == 0:
            return float(self.a)
        elif item == 1:
            return float(self.b)
        elif item == 2:
            return float(self.c)
        elif item == 3:
            return float(self.d)
        elif item == 4:
            return self.e
        elif item == 5:
            return self.f

    def __setitem__(self, key, value):
        if key == 0:
            self.a = value
        elif key == 1:
            self.b = value
        elif key == 2:
            self.c = value
        elif key == 3:
            self.d = value
        elif key == 4:
            self.e = value
        elif key == 5:
            self.f = value

    def __repr__(self):
        return "Matrix(%s, %s, %s, %s, %s, %s)" % (
            Length.str(self.a),
            Length.str(self.b),
            Length.str(self.c),
            Length.str(self.d),
            Length.str(self.e),
            Length.str(self.f),
        )

    def __copy__(self):
        return Matrix(self.a, self.b, self.c, self.d, self.e, self.f)

    def __str__(self):
        """
        Many of SVG's graphics operations utilize 2x3:

        :returns string representation of matrix.
        """
        return "[%3f, %3f,\n %3f, %3f,   %s, %s]" % (
            self.a,
            self.c,
            self.b,
            self.d,
            self.e,
            self.f,
        )

    def parse(self, transform_str):
        """Parses the svg transform string.

        Transforms from SVG 1.1 have a smaller complete set of operations. Whereas in SVG 2.0 they gain
        the CSS transforms and the additional functions and parsing that go with that. This parse is
        compatible with SVG 1.1 and the SVG 2.0 which includes the CSS 2d superset.

        CSS transforms have scalex() scaley() translatex(), translatey(), and skew() (deprecated).
        2D CSS angles haves units: "deg" tau / 360, "rad" tau/tau, "grad" tau/400, "turn" tau.
        2D CSS distances have length/percentages: "px", "cm", "mm", "in", "pt", etc. (+|-)?d+%

        In the case of percentages there must be a known height and width to properly create a matrix out of that.

        """
        if not transform_str:
            return
        if not isinstance(transform_str, str):
            raise TypeError("Must provide a string to parse")

        for sub_element in REGEX_TRANSFORM_TEMPLATE.findall(transform_str.lower()):
            name = sub_element[0]
            params = tuple(REGEX_TRANSFORM_PARAMETER.findall(sub_element[1]))
            params = [mag + units for mag, units in params]
            if SVG_TRANSFORM_MATRIX == name:
                params = map(float, params)
                self.pre_cat(*params)
            elif SVG_TRANSFORM_TRANSLATE == name:
                try:
                    x_param = Length(params[0]).value()
                except IndexError:
                    continue
                try:
                    y_param = Length(params[1]).value()
                    self.pre_translate(x_param, y_param)
                except IndexError:
                    self.pre_translate(x_param)
            elif SVG_TRANSFORM_TRANSLATE_X == name:
                self.pre_translate(Length(params[0]).value(), 0)
            elif SVG_TRANSFORM_TRANSLATE_Y == name:
                self.pre_translate(0, Length(params[0]).value())
            elif SVG_TRANSFORM_SCALE == name:
                params = map(float, params)
                self.pre_scale(*params)
            elif SVG_TRANSFORM_SCALE_X == name:
                self.pre_scale(float(params[0]), 1)
            elif SVG_TRANSFORM_SCALE_Y == name:
                self.pre_scale(1, float(params[0]))
            elif SVG_TRANSFORM_ROTATE == name:
                angle = Angle.parse(params[0])
                try:
                    x_param = Length(params[1]).value()
                except IndexError:
                    self.pre_rotate(angle)
                    continue
                try:
                    y_param = Length(params[2]).value()
                    self.pre_rotate(angle, x_param, y_param)
                except IndexError:
                    self.pre_rotate(angle, x_param)
            elif SVG_TRANSFORM_SKEW == name:
                angle_a = Angle.parse(params[0])
                try:
                    angle_b = Angle.parse(params[1])
                except IndexError:  # this isn't valid.
                    continue
                try:
                    x_param = Length(params[2]).value()
                except IndexError:
                    self.pre_skew(angle_a, angle_b)
                    continue
                try:
                    y_param = Length(params[3]).value()
                    self.pre_skew(angle_a, angle_b, x_param, y_param)
                except IndexError:
                    self.pre_skew(angle_a, angle_b, x_param)
            elif SVG_TRANSFORM_SKEW_X == name:
                angle_a = Angle.parse(params[0])
                try:
                    x_param = Length(params[1]).value()
                except IndexError:
                    self.pre_skew_x(angle_a)
                    continue
                try:
                    y_param = Length(params[2]).value()
                    self.pre_skew_x(angle_a, x_param, y_param)
                except IndexError:
                    self.pre_skew_x(angle_a, x_param)
            elif SVG_TRANSFORM_SKEW_Y == name:
                angle_b = Angle.parse(params[0])
                try:
                    x_param = Length(params[1]).value()
                except IndexError:
                    self.pre_skew_y(angle_b)
                    continue
                try:
                    y_param = Length(params[2]).value()
                    self.pre_skew_y(angle_b, x_param, y_param)
                except IndexError:
                    self.pre_skew_y(angle_b, x_param)
        return self

    def render(
        self,
        ppi=None,
        relative_length=None,
        width=None,
        height=None,
        font_size=None,
        font_height=None,
        viewbox=None,
        **kwargs
    ):
        """
        Provides values to turn trans_x and trans_y values into user units floats rather
        than Lengths by giving the required information to perform the conversions.
        """
        if isinstance(self.e, Length):
            if width is None and relative_length is not None:
                width = relative_length
            self.e = self.e.value(
                ppi=ppi,
                relative_length=width,
                font_size=font_size,
                font_height=font_height,
                viewbox=viewbox,
            )

        if isinstance(self.f, Length):
            if height is None and relative_length is not None:
                height = relative_length
            self.f = self.f.value(
                ppi=ppi,
                relative_length=height,
                font_size=font_size,
                font_height=font_height,
                viewbox=viewbox,
            )
        return self

    @property
    def determinant(self):
        return self.a * self.d - self.c * self.b

    def value_trans_x(self):
        return self.e

    def value_trans_y(self):
        return self.f

    def value_scale_x(self):
        return float(self.a)

    def value_scale_y(self):
        return float(self.d)

    def value_skew_x(self):
        return float(self.b)

    def value_skew_y(self):
        return float(self.c)

    def reset(self):
        """Resets matrix to identity."""
        self.a = 1.0
        self.b = 0.0
        self.c = 0.0
        self.d = 1.0

        self.e = 0.0
        self.f = 0.0

    def inverse(self):
        """
        SVG Matrix:
        [a c e]
        [b d f]
        """
        m00 = self.a
        m01 = self.c
        m02 = self.e
        m10 = self.b
        m11 = self.d
        m12 = self.f
        determinant = m00 * m11 - m01 * m10
        inverse_determinant = 1.0 / determinant
        self.a = m11 * inverse_determinant
        self.c = -m01 * inverse_determinant
        self.b = -m10 * inverse_determinant
        self.d = m00 * inverse_determinant

        self.e = (m01 * m12 - m02 * m11) * inverse_determinant
        self.f = (m10 * m02 - m00 * m12) * inverse_determinant
        return self

    def vector(self):
        """
        provide the matrix suitable for multiplying vectors. This will be the matrix with the same rotation and scale
        aspects but with no translation. This matrix is for multiplying vector elements where the position doesn't
        matter but the scaling and rotation do.
        :return:
        """
        return Matrix(self.a, self.b, self.c, self.d, 0.0, 0.0)

    def is_identity(self):
        return (
            self.a == 1
            and self.b == 0
            and self.c == 0
            and self.d == 1
            and self.e == 0
            and self.f == 0
        )

    def post_cat(self, *components):
        mx = Matrix(*components)
        self.__imatmul__(mx)

    def post_scale(self, sx=1.0, sy=None, x=0.0, y=0.0):
        if sy is None:
            sy = sx
        if x is None:
            x = 0.0
        if y is None:
            y = 0.0
        if x == 0 and y == 0:
            self.post_cat(Matrix.scale(sx, sy))
        else:
            self.post_translate(-x, -y)
            self.post_scale(sx, sy)
            self.post_translate(x, y)

    def post_scale_x(self, sx=1.0, x=0.0, y=0.0):
        self.post_scale(sx, 1, x, y)

    def post_scale_y(self, sy=1.0, x=0.0, y=0.0):
        self.post_scale(1, sy, x, y)

    def post_translate(self, tx=0.0, ty=0.0):
        self.post_cat(Matrix.translate(tx, ty))

    def post_translate_x(self, tx=0.0):
        self.post_translate(tx, 0.0)

    def post_translate_y(self, ty=0.0):
        self.post_translate(0.0, ty)

    def post_rotate(self, angle, x=0.0, y=0.0):
        if x is None:
            x = 0.0
        if y is None:
            y = 0.0
        if x == 0 and y == 0:
            self.post_cat(Matrix.rotate(angle))  # self %= self.get_rotate(theta)
        else:
            matrix = Matrix()
            matrix.post_translate(-x, -y)
            matrix.post_cat(Matrix.rotate(angle))
            matrix.post_translate(x, y)
            self.post_cat(matrix)

    def post_skew(self, angle_a=0.0, angle_b=0.0, x=0.0, y=0.0):
        if x is None:
            x = 0
        if y is None:
            y = 0
        if x == 0 and y == 0:
            self.post_cat(Matrix.skew(angle_a, angle_b))
        else:
            self.post_translate(-x, -y)
            self.post_skew(angle_a, angle_b)
            self.post_translate(x, y)

    def post_skew_x(self, angle_a=0.0, x=0.0, y=0.0):
        self.post_skew(angle_a, 0.0, x, y)

    def post_skew_y(self, angle_b=0.0, x=0.0, y=0.0):
        self.post_skew(0.0, angle_b, x, y)

    def pre_cat(self, *components):
        mx = Matrix(*components)
        self.a, self.b, self.c, self.d, self.e, self.f = Matrix.matrix_multiply(
            mx, self
        )

    def pre_scale(self, sx=1.0, sy=None, x=0.0, y=0.0):
        if sy is None:
            sy = sx
        if x is None:
            x = 0.0
        if y is None:
            y = 0.0
        if x == 0 and y == 0:
            self.pre_cat(Matrix.scale(sx, sy))
        else:
            self.pre_translate(x, y)
            self.pre_scale(sx, sy)
            self.pre_translate(-x, -y)

    def pre_scale_x(self, sx=1.0, x=0.0, y=0.0):
        self.pre_scale(sx, 1, x, y)

    def pre_scale_y(self, sy=1.0, x=0.0, y=0.0):
        self.pre_scale(1, sy, x, y)

    def pre_translate(self, tx=0.0, ty=0.0):
        self.pre_cat(Matrix.translate(tx, ty))

    def pre_translate_x(self, tx=0.0):
        self.pre_translate(tx, 0.0)

    def pre_translate_y(self, ty=0.0):
        self.pre_translate(0.0, ty)

    def pre_rotate(self, angle, x=0.0, y=0.0):
        if x is None:
            x = 0
        if y is None:
            y = 0
        if x == 0 and y == 0:
            self.pre_cat(Matrix.rotate(angle))
        else:
            self.pre_translate(x, y)
            self.pre_rotate(angle)
            self.pre_translate(-x, -y)

    def pre_skew(self, angle_a=0.0, angle_b=0.0, x=0.0, y=0.0):
        if x is None:
            x = 0
        if y is None:
            y = 0
        if x == 0 and y == 0:
            self.pre_cat(Matrix.skew(angle_a, angle_b))
        else:
            self.pre_translate(x, y)
            self.pre_skew(angle_a, angle_b)
            self.pre_translate(-x, -y)

    def pre_skew_x(self, angle_a=0.0, x=0.0, y=0.0):
        self.pre_skew(angle_a, 0, x, y)

    def pre_skew_y(self, angle_b=0.0, x=0.0, y=0.0):
        self.pre_skew(0.0, angle_b, x, y)

    def point_in_inverse_space(self, v0):
        inverse = Matrix(self)
        inverse.inverse()
        return inverse.point_in_matrix_space(v0)

    def point_in_matrix_space(self, v0):
        return Point(
            v0[0] * self.a + v0[1] * self.c + 1 * self.e,
            v0[0] * self.b + v0[1] * self.d + 1 * self.f,
        )

    def transform_point(self, v):
        nx = v[0] * self.a + v[1] * self.c + 1 * self.e
        ny = v[0] * self.b + v[1] * self.d + 1 * self.f
        v[0] = nx
        v[1] = ny
        return v

    def transform_vector(self, v):
        """
        Applies the transformation without the translation.
        """
        nx = v[0] * self.a + v[1] * self.c
        ny = v[0] * self.b + v[1] * self.d
        v[0] = nx
        v[1] = ny
        return v

    @classmethod
    def scale(cls, sx=1.0, sy=None):
        if sy is None:
            sy = sx
        return cls(sx, 0, 0, sy, 0, 0)

    @classmethod
    def scale_x(cls, sx=1.0):
        return cls.scale(sx, 1.0)

    @classmethod
    def scale_y(cls, sy=1.0):
        return cls.scale(1.0, sy)

    @classmethod
    def translate(cls, tx=0.0, ty=0.0):
        """SVG Matrix:
        [a c e]
        [b d f]
        """
        return cls(1.0, 0.0, 0.0, 1.0, tx, ty)

    @classmethod
    def translate_x(cls, tx=0.0):
        return cls.translate(tx, 0)

    @classmethod
    def translate_y(cls, ty=0.0):
        return cls.translate(0.0, ty)

    @classmethod
    def rotate(cls, angle=0.0):
        ct = cos(angle)
        st = sin(angle)
        return cls(ct, st, -st, ct, 0.0, 0.0)

    @classmethod
    def skew(cls, angle_a=0.0, angle_b=0.0):
        aa = tan(angle_a)
        bb = tan(angle_b)
        return cls(1.0, bb, aa, 1.0, 0.0, 0.0)

    @classmethod
    def skew_x(cls, angle=0.0):
        return cls.skew(angle, 0.0)

    @classmethod
    def skew_y(cls, angle=0.0):
        return cls.skew(0.0, angle)

    @classmethod
    def identity(cls):
        """
        1, 0, 0,
        0, 1, 0,
        """
        return cls()

    @staticmethod
    def matrix_multiply(m, s):
        """
        [a c e]      [a c e]   [a b 0]
        [b d f]   %  [b d f] = [c d 0]
        [0 0 1]      [0 0 1]   [e f 1]

        :param m: matrix operand
        :param s: matrix operand
        :return: multiplied matrix.
        """
        r0 = (
            s.a * m.a + s.c * m.b + s.e * 0,
            s.a * m.c + s.c * m.d + s.e * 0,
            s.a * m.e + s.c * m.f + s.e * 1,
        )

        r1 = (
            s.b * m.a + s.d * m.b + s.f * 0,
            s.b * m.c + s.d * m.d + s.f * 0,
            s.b * m.e + s.d * m.f + s.f * 1,
        )
        return float(r0[0]), float(r1[0]), float(r0[1]), float(r1[1]), r0[2], r1[2]


class Viewbox:
    def __init__(self, viewbox, preserve_aspect_ratio=None):
        """
        Viewbox controls the scaling between the drawing size view that is observing that drawing.

        :param viewbox: either values or viewbox attribute or a Viewbox object
        :param preserve_aspect_ratio: preserveAspectRatio
        """
        self.x = None
        self.y = None
        self.width = None
        self.height = None
        self.preserve_aspect_ratio = preserve_aspect_ratio
        if isinstance(viewbox, dict):
            self.property_by_values(viewbox)
        elif isinstance(viewbox, Viewbox):
            self.property_by_object(viewbox)
        else:
            self.set_viewbox(viewbox)

    def __str__(self):
        return "%s %s %s %s" % (
            Length.str(self.x),
            Length.str(self.y),
            Length.str(self.width),
            Length.str(self.height),
        )

    def property_by_object(self, obj):
        self.x = obj.x
        self.y = obj.y
        self.width = obj.width
        self.height = obj.height
        self.preserve_aspect_ratio = obj.preserve_aspect_ratio

    def property_by_values(self, values):
        viewbox = values.get(SVG_ATTR_VIEWBOX)
        if viewbox is not None:
            self.set_viewbox(viewbox)
        if SVG_ATTR_PRESERVEASPECTRATIO in values:
            self.preserve_aspect_ratio = values[SVG_ATTR_PRESERVEASPECTRATIO]

    def set_viewbox(self, viewbox):
        if viewbox is not None:
            dims = list(REGEX_FLOAT.findall(viewbox))
            try:
                self.x = float(dims[0])
                self.y = float(dims[1])
                self.width = float(dims[2])
                self.height = float(dims[3])
            except IndexError:
                pass

    def transform(self, element):
        return Viewbox.viewbox_transform(
            element.x,
            element.y,
            element.width,
            element.height,
            self.x,
            self.y,
            self.width,
            self.height,
            self.preserve_aspect_ratio,
        )

    @staticmethod
    def viewbox_transform(
        e_x, e_y, e_width, e_height, vb_x, vb_y, vb_width, vb_height, aspect
    ):
        """
        SVG 1.1 7.2, SVG 2.0 8.2 equivalent transform of an SVG viewport.
        With regards to https://github.com/w3c/svgwg/issues/215 use 8.2 version.

        It creates transform commands equal to that viewport expected.

        Let e-x, e-y, e-width, e-height be the position and size of the element respectively.
        Let vb-x, vb-y, vb-width, vb-height be the min-x, min-y, width and height values of the viewBox attribute
        respectively.

        Let align be the align value of preserveAspectRatio, or 'xMidYMid' if preserveAspectRatio is not defined.
        Let meetOrSlice be the meetOrSlice value of preserveAspectRatio, or 'meet' if preserveAspectRatio is not defined
        or if meetOrSlice is missing from this value.

        :param e_x: element_x value
        :param e_y: element_y value
        :param e_width: element_width value
        :param e_height: element_height value
        :param vb_x: viewbox_x value
        :param vb_y: viewbox_y value
        :param vb_width: viewbox_width value
        :param vb_height: viewbox_height value
        :param aspect: preserve aspect ratio value
        :return: string of the SVG transform commands to account for the viewbox.
        """
        if (
            e_x is None
            or e_y is None
            or e_width is None
            or e_height is None
            or vb_x is None
            or vb_y is None
            or vb_width is None
            or vb_height is None
        ):
            return ""
        if aspect is not None:
            aspect_slice = aspect.split(" ")
            try:
                align = aspect_slice[0]
            except IndexError:
                align = "xMidyMid"
            try:
                meet_or_slice = aspect_slice[1]
            except IndexError:
                meet_or_slice = "meet"
        else:
            align = "xMidyMid"
            meet_or_slice = "meet"
        # Initialize scale-x to e-width/vb-width.
        scale_x = e_width / vb_width
        # Initialize scale-y to e-height/vb-height.
        scale_y = e_height / vb_height

        # If align is not 'none' and meetOrSlice is 'meet', set the larger of scale-x and scale-y to the smaller.
        if align != SVG_VALUE_NONE and meet_or_slice == "meet":
            scale_x = scale_y = min(scale_x, scale_y)
        # Otherwise, if align is not 'none' and meetOrSlice is 'slice', set the smaller of scale-x and scale-y to the larger
        elif align != SVG_VALUE_NONE and meet_or_slice == "slice":
            scale_x = scale_y = max(scale_x, scale_y)
        # Initialize translate-x to e-x - (vb-x * scale-x).
        translate_x = e_x - (vb_x * scale_x)
        # Initialize translate-y to e-y - (vb-y * scale-y)
        translate_y = e_y - (vb_y * scale_y)
        # If align contains 'xMid', add (e-width - vb-width * scale-x) / 2 to translate-x.
        align = align.lower()
        if "xmid" in align:
            translate_x += (e_width - vb_width * scale_x) / 2.0
        # If align contains 'xMax', add (e-width - vb-width * scale-x) to translate-x.
        if "xmax" in align:
            translate_x += e_width - vb_width * scale_x
        # If align contains 'yMid', add (e-height - vb-height * scale-y) / 2 to translate-y.
        if "ymid" in align:
            translate_y += (e_height - vb_height * scale_y) / 2.0
        # If align contains 'yMax', add (e-height - vb-height * scale-y) to translate-y.
        if "ymax" in align:
            translate_y += e_height - vb_height * scale_y
        # The transform applied to content contained by the element is given by:
        # translate(translate-x, translate-y) scale(scale-x, scale-y)
        if isinstance(scale_x, Length) or isinstance(scale_y, Length):
            raise ValueError
        if translate_x == 0 and translate_y == 0:
            if scale_x == 1 and scale_y == 1:
                return ""  # Nothing happens.
            else:
                return "scale(%s, %s)" % (Length.str(scale_x), Length.str(scale_y))
        else:
            if scale_x == 1 and scale_y == 1:
                return "translate(%s, %s)" % (
                    Length.str(translate_x),
                    Length.str(translate_y),
                )
            else:
                return "translate(%s, %s) scale(%s, %s)" % (
                    Length.str(translate_x),
                    Length.str(translate_y),
                    Length.str(scale_x),
                    Length.str(scale_y),
                )
