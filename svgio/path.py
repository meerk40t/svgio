from .constants import *
from .basic import *

try:
    from math import tau
except ImportError:
    from math import pi

    tau = pi * 2

try:
    from collections.abc import MutableSequence  # noqa
except ImportError:
    from collections import MutableSequence  # noqa


########################
#  Path Lexical Parsing
########################


class PathLexicalParser:
    def __init__(self):
        self.parser = None
        self.pathd = None
        self.pos = 0
        self.limit = 0
        self.inline_close = None

    def _command(self):
        while self.pos < self.limit:
            match = svg_re.match(self.pathd, self.pos)
            if match is None:
                return None  # Did not match at command sequence.
            self.pos = match.end()
            kind = match.lastgroup
            if kind == "SKIP":
                continue
            return match.group()
        return None

    def _more(self):
        while self.pos < self.limit:
            match = num_re.match(self.pathd, self.pos)
            if match is None:
                return False
            kind = match.lastgroup
            if kind == "CLOSE":
                self.inline_close = match.group()
                return False
            if kind == "SKIP":
                # move skipped elements forward.
                self.pos = match.end()
                continue
            return True
        return None

    def _number(self):
        while self.pos < self.limit:
            match = num_re.match(self.pathd, self.pos)
            if match is None:
                break  # No more matches.
            kind = match.lastgroup
            if kind == "CLOSE":
                # Inline Close
                self.inline_close = match.group()
                return None
            self.pos = match.end()
            if kind == "SKIP":
                continue
            return float(match.group())
        return None

    def _flag(self):
        while self.pos < self.limit:
            match = flag_re.match(self.pathd, self.pos)
            if match is None:
                break  # No more matches.
            self.pos = match.end()
            kind = match.lastgroup
            if kind == "SKIP":
                continue
            return bool(int(match.group()))
        return None

    def _coord(self):
        x = self._number()
        if x is None:
            return None
        y = self._number()
        if y is None:
            raise ValueError
        return x, y

    def _rcoord(self):
        position = self._coord()
        if position is None:
            return None
        current_pos = self.parser.current_point
        if current_pos is None:
            return position
        return position[0] + current_pos.x, position[1] + current_pos.y

    def parse(self, parser, pathd):
        self.parser = parser
        self.parser.start()
        self.pathd = pathd
        self.pos = 0
        self.limit = len(pathd)
        while True:
            cmd = self._command()
            if cmd is None:
                return
            elif cmd == "z" or cmd == "Z":
                if self._more():
                    raise ValueError
                self.parser.closed(relative=cmd.islower())
                self.inline_close = None
                continue
            elif cmd == "m":
                if not self._more():
                    raise ValueError
                coord = self._rcoord()
                self.parser.move(coord, relative=True)
                while self._more():
                    coord = self._rcoord()
                    self.parser.line(coord, relative=True)
            elif cmd == "M":
                if not self._more():
                    raise ValueError
                coord = self._coord()
                self.parser.move(coord, relative=False)
                while self._more():
                    coord = self._coord()
                    self.parser.line(coord, relative=False)
            elif cmd == "l":
                while True:
                    coord = self._rcoord()
                    if coord is None:
                        coord = self.inline_close
                        if coord is None:
                            raise ValueError
                    self.parser.line(coord, relative=True)
                    if not self._more():
                        break
            elif cmd == "L":
                while True:
                    coord = self._coord()
                    if coord is None:
                        coord = self.inline_close
                        if coord is None:
                            raise ValueError
                    self.parser.line(coord, relative=False)
                    if not self._more():
                        break
            elif cmd == "t":
                while True:
                    coord = self._rcoord()
                    if coord is None:
                        coord = self.inline_close
                        if coord is None:
                            raise ValueError
                    self.parser.smooth_quad(coord, relative=True)
                    if not self._more():
                        break
            elif cmd == "T":
                while True:
                    coord = self._coord()
                    if coord is None:
                        coord = self.inline_close
                        if coord is None:
                            raise ValueError
                    self.parser.smooth_quad(coord, relative=False)
                    if not self._more():
                        break
            elif cmd == "h":
                while True:
                    value = self._number()
                    self.parser.horizontal(value, relative=True)
                    if not self._more():
                        break
            elif cmd == "H":
                while True:
                    value = self._number()
                    self.parser.horizontal(value, relative=False)
                    if not self._more():
                        break
            elif cmd == "v":
                while True:
                    value = self._number()
                    self.parser.vertical(value, relative=True)
                    if not self._more():
                        break
            elif cmd == "V":
                while self._more():
                    value = self._number()
                    self.parser.vertical(value, relative=False)
            elif cmd == "c":
                while True:
                    coord1, coord2, coord3 = (
                        self._rcoord(),
                        self._rcoord(),
                        self._rcoord(),
                    )
                    if coord1 is None:
                        coord1 = self.inline_close
                        if coord1 is None:
                            raise ValueError
                    if coord2 is None:
                        coord2 = self.inline_close
                        if coord2 is None:
                            raise ValueError
                    if coord3 is None:
                        coord3 = self.inline_close
                        if coord3 is None:
                            raise ValueError
                    self.parser.cubic(coord1, coord2, coord3, relative=True)
                    if not self._more():
                        break
            elif cmd == "C":
                while True:
                    coord1, coord2, coord3 = self._coord(), self._coord(), self._coord()
                    if coord1 is None:
                        coord1 = self.inline_close
                        if coord1 is None:
                            raise ValueError
                    if coord2 is None:
                        coord2 = self.inline_close
                        if coord2 is None:
                            raise ValueError
                    if coord3 is None:
                        coord3 = self.inline_close
                        if coord3 is None:
                            raise ValueError
                    self.parser.cubic(coord1, coord2, coord3, relative=False)
                    if not self._more():
                        break
            elif cmd == "q":
                while True:
                    coord1, coord2 = self._rcoord(), self._rcoord()
                    if coord1 is None:
                        coord1 = self.inline_close
                        if coord1 is None:
                            raise ValueError
                    if coord2 is None:
                        coord2 = self.inline_close
                        if coord2 is None:
                            raise ValueError
                    self.parser.quad(coord1, coord2, relative=True)
                    if not self._more():
                        break
            elif cmd == "Q":
                while True:
                    coord1, coord2 = self._coord(), self._coord()
                    if coord1 is None:
                        coord1 = self.inline_close
                        if coord1 is None:
                            raise ValueError
                    if coord2 is None:
                        coord2 = self.inline_close
                        if coord2 is None:
                            raise ValueError
                    self.parser.quad(coord1, coord2, relative=False)
                    if not self._more():
                        break
            elif cmd == "s":
                while True:
                    coord1, coord2 = self._rcoord(), self._rcoord()
                    if coord1 is None:
                        coord1 = self.inline_close
                        if coord1 is None:
                            raise ValueError
                    if coord2 is None:
                        coord2 = self.inline_close
                        if coord2 is None:
                            raise ValueError
                    self.parser.smooth_cubic(coord1, coord2, relative=True)
                    if not self._more():
                        break
            elif cmd == "S":
                while True:
                    coord1, coord2 = self._coord(), self._coord()
                    if coord1 is None:
                        coord1 = self.inline_close
                        if coord1 is None:
                            raise ValueError
                    if coord2 is None:
                        coord2 = self.inline_close
                        if coord2 is None:
                            raise ValueError
                    self.parser.smooth_cubic(coord1, coord2, relative=False)
                    if not self._more():
                        break
            elif cmd == "a":
                while self._more():
                    rx, ry, rotation, arc, sweep, coord = (
                        self._number(),
                        self._number(),
                        self._number(),
                        self._flag(),
                        self._flag(),
                        self._rcoord(),
                    )
                    if sweep is None:
                        raise ValueError
                    if coord is None:
                        coord = self.inline_close
                        if coord is None:
                            raise ValueError
                    self.parser.arc(rx, ry, rotation, arc, sweep, coord, relative=True)
            elif cmd == "A":
                while self._more():
                    rx, ry, rotation, arc, sweep, coord = (
                        self._number(),
                        self._number(),
                        self._number(),
                        self._flag(),
                        self._flag(),
                        self._coord(),
                    )
                    if coord is None:
                        coord = self.inline_close
                        if coord is None:
                            raise ValueError
                    self.parser.arc(rx, ry, rotation, arc, sweep, coord, relative=False)
        self.parser.end()


########################
#  SVG Path Segments
########################


class PathSegment:
    """
    Path Segments are the base class for all the segment within a Path.
    These are defined in SVG 1.1 8.3 and SVG 2.0 9.3
    https://www.w3.org/TR/SVG11/paths.html#PathData
    https://www.w3.org/TR/SVG2/paths.html#PathElement

    These segments define a 1:1 relationship with the path_d or path data attribute, denoted in
    SVG by the 'd' attribute. These are moveto, closepath, lineto, and the curves which are cubic
    bezier curves, quadratic bezier curves, and elliptical arc. These are classed as Move, Close,
    Line, CubicBezier, QuadraticBezier, and Arc. And in path_d are denoted as M, Z, L, C, Q, A.

    There are lowercase versions of these commands. And for C, and Q there are S and T which are
    smooth versions. For lines there are also V and H commands which denote vertical and horizontal
    versions of the line command.

    The major difference between paths in 1.1 and 2.0 is the use of Z to truncate a command to close.
    "M0,0C 0,100 100,0 z is valid in 2.0 since the last z replaces the 0,0. These are read by
    svg.elements but they are not written.
    """

    def __init__(self, **kwargs):
        try:
            self.relative = bool(kwargs["relative"])
        except (KeyError, ValueError):
            self.relative = False
        try:
            self.smooth = bool(kwargs["smooth"])
        except (KeyError, ValueError):
            self.smooth = True
        self.start = None
        self.end = None

    def __mul__(self, other):
        if isinstance(other, (Matrix, str)):
            n = copy(self)
            n *= other
            return n
        return NotImplemented

    __rmul__ = __mul__

    def __iadd__(self, other):
        if isinstance(other, PathSegment):
            path = Path(self, other)
            return path
        elif isinstance(other, str):
            path = Path(self) + other
            return path
        return NotImplemented

    __add__ = __iadd__

    def __str__(self):
        """
        This defines an individual path segment string. Since this isn't part of a Path it appends a pseudo-Move
        command to correctly provide the starting position.
        :return: string representation of the object.
        """
        d = self.d()
        if self.start is not None:
            if self.relative:
                return "m %s %s" % (self.start, d)
            else:
                return "M %s %s" % (self.start, d)
        return d

    def __iter__(self):
        self.n = -1
        return self

    def __next__(self):
        self.n += 1
        try:
            val = self[self.n]
            if val is None:
                self.n += 1
                val = self[self.n]
            return val
        except IndexError:
            raise StopIteration

    next = __next__

    @staticmethod
    def segment_length(
        curve,
        start=0.0,
        end=1.0,
        start_point=None,
        end_point=None,
        error=ERROR,
        min_depth=MIN_DEPTH,
        depth=0,
    ):
        """Recursively approximates the length by straight lines"""
        if start_point is None:
            start_point = curve.point(start)
        if end_point is None:
            end_point = curve.point(end)
        mid = (start + end) / 2.0
        mid_point = curve.point(mid)
        length = abs(end_point - start_point)
        first_half = abs(mid_point - start_point)
        second_half = abs(end_point - mid_point)

        length2 = first_half + second_half
        if (length2 - length > error) or (depth < min_depth):
            # Calculate the length of each segment:
            depth += 1
            return PathSegment.segment_length(
                curve, start, mid, start_point, mid_point, error, min_depth, depth
            ) + PathSegment.segment_length(
                curve, mid, end, mid_point, end_point, error, min_depth, depth
            )
        # This is accurate enough.
        return length2

    def _line_length(self, start=0.0, end=1.0, error=ERROR, min_depth=MIN_DEPTH):
        return PathSegment.segment_length(
            self, start, end, error=error, min_depth=min_depth
        )

    def bbox(self):
        """returns the bounding box for the segment.
        xmin, ymin, xmax, ymax
        """
        xs = [p.x for p in self if p is not None]
        ys = [p.y for p in self if p is not None]
        xmin = min(xs)
        xmax = max(xs)
        ymin = min(ys)
        ymax = max(ys)
        return xmin, ymin, xmax, ymax

    def reverse(self):
        """
        Reverses the current path segment.
        """
        end = self.end
        self.end = self.start
        self.start = end

    def point(self, position):
        """
        Returns the point at a given amount through the path segment.
        :param position:  t value between 0 and 1
        :return: Point instance
        """
        return Point(self.npoint([position])[0])

    def npoint(self, positions):
        """
        Returns the points at given positions along the path segment
        :param positions: N-sized sequence of t value between 0 and 1
        :return: N-sized sequence of 2-sized sequence of float
        """
        return [self.end] * len(positions)

    def length(self, error=ERROR, min_depth=MIN_DEPTH):
        """
        Returns the length of this path segment.

        :param error:
        :param min_depth:
        :return:
        """
        return 0

    def d(self, current_point=None, relative=None, smooth=None):
        """Returns the fragment path_d value for the current path segment.

        For a relative segment the current_point must be provided. If it is omitted then only an absolute segment
        can be returned."""
        raise NotImplementedError


class Move(PathSegment):
    """Represents move commands. Moves to a new location without any path distance.
    Paths that consist of only move commands, are valid.

    Move serve to make discontinuous paths into continuous linked paths segments
    with non-drawn sections.
    """

    def __init__(self, *args, **kwargs):
        """
        Move commands most importantly go to a place. So if one location is given, that's the end point.
        If two locations are given then first is the start location.

        For many Move commands it is not necessary to have an original start location. The start point provides a
        linked locations for some elements that may require it. If known it can be provided.

        Move(p) where p is the End point.
        Move(s,e) where s is the Start point, e is the End point.
        Move(p, start=s) where p is End point, s is the Start point.
        Move(p, end=e) where p is the Start point, e is the End point.
        Move(start=s, end=e) where s is the Start point, e is the End point.
        """
        PathSegment.__init__(self, **kwargs)
        self.end = None
        self.start = None
        if len(args) == 0:
            if "end" in kwargs:
                self.end = kwargs["end"]
            if "start" in kwargs:
                self.start = kwargs["start"]
        elif len(args) == 1:
            if len(kwargs) == 0:
                self.end = args[0]
            else:
                if "end" in kwargs:
                    self.start = args[0]
                    self.end = kwargs["end"]
                elif "start" in kwargs:
                    self.start = kwargs["start"]
                    self.end = args[0]
        elif len(args) == 2:
            self.start = args[0]
            self.end = args[1]
        if self.start is not None:
            self.start = Point(self.start)
        if self.end is not None:
            self.end = Point(self.end)

    def __imul__(self, other):
        if isinstance(other, str):
            other = Matrix(other)
        if isinstance(other, Matrix):
            if self.start is not None:
                self.start *= other
            if self.end is not None:
                self.end *= other
        return self

    def __repr__(self):
        if self.start is None:
            return "Move(end=%s)" % repr(self.end)
        else:
            return "Move(start=%s, end=%s)" % (repr(self.start), repr(self.end))

    def __copy__(self):
        return Move(self.start, self.end, relative=self.relative)

    def __eq__(self, other):
        if not isinstance(other, Move):
            return NotImplemented
        return self.start == other.start and self.end == other.end

    def __ne__(self, other):
        if not isinstance(other, Move):
            return NotImplemented
        return not self == other

    def __len__(self):
        return 2

    def __getitem__(self, item):
        if item == 0:
            return self.start
        elif item == 1:
            return self.end
        else:
            raise IndexError

    def d(self, current_point=None, relative=None, smooth=None):
        if (
            current_point is None
            or (relative is None and self.relative)
            or (relative is not None and not relative)
        ):
            return "M %s" % self.end
        return "m %s" % (self.end - current_point)


class Curve(PathSegment):
    """Represents curve commands"""

    def __init__(self, start=None, end=None, **kwargs):
        PathSegment.__init__(self, **kwargs)
        self.start = Point(start) if start is not None else None
        self.end = Point(end) if end is not None else None


class Linear(PathSegment):
    """Represents line commands."""

    def __init__(self, start=None, end=None, **kwargs):
        PathSegment.__init__(self, **kwargs)
        self.start = Point(start) if start is not None else None
        self.end = Point(end) if end is not None else None

    def __copy__(self):
        return self.__class__(self.start, self.end, relative=self.relative)

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self.start == other.start and self.end == other.end

    def __ne__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return not self == other

    def __imul__(self, other):
        if isinstance(other, str):
            other = Matrix(other)
        if isinstance(other, Matrix):
            if self.start is not None:
                self.start *= other
            if self.end is not None:
                self.end *= other
        return self

    def __len__(self):
        return 2

    def __getitem__(self, item):
        if item == 0:
            return self.start
        elif item == 1:
            return self.end
        else:
            raise IndexError

    def npoint(self, positions):
        try:
            import numpy as np

            xy = np.empty(shape=(len(positions), 2), dtype=float)
            xy[:, 0] = np.interp(positions, [0, 1], [self.start.x, self.end.x])
            xy[:, 1] = np.interp(positions, [0, 1], [self.start.y, self.end.y])
            return xy
        except ImportError:
            return [Point.towards(self.start, self.end, pos) for pos in positions]

    def length(self, error=None, min_depth=None):
        if self.start is not None and self.end is not None:
            return Point.distance(self.end, self.start)
        else:
            return 0

    def closest_segment_point(self, p, respect_bounds=True):
        """ Gives the point on the line closest to the given point. """
        a = self.start
        b = self.end
        v_ap_x = p[0] - a.x
        v_ap_y = p[1] - a.y
        v_ab_x = b.x - a.x
        v_ab_y = b.y - a.y
        sq_distance_ab = v_ab_x * v_ab_x + v_ab_y * v_ab_y
        ab_ap_product = v_ab_x * v_ap_x + v_ab_y * v_ap_y
        if sq_distance_ab == 0:
            return 0  # Line is point.
        amount = ab_ap_product / float(sq_distance_ab)
        if respect_bounds:
            if amount > 1:
                amount = 1
            if amount < 0:
                amount = 0
        return self.point(amount)

    def d(self, current_point=None, relative=None, smooth=None):
        raise NotImplementedError


class Close(Linear):
    """Represents close commands. If this exists at the end of the shape then the shape is closed.
    the methodology of a single flag close fails in a couple ways. You can have multi-part shapes
    which can close or not close several times.
    """

    def __repr__(self):
        if self.start is None and self.end is None:
            return "Close()"
        s = self.start
        if s is not None:
            s = repr(s)
        e = self.end
        if e is not None:
            e = repr(e)
        return "Close(start=%s, end=%s)" % (s, e)

    def d(self, current_point=None, relative=None, smooth=None):
        if (
            current_point is None
            or (relative is None and self.relative)
            or (relative is not None and not relative)
        ):
            return "Z"
        else:
            return "z"


class Line(Linear):
    """Represents line commands."""

    def __repr__(self):
        if self.start is None:
            return "Line(end=%s)" % (repr(self.end))
        return "Line(start=%s, end=%s)" % (repr(self.start), repr(self.end))

    def d(self, current_point=None, relative=None, smooth=None):
        if (
            current_point is None
            or (relative is None and self.relative)
            or (relative is not None and not relative)
        ):
            return "L %s" % self.end
        else:
            return "l %s" % (self.end - current_point)


class QuadraticBezier(Curve):
    """Represents Quadratic Bezier commands."""

    def __init__(self, start, control, end, **kwargs):
        Curve.__init__(self, start, end, **kwargs)
        self.control = Point(control) if control is not None else None

    def __repr__(self):
        return "QuadraticBezier(start=%s, control=%s, end=%s)" % (
            repr(self.start),
            repr(self.control),
            repr(self.end),
        )

    def __copy__(self):
        return QuadraticBezier(
            self.start,
            self.control,
            self.end,
            relative=self.relative,
            smooth=self.smooth,
        )

    def __eq__(self, other):
        if not isinstance(other, QuadraticBezier):
            return NotImplemented
        return (
            self.start == other.start
            and self.end == other.end
            and self.control == other.control
        )

    def __ne__(self, other):
        if not isinstance(other, QuadraticBezier):
            return NotImplemented
        return not self == other

    def __imul__(self, other):
        if isinstance(other, str):
            other = Matrix(other)
        if isinstance(other, Matrix):
            if self.start is not None:
                self.start *= other
            if self.control is not None:
                self.control *= other
            if self.end is not None:
                self.end *= other
        return self

    def __len__(self):
        return 3

    def __getitem__(self, item):
        if item == 0:
            return self.start
        elif item == 1:
            return self.control
        elif item == 2:
            return self.end
        raise IndexError

    def npoint(self, positions):
        """Calculate the x,y position at a certain position of the path. `pos` may be a
        float or a NumPy array."""
        x0, y0 = self.start
        x1, y1 = self.control
        x2, y2 = self.end

        def _compute_point(position):
            # compute factors
            n_pos = 1 - position
            pos_2 = position ** 2
            n_pos_2 = n_pos ** 2
            n_pos_pos = n_pos * position

            return (
                n_pos_2 * x0 + 2 * n_pos_pos * x1 + pos_2 * x2,
                n_pos_2 * y0 + 2 * n_pos_pos * y1 + pos_2 * y2,
            )

        try:
            import numpy as np

            xy = np.empty(shape=(len(positions), 2))
            xy[:, 0], xy[:, 1] = _compute_point(np.array(positions))
            return xy
        except ImportError:
            return [Point(*_compute_point(position)) for position in positions]

    def bbox(self):
        """
        Returns the bounding box for the quadratic bezier curve.
        """
        n = self.start.x - self.control.x
        d = self.start.x - 2 * self.control.x + self.end.x
        if d != 0:
            t = n / float(d)
        else:
            t = 0.5
        if 0 < t < 1:
            x_values = [self.start.x, self.end.x, self.point(t).x]
        else:
            x_values = [self.start.x, self.end.x]
        n = self.start.y - self.control.y
        d = self.start.y - 2 * self.control.y + self.end.y
        if d != 0:
            t = n / float(d)
        else:
            t = 0.5
        if 0 < t < 1:
            y_values = [self.start.y, self.end.y, self.point(t).y]
        else:
            y_values = [self.start.y, self.end.y]
        return min(x_values), min(y_values), max(x_values), max(y_values)

    def length(self, error=None, min_depth=None):
        """Calculate the length of the path up to a certain position"""
        a = self.start - 2 * self.control + self.end
        b = 2 * (self.control - self.start)
        try:
            # For an explanation of this case, see
            # http://www.malczak.info/blog/quadratic-bezier-curve-length/
            A = 4 * (a.real ** 2 + a.imag ** 2)
            B = 4 * (a.real * b.real + a.imag * b.imag)
            C = b.real ** 2 + b.imag ** 2

            Sabc = 2 * sqrt(A + B + C)
            A2 = sqrt(A)
            A32 = 2 * A * A2
            C2 = 2 * sqrt(C)
            BA = B / A2

            s = (
                A32 * Sabc
                + A2 * B * (Sabc - C2)
                + (4 * C * A - B ** 2) * log((2 * A2 + BA + Sabc) / (BA + C2))
            ) / (4 * A32)
        except (ZeroDivisionError, ValueError):
            # a_dot_b = a.real * b.real + a.imag * b.imag
            if abs(a) < 1e-10:
                s = abs(b)
            else:
                k = abs(b) / abs(a)
                if k >= 2:
                    s = abs(b) - abs(a)
                else:
                    s = abs(a) * (k ** 2 / 2 - k + 1)
        return s

    def is_smooth_from(self, previous):
        """Checks if this segment would be a smooth segment following the previous"""
        if isinstance(previous, QuadraticBezier):
            return self.start == previous.end and (self.control - self.start) == (
                previous.end - previous.control
            )
        else:
            return self.control == self.start

    def d(self, current_point=None, relative=None, smooth=None):
        if (smooth is None and self.smooth) or (smooth is not None and smooth):
            if (
                current_point is None
                or (relative is None and self.relative)
                or (relative is not None and not relative)
            ):
                return "T %s" % self.end
            else:
                return "t %s" % (self.end - current_point)
        else:
            if (
                current_point is None
                or (relative is None and self.relative)
                or (relative is not None and not relative)
            ):
                return "Q %s %s" % (self.control, self.end)
            else:
                return "q %s %s" % (
                    self.control - current_point,
                    self.end - current_point,
                )


class CubicBezier(Curve):
    """Represents Cubic Bezier commands."""

    def __init__(self, start, control1, control2, end, **kwargs):
        Curve.__init__(self, start, end, **kwargs)
        self.control1 = Point(control1) if control1 is not None else None
        self.control2 = Point(control2) if control1 is not None else None

    def __repr__(self):
        return "CubicBezier(start=%s, control1=%s, control2=%s, end=%s)" % (
            repr(self.start),
            repr(self.control1),
            repr(self.control2),
            repr(self.end),
        )

    def __copy__(self):
        return CubicBezier(
            self.start,
            self.control1,
            self.control2,
            self.end,
            relative=self.relative,
            smooth=self.smooth,
        )

    def __eq__(self, other):
        if not isinstance(other, CubicBezier):
            return NotImplemented
        return (
            self.start == other.start
            and self.end == other.end
            and self.control1 == other.control1
            and self.control2 == other.control2
        )

    def __ne__(self, other):
        if not isinstance(other, CubicBezier):
            return NotImplemented
        return not self == other

    def __imul__(self, other):
        if isinstance(other, str):
            other = Matrix(other)
        if isinstance(other, Matrix):
            if self.start is not None:
                self.start *= other
            if self.control1 is not None:
                self.control1 *= other
            if self.control2 is not None:
                self.control2 *= other
            if self.end is not None:
                self.end *= other
        return self

    def __len__(self):
        return 4

    def __getitem__(self, item):
        if item == 0:
            return self.start
        elif item == 1:
            return self.control1
        elif item == 2:
            return self.control2
        elif item == 3:
            return self.end
        else:
            raise IndexError

    def reverse(self):
        PathSegment.reverse(self)
        c2 = self.control2
        self.control2 = self.control1
        self.control1 = c2

    def npoint(self, positions):
        """Calculate the x,y position at a certain position of the path. `pos` may be a
        float or a NumPy array."""
        x0, y0 = self.start
        x1, y1 = self.control1
        x2, y2 = self.control2
        x3, y3 = self.end

        def _compute_point(position):
            # compute factors
            pos_3 = position ** 3
            n_pos = 1 - position
            n_pos_3 = n_pos ** 3
            pos_2_n_pos = position * position * n_pos
            n_pos_2_pos = n_pos * n_pos * position
            return (
                n_pos_3 * x0 + 3 * (n_pos_2_pos * x1 + pos_2_n_pos * x2) + pos_3 * x3,
                n_pos_3 * y0 + 3 * (n_pos_2_pos * y1 + pos_2_n_pos * y2) + pos_3 * y3,
            )

        try:
            import numpy as np

            xy = np.empty(shape=(len(positions), 2))
            xy[:, 0], xy[:, 1] = _compute_point(np.array(positions))
            return xy
        except ImportError:
            return [Point(*_compute_point(position)) for position in positions]

    def bbox(self):
        """returns the tight fitting bounding box of the bezier curve.
        Code by:
        https://github.com/mathandy/svgpathtools
        """
        xmin, xmax = self._real_minmax(0)
        ymin, ymax = self._real_minmax(1)
        return xmin, ymin, xmax, ymax

    def _real_minmax(self, v):
        """returns the minimum and maximum for a real cubic bezier, with a non-zero denom
        Code by:
        https://github.com/mathandy/svgpathtools
        """
        local_extremizers = [0, 1]
        a = [c[v] for c in self]
        denom = a[0] - 3 * a[1] + 3 * a[2] - a[3]
        if abs(denom) >= 1e-12:
            delta = a[1] ** 2 - (a[0] + a[1]) * a[2] + a[2] ** 2 + (a[0] - a[1]) * a[3]
            if delta >= 0:  # otherwise no local extrema
                sqdelta = sqrt(delta)
                tau = a[0] - 2 * a[1] + a[2]
                r1 = (tau + sqdelta) / denom
                r2 = (tau - sqdelta) / denom
                if 0 < r1 < 1:
                    local_extremizers.append(r1)
                if 0 < r2 < 1:
                    local_extremizers.append(r2)
        else:
            local_extremizers.append(0.5)
        local_extrema = [self.point(t)[v] for t in local_extremizers]
        return min(local_extrema), max(local_extrema)

    def _length_scipy(self, error=ERROR):
        from scipy.integrate import quad

        p0 = complex(*self.start)
        p1 = complex(*self.control1)
        p2 = complex(*self.control2)
        p3 = complex(*self.end)

        def _abs_derivative(t):
            return abs(
                3 * (p1 - p0) * (1 - t) ** 2
                + 6 * (p2 - p1) * (1 - t) * t
                + 3 * (p3 - p2) * t ** 2
            )

        return quad(_abs_derivative, 0.0, 1.0, epsabs=error, limit=1000)[0]

    def _length_default(self, error=ERROR, min_depth=MIN_DEPTH):
        return self._line_length(0, 1, error, min_depth)

    def length(self, error=ERROR, min_depth=MIN_DEPTH):
        """Calculate the length of the path up to a certain position"""
        try:
            return self._length_scipy(error)
        except:  # Fallback on any failure
            return self._length_default(error, min_depth)

    def is_smooth_from(self, previous):
        """Checks if this segment would be a smooth segment following the previous"""
        if isinstance(previous, CubicBezier):
            return self.start == previous.end and (self.control1 - self.start) == (
                previous.end - previous.control2
            )
        else:
            return self.control1 == self.start

    def d(self, current_point=None, relative=None, smooth=None):
        if (smooth is None and self.smooth) or (smooth is not None and smooth):
            if (
                current_point is None
                or (relative is None and self.relative)
                or (relative is not None and not relative)
            ):
                return "S %s %s" % (self.control2, self.end)
            else:
                return "s %s %s" % (
                    self.control2 - current_point,
                    self.end - current_point,
                )
        else:
            if (
                current_point is None
                or (relative is None and self.relative)
                or (relative is not None and not relative)
            ):
                return "C %s %s %s" % (self.control1, self.control2, self.end)
            else:
                return "c %s %s %s" % (
                    self.control1 - current_point,
                    self.control2 - current_point,
                    self.end - current_point,
                )


class Arc(Curve):
    def __init__(self, *args, **kwargs):
        """
        Represents Arc commands.

        Arc objects can take different parameters to create arcs.
        Since we expect taking in SVG parameters. We accept SVG parameterization which is:
        start, rx, ry, rotation, arc_flag, sweep_flag, end.

        To do matrix transitions, the native parameterization is start, end, center, prx, pry, sweep

        'start, end, center, prx, pry' are points and sweep amount is a t value in tau radians.
        If points are modified by an affine transformation, the arc is transformed.
        There is a special case for when the scale factor inverts, it inverts the sweep.

        Note: t-values are not angles from center in elliptical arcs. These are the same thing in
        circular arcs. But, here t is a parameterization around the ellipse, as if it were a circle.
        The position on the arc is (a * cos(t), b * sin(t)). If r-major was 0 for example. The
        positions would all fall on the x-axis. And the angle from center would all be either 0 or
        tau/2. However, since t is the parameterization we can conceptualize it as a position on a
        circle which is then scaled and rotated by a matrix.

        prx is the point at t 0 in the ellipse.
        pry is the point at t tau/4 in the ellipse.
        prx -> center -> pry should form a right triangle.

        The rotation can be defined as the angle from center to prx. Since prx is located at
        t(0) its deviation can only be the result of a rotation.

        Sweep is a value in t.
        The sweep angle can be a value greater than tau and less than -tau.
        However if this is the case, conversion back to Path.d() is expected to fail.
        We can denote these arc events but not as a single command.

        start_t + sweep = end_t
        """
        Curve.__init__(self, **kwargs)
        self.center = None
        self.prx = None
        self.pry = None
        self.sweep = None
        if len(args) == 6 and isinstance(args[1], complex):
            self._svg_complex_parameterize(*args)
            return
        elif len(kwargs) == 6 and "rotation" in kwargs:
            self._svg_complex_parameterize(**kwargs)
            return
        elif len(args) == 7:
            # This is an svg parameterized call.
            # A: rx ry x-axis-rotation large-arc-flag sweep-flag x y
            self._svg_parameterize(
                args[0], args[1], args[2], args[3], args[4], args[5], args[6]
            )
            return
        if (
            "left" in kwargs
            and "right" in kwargs
            and "top" in kwargs
            and "bottom" in kwargs
        ):
            left = kwargs["left"]
            right = kwargs["right"]
            top = kwargs["top"]
            bottom = kwargs["bottom"]
            self.center = Point((left + right) / 2.0, (top + bottom) / 2.0)
            rx = (right - left) / 2.0
            ry = (bottom - top) / 2.0
            self.prx = Point(self.center.x + rx, self.center.y)
            self.pry = Point(self.center.x, self.center.y + ry)
        len_args = len(args)
        if len_args > 0:
            if args[0] is not None:
                self.start = Point(args[0])
        if len_args > 1:
            if args[1] is not None:
                self.end = Point(args[1])
        if len_args > 2:
            if args[2] is not None:
                self.center = Point(args[2])
        if len_args > 3:
            if args[3] is not None:
                self.prx = Point(args[3])
        if len_args > 4:
            if args[4] is not None:
                self.pry = Point(args[4])
        if len_args > 5:
            self.sweep = args[5]
            return  # The args gave us everything.
        if "start" in kwargs:
            self.start = Point(kwargs["start"])
        if "end" in kwargs:
            self.end = Point(kwargs["end"])
        if "center" in kwargs:
            self.center = Point(kwargs["center"])
        if "prx" in kwargs:
            self.prx = Point(kwargs["prx"])
        if "pry" in kwargs:
            self.pry = Point(kwargs["pry"])
        if "sweep" in kwargs:
            self.sweep = kwargs["sweep"]
        cw = True  # Clockwise default. (sometimes needed)
        if self.start is not None and self.end is not None and self.center is None:
            # Start and end, but no center.
            # Solutions require a radius, a control point, or a bulge
            control = None
            sagitta = None
            if "bulge" in kwargs:
                bulge = float(kwargs["bulge"])
                sagitta = bulge * self.start.distance_to(self.end) / 2.0
            elif "sagitta" in kwargs:
                sagitta = float(kwargs["sagitta"])
            if sagitta is not None:
                control = Point.towards(self.start, self.end, 0.5)
                angle = self.start.angle_to(self.end)
                control = control.polar_to(angle - tau / 4.0, sagitta)
            if "control" in kwargs:  # Control is any additional point on the arc.
                control = Point(kwargs["control"])
            if control is not None:
                delta_a = control - self.start
                delta_b = self.end - control
                try:
                    slope_a = delta_a.y / delta_a.x
                except ZeroDivisionError:
                    slope_a = float("inf")
                try:
                    slope_b = delta_b.y / delta_b.x
                except ZeroDivisionError:
                    slope_b = float("inf")
                ab_mid = Point.towards(self.start, control, 0.5)
                bc_mid = Point.towards(control, self.end, 0.5)
                if delta_a.y == 0:  # slope_a == 0
                    cx = ab_mid.x
                    if delta_b.x == 0:  # slope_b == inf
                        cy = bc_mid.y
                    else:
                        cy = bc_mid.y + (bc_mid.x - cx) / slope_b
                elif delta_b.y == 0:  # slope_b == 0
                    cx = bc_mid.x
                    if delta_a.y == 0:  # slope_a == inf
                        cy = ab_mid.y
                    else:
                        cy = ab_mid.y + (ab_mid.x - cx) / slope_a
                elif delta_a.x == 0:  # slope_a == inf
                    cy = ab_mid.y
                    cx = slope_b * (bc_mid.y - cy) + bc_mid.x
                elif delta_b.x == 0:  # slope_b == inf
                    cy = bc_mid.y
                    cx = slope_a * (ab_mid.y - cy) + ab_mid.x
                elif slope_a == slope_b:
                    cx = ab_mid.x
                    cy = ab_mid.y
                else:
                    cx = (
                        slope_a * slope_b * (ab_mid.y - bc_mid.y)
                        - slope_a * bc_mid.x
                        + slope_b * ab_mid.x
                    ) / (slope_b - slope_a)
                    cy = ab_mid.y - (cx - ab_mid.x) / slope_a
                self.center = Point(cx, cy)
                cw = bool(Point.orientation(self.start, control, self.end) == 2)
            elif "r" in kwargs:
                r = kwargs["r"]
                mid = Point(
                    (self.start.x + self.end.x) / 2.0, (self.start.y + self.end.y) / 2.0
                )
                q = Point.distance(self.start, self.end)
                hq = q / 2.0
                if r < hq:
                    kwargs["r"] = r = hq  # Correct potential math domain error.
                self.center = Point(
                    mid.x + sqrt(r ** 2 - hq ** 2) * (self.start.y - self.end.y) / q,
                    mid.y + sqrt(r ** 2 - hq ** 2) * (self.end.x - self.start.x) / q,
                )
                cw = bool(Point.orientation(self.start, self.center, self.end) == 1)
                if "ccw" in kwargs and kwargs["ccw"] and cw or not cw:
                    # ccw arg exists, is true, and we found the cw center, or we didn't find the cw center.
                    self.center = Point(
                        mid.x
                        - sqrt(r ** 2 - hq ** 2) * (self.start.y - self.end.y) / q,
                        mid.y
                        - sqrt(r ** 2 - hq ** 2) * (self.end.x - self.start.x) / q,
                    )
            elif "rx" in kwargs and "ry" in kwargs:
                # This formulation will assume p1 and p2 are both axis aligned.
                # rx = kwargs["rx"]
                # ry = kwargs["ry"]
                # We will assume rx == abs(self.start.x - self.end.x)
                self.sweep = tau / 4.0
                self.center = Point(self.start.x, self.end.y)
                cw = bool(Point.orientation(self.start, self.center, self.end) == 1)
                if "scooped" in kwargs and kwargs["scooped"]:
                    self.sweep = -self.sweep
                    cw = not cw
                if ("ccw" in kwargs and kwargs["ccw"] and cw) or not cw:
                    self.center = Point(self.end.x, self.start.y)

        if self.center is None:
            raise ValueError("Not enough values to solve for center.")
        if "r" in kwargs:
            r = kwargs["r"]
            if self.prx is None:
                self.prx = Point(self.center.x + r, self.center.y)
            if self.pry is None:
                self.pry = Point(self.center.x, self.center.y + r)
        if "rx" in kwargs:
            rx = kwargs["rx"]
            if self.prx is None:
                if "rotation" in kwargs:
                    theta = kwargs["rotation"]
                    self.prx = Point.polar(self.center, theta, rx)
                else:
                    self.prx = Point(self.center.x + rx, self.center.y)
        if "ry" in kwargs:
            ry = kwargs["ry"]
            if self.pry is None:
                if "rotation" in kwargs:
                    theta = kwargs["rotation"]
                    theta += tau / 4.0
                    self.pry = Point.polar(self.center, theta, ry)
                else:
                    self.pry = Point(self.center.x, self.center.y + ry)
        if self.start is not None and (self.prx is None or self.pry is None):
            radius_s = Point.distance(self.center, self.start)
            self.prx = Point(self.center.x + radius_s, self.center.y)
            self.pry = Point(self.center.x, self.center.y + radius_s)
        if self.end is not None and (self.prx is None or self.pry is None):
            radius_e = Point.distance(self.center, self.end)
            self.prx = Point(self.center.x + radius_e, self.center.y)
            self.pry = Point(self.center.x, self.center.y + radius_e)
        if self.sweep is None and self.start is not None and self.end is not None:
            start_t = self.get_start_t()
            end_t = self.get_end_t()
            self.sweep = end_t - start_t
            if "ccw" in kwargs:
                cw = not bool(kwargs["ccw"])
            if cw and self.sweep < 0:
                self.sweep += tau
            if not cw and self.sweep > 0:
                self.sweep -= tau
        if self.sweep is not None and self.start is not None and self.end is None:
            start_t = self.get_start_t()
            end_t = start_t + self.sweep
            self.end = self.point_at_t(end_t)
        if self.sweep is not None and self.start is None and self.end is not None:
            end_t = self.get_end_t()
            start_t = end_t - self.sweep
            self.end = self.point_at_t(start_t)

    def __repr__(self):
        return "Arc(%s, %s, %s, %s, %s, %s)" % (
            repr(self.start),
            repr(self.end),
            repr(self.center),
            repr(self.prx),
            repr(self.pry),
            self.sweep,
        )

    def __copy__(self):
        return Arc(
            self.start,
            self.end,
            self.center,
            self.prx,
            self.pry,
            self.sweep,
            relative=self.relative,
        )

    def __eq__(self, other):
        if not isinstance(other, Arc):
            return NotImplemented
        return (
            self.start == other.start
            and self.end == other.end
            and self.prx == other.prx
            and self.pry == other.pry
            and self.center == other.center
            and self.sweep == other.sweep
        )

    def __ne__(self, other):
        if not isinstance(other, Arc):
            return NotImplemented
        return not self == other

    def __imul__(self, other):
        if isinstance(other, str):
            other = Matrix(other)
        if isinstance(other, Matrix):
            if self.start is not None:
                self.start *= other
            if self.center is not None:
                self.center *= other
            if self.end is not None:
                self.end *= other
            if self.prx is not None:
                self.prx *= other
            if self.pry is not None:
                self.pry *= other
            if other.value_scale_x() < 0:
                self.sweep = -self.sweep
            if other.value_scale_y() < 0:
                self.sweep = -self.sweep
        return self

    def __len__(self):
        return 5

    def __getitem__(self, item):
        if item == 0:
            return self.start
        elif item == 1:
            return self.end
        elif item == 2:
            return self.center
        elif item == 3:
            return self.prx
        elif item == 4:
            return self.pry
        raise IndexError

    @property
    def theta(self):
        """legacy property"""
        return Angle.radians(self.get_start_t()).as_positive_degrees

    @property
    def delta(self):
        """legacy property"""
        return Angle.radians(self.sweep).as_degrees

    def reverse(self):
        PathSegment.reverse(self)
        self.sweep = -self.sweep

    def npoint(self, positions):
        try:
            import numpy as np

            return self._points_numpy(np.array(positions))
        except ImportError:
            if self.start == self.end and self.sweep == 0:
                # This is equivalent of omitting the segment
                return [self.start] * len(positions)

            start_t = self.get_start_t()
            return [
                self.start
                if pos == 0
                else self.end
                if pos == 1
                else self.point_at_t(start_t + self.sweep * pos)
                for pos in positions
            ]

    def _points_numpy(self, positions):
        """Vectorized version of `point()`.

        :param positions: 1D numpy array of float in [0, 1]
        :return: 1D numpy array of complex
        """
        import numpy as np

        xy = np.empty((len(positions), 2), dtype=float)

        if self.start == self.end and self.sweep == 0:
            xy[:, 0], xy[:, 1] = self.start
        else:
            t = self.get_start_t() + self.sweep * positions

            rotation = self.get_rotation()
            a = self.rx
            b = self.ry
            cx = self.center.x
            cy = self.center.y
            cos_rot = cos(rotation)
            sin_rot = sin(rotation)
            cos_t = np.cos(t)
            sin_t = np.sin(t)
            xy[:, 0] = cx + a * cos_t * cos_rot - b * sin_t * sin_rot
            xy[:, 1] = cy + a * cos_t * sin_rot + b * sin_t * cos_rot

            # ensure clean endings
            xy[positions == 0, :] = list(self.start)
            xy[positions == 1, :] = list(self.end)

        return xy

    def _integral_length(self):
        def ellipse_part_integral(t1, t2, a, b, n=100000):
            # function to integrate
            def f(t):
                return sqrt(1 - (1 - a ** 2 / b ** 2) * sin(t) ** 2)

            start = min(t1, t2)
            seg_len = abs(t1 - t2) / n
            return b * sum(f(start + seg_len * i) * seg_len for i in range(1, n + 1))

        start_angle = self.get_start_t()
        end_angle = start_angle + self.sweep
        return ellipse_part_integral(start_angle, end_angle, self.rx, self.ry)

    def _exact_length(self):
        """scipy is not a dependency. However, if scipy exists this function will find the
        exact arc length. By default .length() delegates to here and on failure uses the
        fallback method."""
        from scipy.special import ellipeinc

        a = self.rx
        b = self.ry
        phi = self.get_start_t()
        m = 1 - (a / b) ** 2
        d1 = ellipeinc(phi, m)
        phi = phi + self.sweep
        m = 1 - (a / b) ** 2
        d2 = ellipeinc(phi, m)
        return b * abs(d2 - d1)

    def length(self, error=ERROR, min_depth=MIN_DEPTH):
        """The length of an elliptical arc segment requires numerical
        integration, and in that case it's simpler to just do a geometric
        approximation, as for cubic bezier curves.
        """
        if self.sweep == 0:
            return 0
        if self.start == self.end and self.sweep == 0:
            # This is equivalent of omitting the segment
            return 0
        a = self.rx
        b = self.ry
        d = abs(a - b)

        if d < ERROR:  # This is a circle.
            return abs(self.rx * self.sweep)
        try:
            return self._exact_length()
        except:  # Fallback on any failure
            return self._line_length(error=error, min_depth=min_depth)

    def _svg_complex_parameterize(
        self, start, radius, rotation, arc_flag, sweep_flag, end
    ):
        """Parameterization with complex radius and having rotation factors."""
        self._svg_parameterize(
            Point(start),
            radius.real,
            radius.imag,
            rotation,
            arc_flag,
            sweep_flag,
            Point(end),
        )

    def _svg_parameterize(
        self, start, rx, ry, rotation, large_arc_flag, sweep_flag, end
    ):
        """Conversion from svg parameterization, our chosen native native form.
        http://www.w3.org/TR/SVG/implnote.html#ArcImplementationNotes"""

        large_arc_flag = bool(large_arc_flag)
        sweep_flag = bool(sweep_flag)
        start = Point(start)
        self.start = start
        end = Point(end)
        self.end = end
        if start == end or rx == 0 or ry == 0:
            # If start is equal to end, there are infinite number of circles so these void out.
            # We still permit this kind of arc, but SVG parameterization cannot be used to achieve it.
            self.sweep = 0
            self.prx = Point(start)
            self.pry = Point(start)
            self.center = Point(start)
            return
        cosr = cos(radians(rotation))
        sinr = sin(radians(rotation))
        dx = (start.real - end.real) / 2
        dy = (start.imag - end.imag) / 2
        x1prim = cosr * dx + sinr * dy
        x1prim_sq = x1prim * x1prim
        y1prim = -sinr * dx + cosr * dy
        y1prim_sq = y1prim * y1prim

        rx_sq = rx * rx
        ry_sq = ry * ry

        # Correct out of range radii
        radius_check = (x1prim_sq / rx_sq) + (y1prim_sq / ry_sq)
        if radius_check > 1:
            rx *= sqrt(radius_check)
            ry *= sqrt(radius_check)
            rx_sq = rx * rx
            ry_sq = ry * ry

        t1 = rx_sq * y1prim_sq
        t2 = ry_sq * x1prim_sq
        c = sqrt(abs((rx_sq * ry_sq - t1 - t2) / (t1 + t2)))

        if large_arc_flag == sweep_flag:
            c = -c
        cxprim = c * rx * y1prim / ry
        cyprim = -c * ry * x1prim / rx

        center = Point(
            (cosr * cxprim - sinr * cyprim) + ((start.real + end.real) / 2),
            (sinr * cxprim + cosr * cyprim) + ((start.imag + end.imag) / 2),
        )

        ux = (x1prim - cxprim) / rx
        uy = (y1prim - cyprim) / ry
        vx = (-x1prim - cxprim) / rx
        vy = (-y1prim - cyprim) / ry
        n = sqrt(ux * ux + uy * uy)
        p = ux
        # theta = degrees(acos(p / n))
        # if uy < 0:
        #     theta = -theta
        # theta = theta % 360

        n = sqrt((ux * ux + uy * uy) * (vx * vx + vy * vy))
        p = ux * vx + uy * vy
        d = p / n
        # In certain cases the above calculation can through inaccuracies
        # become just slightly out of range, f ex -1.0000000000000002.
        if d > 1.0:
            d = 1.0
        elif d < -1.0:
            d = -1.0
        delta = degrees(acos(d))
        if (ux * vy - uy * vx) < 0:
            delta = -delta
        delta = delta % 360
        if not sweep_flag:
            delta -= 360
        # built parameters, delta, theta, center

        rotate_matrix = Matrix()
        rotate_matrix.post_rotate(
            Angle.degrees(rotation).as_radians, center.x, center.y
        )

        self.center = center
        self.prx = Point(center.x + rx, center.y)
        self.pry = Point(center.x, center.y + ry)

        self.prx.matrix_transform(rotate_matrix)
        self.pry.matrix_transform(rotate_matrix)
        self.sweep = Angle.degrees(delta).as_radians

    def as_quad_curves(self, arc_required=None):
        if arc_required is None:
            sweep_limit = tau / 12.0
            arc_required = int(ceil(abs(self.sweep) / sweep_limit))
        if arc_required == 0:
            return
        t_slice = self.sweep / float(arc_required)

        current_t = self.get_start_t()
        p_start = self.start

        theta = self.get_rotation()
        cos_theta = cos(theta)
        sin_theta = sin(theta)

        a = self.rx
        b = self.ry
        cx = self.center.x
        cy = self.center.y

        for i in range(0, arc_required):
            next_t = current_t + t_slice
            mid_t = (next_t + current_t) / 2
            p_end = self.point_at_t(next_t)
            if i == arc_required - 1:
                p_end = self.end
            cos_mid_t = cos(mid_t)
            sin_mid_t = sin(mid_t)
            alpha = (4.0 - cos(t_slice)) / 3.0
            px = cx + alpha * (a * cos_mid_t * cos_theta - b * sin_mid_t * sin_theta)
            py = cy + alpha * (a * cos_mid_t * sin_theta + b * sin_mid_t * cos_theta)
            yield QuadraticBezier(p_start, (px, py), p_end)
            p_start = p_end
            current_t = next_t

    def as_cubic_curves(self, arc_required=None):
        if arc_required is None:
            sweep_limit = tau / 12.0
            arc_required = int(ceil(abs(self.sweep) / sweep_limit))
        if arc_required == 0:
            return
        t_slice = self.sweep / float(arc_required)

        theta = self.get_rotation()
        rx = self.rx
        ry = self.ry
        p_start = self.start
        current_t = self.get_start_t()
        x0 = self.center.x
        y0 = self.center.y
        cos_theta = cos(theta)
        sin_theta = sin(theta)

        for i in range(0, arc_required):
            next_t = current_t + t_slice

            alpha = sin(t_slice) * (sqrt(4 + 3 * pow(tan(t_slice / 2.0), 2)) - 1) / 3.0

            cos_start_t = cos(current_t)
            sin_start_t = sin(current_t)

            ePrimen1x = -rx * cos_theta * sin_start_t - ry * sin_theta * cos_start_t
            ePrimen1y = -rx * sin_theta * sin_start_t + ry * cos_theta * cos_start_t

            cos_end_t = cos(next_t)
            sin_end_t = sin(next_t)

            p2En2x = x0 + rx * cos_end_t * cos_theta - ry * sin_end_t * sin_theta
            p2En2y = y0 + rx * cos_end_t * sin_theta + ry * sin_end_t * cos_theta
            p_end = (p2En2x, p2En2y)
            if i == arc_required - 1:
                p_end = self.end

            ePrimen2x = -rx * cos_theta * sin_end_t - ry * sin_theta * cos_end_t
            ePrimen2y = -rx * sin_theta * sin_end_t + ry * cos_theta * cos_end_t

            p_c1 = (p_start[0] + alpha * ePrimen1x, p_start[1] + alpha * ePrimen1y)
            p_c2 = (p_end[0] - alpha * ePrimen2x, p_end[1] - alpha * ePrimen2y)

            yield CubicBezier(p_start, p_c1, p_c2, p_end)
            p_start = Point(p_end)
            current_t = next_t

    def is_circular(self):
        a = self.rx
        b = self.ry
        return a == b

    @property
    def radius(self):
        """Legacy complex radius property

        Point will work like a complex for legacy reasons.
        """
        return Point(self.rx, self.ry)

    @property
    def rx(self):
        return Point.distance(self.center, self.prx)

    @property
    def ry(self):
        return Point.distance(self.center, self.pry)

    def get_rotation(self):
        return Point.angle(self.center, self.prx)

    def get_start_angle(self):
        """
        :return: Angle from the center point to start point.
        """
        return self.angle_at_point(self.start)

    def get_end_angle(self):
        """
        :return: Angle from the center point to end point.
        """
        return self.angle_at_point(self.end)

    def get_start_t(self):
        """
        start t value in the ellipse.

        :return: t parameter of start point.
        """
        return self.t_at_point(self.point_at_angle(self.get_start_angle()))

    def get_end_t(self):
        """
        end t value in the ellipse.

        :return: t parameter of start point.
        """
        return self.t_at_point(self.point_at_angle(self.get_end_angle()))

    def point_at_angle(self, angle):
        """
        find the point on the ellipse from the center at the given angle.
        Note: For non-circular arcs this is different than point(t).

        :param angle: angle from center to find point
        :return: point found
        """
        angle -= self.get_rotation()
        a = self.rx
        b = self.ry
        if a == b:
            return self.point_at_t(angle)
        t = atan2(a * tan(angle), b)
        tau_1_4 = tau / 4.0
        tau_3_4 = 3 * tau_1_4
        if tau_3_4 >= abs(angle) % tau > tau_1_4:
            t += tau / 2.0
        return self.point_at_t(t)

    def angle_at_point(self, p):
        """
        find the angle to the point.

        :param p: point
        :return: angle to given point.
        """
        return self.center.angle_to(p)

    def t_at_point(self, p):
        """
        find the t parameter to at the point.

        :param p: point
        :return: t parameter to the given point.
        """
        angle = self.angle_at_point(p)
        angle -= self.get_rotation()
        a = self.rx
        b = self.ry
        t = atan2(a * tan(angle), b)
        tau_1_4 = tau / 4.0
        tau_3_4 = 3 * tau_1_4
        if tau_3_4 >= abs(angle) % tau > tau_1_4:
            t += tau / 2.0
        return t

    def point_at_t(self, t):
        """
        find the point that corresponds to given value t.
        Where t=0 is the first point and t=tau is the final point.

        In the case of a circle: t = angle.

        :param t:
        :return:
        """
        rotation = self.get_rotation()
        a = self.rx
        b = self.ry
        cx = self.center.x
        cy = self.center.y
        cos_rot = cos(rotation)
        sin_rot = sin(rotation)
        cos_t = cos(t)
        sin_t = sin(t)
        px = cx + a * cos_t * cos_rot - b * sin_t * sin_rot
        py = cy + a * cos_t * sin_rot + b * sin_t * cos_rot
        return Point(px, py)

    def get_ellipse(self):
        return Ellipse(self.center, self.rx, self.ry, self.get_rotation())

    def bbox(self):
        """Find the bounding box of a arc.
        Code from: https://github.com/mathandy/svgpathtools
        """
        if self.sweep == 0:
            return self.start.x, self.start.y, self.end.x, self.end.y
        phi = self.get_rotation().as_radians
        if cos(phi) == 0:
            atan_x = tau / 4.0
            atan_y = 0
        elif sin(phi) == 0:
            atan_x = 0
            atan_y = tau / 4.0
        else:
            rx, ry = self.rx, self.ry
            atan_x = atan(-(ry / rx) * tan(phi))
            atan_y = atan((ry / rx) / tan(phi))

        def angle_inv(ang, k):  # inverse of angle from Arc.derivative()
            return ((ang + (tau / 2.0) * k) * (360 / tau) - self.theta) / self.delta

        xtrema = [self.start.x, self.end.x]
        ytrema = [self.start.y, self.end.y]

        for k in range(-4, 5):
            tx = angle_inv(atan_x, k)
            ty = angle_inv(atan_y, k)
            if 0 <= tx <= 1:
                xtrema.append(self.point(tx).x)
            if 0 <= ty <= 1:
                ytrema.append(self.point(ty).y)

        return min(xtrema), min(ytrema), max(xtrema), max(ytrema)

    def d(self, current_point=None, relative=None, smooth=None):
        if (
            current_point is None
            or (relative is None and self.relative)
            or (relative is not None and not relative)
        ):
            return "A %G,%G %G %d,%d %s" % (
                self.rx,
                self.ry,
                self.get_rotation().as_degrees,
                int(abs(self.sweep) > (tau / 2.0)),
                int(self.sweep >= 0),
                self.end,
            )
        else:
            return "a %G,%G %G %d,%d %s" % (
                self.rx,
                self.ry,
                self.get_rotation().as_degrees,
                int(abs(self.sweep) > (tau / 2.0)),
                int(self.sweep >= 0),
                self.end - current_point,
            )


class Path(MutableSequence):
    """
    A Path is a Mutable sequence of path segments

    It is a generalized shape which can map out all the other shapes.

    Each PathSegment object maps a particular command. Each one exists only once in each path and every point contained
    within the object is also unique. We attempt to internally maintain some validity. Each end point should link
    to the following segments start point. And each close point should connect from the preceding segments endpoint to
    the last Move command.

    These are soft checks made only at the time of addition and some manipulations. Modifying the points of the segments
    can and will cause path invalidity. Some SVG invalid operations are permitted such as arcs longer than tau radians
    or beginning sequences without a move. The expectation is that these will eventually be used as part of a valid path
    so these fragment paths are permitted. In some cases these invalid paths will still have consistent path_d values,
    in other cases, there will be no valid methods to reproduce these.
    """

    def __init__(self, *args, **kwargs):
        self._length = None
        self._lengths = None
        self._segments = list()
        if len(args) != 1:
            self._segments.extend(args)
        else:
            s = args[0]
            if isinstance(s, Subpath):
                self._segments.extend(s.segments(transformed=False))
            elif isinstance(s, str):
                self._segments = list()
                self.parse(s)
            elif isinstance(s, tuple):
                # We have no guarantee of the validity of the source data
                self._segments.extend(s)
                self.validate_connections()
            elif isinstance(s, list):
                # We have no guarantee of the validity of the source data
                self._segments.extend(s)
                self.validate_connections()
            elif isinstance(s, PathSegment):
                self._segments.append(s)

    def __copy__(self):
        path = Path(self)
        segs = path._segments
        for i in range(0, len(segs)):
            segs[i] = copy(segs[i])
        return path

    def __getitem__(self, index):
        return self._segments[index]

    #
    # def __eq__(self, other):
    #     if not isinstance(other, Path):
    #         return NotImplemented
    #     if self.fill != other.fill or self.stroke != other.stroke:
    #         return False
    #     first = self
    #     if not isinstance(first, Path):
    #         first = Path(first)
    #     second = other
    #     if not isinstance(second, Path):
    #         second = Path(second)
    #     return first == second
    #
    # def __ne__(self, other):
    #     if not isinstance(other, Path):
    #         return NotImplemented
    #     return not self == other

    # def __iadd__(self, other):
    #     if isinstance(other, Path):
    #         return Path(self) + Path(other)
    #     return NotImplemented
    #
    # __add__ = __iadd__

    def __matmul__(self, other):
        m = copy(self)
        m.__imatmul__(other)
        return m

    def __rmatmul__(self, other):
        m = copy(other)
        m.__imatmul__(self)
        return m

    def __imatmul__(self, other):
        """
        The % operation with a matrix works much like multiplication except that it automatically reifies the shape.
        """
        if isinstance(other, str):
            other = Matrix(other)
        if isinstance(other, Matrix):
            self.transform *= other
        return self

    def __setitem__(self, index, new_element):
        if isinstance(new_element, str):
            new_element = Path(new_element)
            if len(new_element) == 0:
                return
            new_element = new_element.segments()
            if isinstance(index, int):
                if len(new_element) > 1:
                    raise ValueError  # Cannot insert multiple items into a single space. Requires slice.
                new_element = new_element[0]
        self._segments[index] = new_element
        self._length = None
        self._lengths = None
        if isinstance(index, slice):
            self.validate_connections()
        else:
            self._validate_connection(index - 1)
            self._validate_connection(index)
            if isinstance(new_element, Move):
                self._validate_move(index)
            if isinstance(new_element, Close):
                self._validate_close(index)

    def __delitem__(self, index):
        original_element = self._segments[index]
        del self._segments[index]
        self._length = None
        if isinstance(index, slice):
            self.validate_connections()
        else:
            self._validate_connection(index - 1)
            if isinstance(original_element, (Close, Move)):
                self._validate_subpath(index)

    def __iadd__(self, other):
        if isinstance(other, str):
            self.parse(other)
        elif isinstance(other, (Path, Subpath)):
            self.extend(map(copy, list(other)))
        elif isinstance(other, PathSegment):
            self.append(other)
        else:
            return NotImplemented
        return self

    def __add__(self, other):
        if isinstance(other, (str, Path, Subpath, PathSegment)):
            n = copy(self)
            n += other
            return n
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, str):
            path = Path(other)
            path.extend(map(copy, self._segments))
            return path
        elif isinstance(other, PathSegment):
            path = copy(self)
            path.insert(0, other)
            return path
        else:
            return NotImplemented

    def __len__(self):
        return len(self._segments)

    def __str__(self):
        return self.d()

    def __repr__(self):
        return "%s(%s)" % (self.__class__.__name__, self.d())

    def __eq__(self, other):
        if isinstance(other, str):
            return self.__eq__(Path(other))
        if not isinstance(other, Path):
            return NotImplemented
        if len(self) != len(other):
            return False
        p = self
        q = other
        for s, o in zip(q._segments, p._segments):
            if not s == o:
                return False
        return True

    def __ne__(self, other):
        if not isinstance(other, (Path, str)):
            return NotImplemented
        return not self == other

    def _calc_lengths(self, error=ERROR, min_depth=MIN_DEPTH, segments=None):
        """
        Calculate the length values for the segments of the Shape.

        :param error: error permitted for length calculations.
        :param min_depth: minimum depth for the length calculation.
        :param segments: optional segments to use.
        :return:
        """
        if segments is None:
            segments = self.segments(False)
        if self._length is not None:
            return
        lengths = [each.length(error=error, min_depth=min_depth) for each in segments]
        self._length = sum(lengths)
        if self._length == 0:
            self._lengths = lengths
        else:
            self._lengths = [each / self._length for each in lengths]

    def npoint(self, positions, error=ERROR):
        """
        Find a points between 0 and 1 within the shape. Numpy acceleration allows points to be an array of floats.
        """
        try:
            import numpy as np
        except ImportError:
            return [self.point(pos) for pos in positions]

        segments = self.segments(False)
        if len(segments) == 0:
            return None
        # Shortcuts
        if self._length is None:
            self._calc_lengths(error=error, segments=segments)
        xy = np.empty((len(positions), 2), dtype=float)
        if self._length == 0:
            i = int(round(positions * (len(segments) - 1)))
            point = segments[i].point(0.0)
            xy[:] = point
            return xy

        # Find which segment the point we search for is located on:
        segment_start = 0
        for index, segment in enumerate(segments):
            segment_end = segment_start + self._lengths[index]
            position_subset = (segment_start <= positions) & (positions < segment_end)
            v0 = positions[position_subset]
            if not len(v0):
                continue  # Nothing matched.
            d = segment_end - segment_start
            if d == 0:  # This segment is 0 length.
                segment_pos = 0.0
            else:
                segment_pos = (v0 - segment_start) / d
            c = segment.npoint(segment_pos)
            xy[position_subset] = c[:]
            segment_start = segment_end

        # the loop above will miss position == 1
        xy[positions == 1] = np.array(list(segments[-1].end))
        return xy

    def point(self, position, error=ERROR):
        """
        Find a point between 0 and 1 within the Shape, going through the shape with regard to position.

        :param position: value between 0 and 1 within the shape.
        :param error: Length error permitted.
        :return: Point at the given location.
        """
        segments = self.segments(False)
        if len(segments) == 0:
            return None
        # Shortcuts
        try:
            if position <= 0.0:
                return segments[0].point(position)
            if position >= 1.0:
                return segments[-1].point(position)
        except ValueError:
            return self.npoint([position], error=error)[0]

        if self._length is None:
            self._calc_lengths(error=error, segments=segments)

        if self._length == 0:
            i = int(round(position * (len(segments) - 1)))
            return segments[i].point(0.0)
        # Find which segment the point we search for is located on:
        segment_start = 0
        segment_pos = 0
        segment = segments[0]
        for index, segment in enumerate(segments):
            segment_end = segment_start + self._lengths[index]
            if segment_end >= position:
                # This is the segment! How far in on the segment is the point?
                segment_pos = (position - segment_start) / (segment_end - segment_start)
                break
            segment_start = segment_end
        return segment.point(segment_pos)

    def length(self, error=ERROR, min_depth=MIN_DEPTH):
        self._calc_lengths(error, min_depth)
        return self._length

    def _validate_subpath(self, index):
        """ensure the subpath containing this index is valid."""
        if index < 0 or index + 1 >= len(self._segments):
            return  # This connection doesn't exist.
        for j in range(index, len(self._segments)):
            close_search = self._segments[j]
            if isinstance(close_search, Move):
                return  # Not a closed path, subpath is valid.
            if isinstance(close_search, Close):
                for k in range(index, -1, -1):
                    move_search = self._segments[k]
                    if isinstance(move_search, Move):
                        self._segments[j].end = Point(move_search.end)
                        return
                self._segments[j].end = Point(self._segments[0].end)
                return

    def _validate_move(self, index):
        """ensure the next closed point from this index points to a valid location."""
        for i in range(index + 1, len(self._segments)):
            segment = self._segments[i]
            if isinstance(segment, Move):
                return  # Not a closed path, the move is valid.
            if isinstance(segment, Close):
                segment.end = Point(self._segments[index].end)
                return

    def _validate_close(self, index):
        """ensure the close element at this position correctly links to the previous move"""
        for i in range(index, -1, -1):
            segment = self._segments[i]
            if isinstance(segment, Move):
                self._segments[index].end = Point(segment.end)
                return
        self._segments[index].end = (
            Point(self._segments[0].end) if self._segments[0].end is not None else None
        )
        # If move is never found, just the end point of the first element. Unless that's not a thing.

    def _validate_connection(self, index, prefer_second=False):
        """
        Validates the connection at the index.
        Connection 0 is the connection between getitem(0) and getitem(1)

        prefer_second is for those cases where failing the connection requires replacing
        a existing value. It will prefer the authority of right side, second value.
        """
        if index < 0 or index + 1 >= len(self._segments):
            return  # This connection doesn't exist.
        first = self._segments[index]
        second = self._segments[index + 1]
        if first.end is not None and second.start is None:
            second.start = Point(first.end)
        elif first.end is None and second.start is not None:
            first.end = Point(second.start)
        elif first.end != second.start:
            # The two values exist but are not equal. One must replace the other.
            if prefer_second:
                first.end = Point(second.start)
            else:
                second.start = Point(first.end)

    def parse(self, pathdef):
        """Parses the SVG path."""
        tokens = PathLexicalParser()
        tokens.parse(self, pathdef)

    def validate_connections(self):
        """
        Force validate all connections.

        This will scan path connections and link any adjacent elements together by replacing any None points or causing
        the start position of the next element to equal the end position of the previous. This should only be needed
        when combining paths and elements together. Close elements are always connected to the last Move element or to
        the end position of the first element in the list. The start element of the first segment may or may not be
        None.
        """
        zpoint = None
        last_segment = None
        for segment in self._segments:
            if zpoint is None or isinstance(segment, Move):
                zpoint = segment.end
            if last_segment is not None:
                if segment.start is None and last_segment.end is not None:
                    segment.start = Point(last_segment.end)
                elif last_segment.end is None and segment.start is not None:
                    last_segment.end = Point(segment.start)
                elif last_segment.end != segment.start:
                    segment.start = Point(last_segment.end)
            if (
                isinstance(segment, Close)
                and zpoint is not None
                and segment.end != zpoint
            ):
                segment.end = Point(zpoint)
            last_segment = segment

    def _is_valid(self):
        """
        Checks validation of all connections.

        Paths are valid if all end points match the start of the next point and all close
        commands return to the last valid move command.

        This does not check for incongruent path validity. Path fragments without initial moves
        double closed paths, may all pass this check.
        """
        zpoint = None
        last_segment = None
        for segment in self._segments:
            if zpoint is None or isinstance(segment, Move):
                zpoint = segment.end
            if last_segment is not None:
                if segment.start is None:
                    return False
                elif last_segment.end is None:
                    return False
                elif last_segment.end != segment.start:
                    return False
            if (
                isinstance(segment, Close)
                and zpoint is not None
                and segment.end != zpoint
            ):
                return False
            last_segment = segment
        return True

    @property
    def first_point(self):
        """First point along the Path. This is the start point of the first segment unless it starts
        with a Move command with a None start in which case first point is that Move's destination."""
        if len(self._segments) == 0:
            return None
        if self._segments[0].start is not None:
            return Point(self._segments[0].start)
        return (
            Point(self._segments[0].end) if self._segments[0].end is not None else None
        )

    @property
    def current_point(self):
        if len(self._segments) == 0:
            return None
        return (
            Point(self._segments[-1].end)
            if self._segments[-1].end is not None
            else None
        )

    @property
    def z_point(self):
        """
        Z is the destination of the last Move. It can mean, but doesn't necessarily mean the first_point in the path.
        This behavior of Z is defined in svg spec:
        http://www.w3.org/TR/SVG/paths.html#PathDataClosePathCommand
        """
        end_pos = None
        for segment in reversed(self._segments):
            if isinstance(segment, Move):
                end_pos = segment.end
                break
        if end_pos is None:
            try:
                end_pos = self._segments[0].end
            except IndexError:
                pass  # entire path is "z".
        return end_pos

    @property
    def smooth_point(self):
        """Returns the smoothing control point for the smooth commands.
        With regards to the SVG standard if the last command was a curve the smooth
        control point is the reflection of the previous control point.

        If the last command was not a curve, the smooth_point is coincident with the current.
        https://www.w3.org/TR/SVG/paths.html#PathDataCubicBezierCommands
        """

        if len(self._segments) == 0:
            return None
        start_pos = self.current_point
        last_segment = self._segments[-1]
        if isinstance(last_segment, QuadraticBezier):
            previous_control = last_segment.control
            return previous_control.reflected_across(start_pos)
        elif isinstance(last_segment, CubicBezier):
            previous_control = last_segment.control2
            return previous_control.reflected_across(start_pos)
        return start_pos

    def start(self):
        pass

    def end(self):
        pass

    def move(self, *points, **kwargs):
        relative = kwargs["relative"] if "relative" in kwargs else False
        start_pos = self.current_point
        end_pos = points[0]
        if end_pos in ("z", "Z"):
            end_pos = self.z_point
        segment = Move(start_pos, end_pos)
        segment.relative = relative
        self.append(segment)
        if len(points) > 1:
            self.line(*points[1:], relative=relative)
        return self

    def line(self, *points, **kwargs):
        relative = kwargs["relative"] if "relative" in kwargs else False
        start_pos = self.current_point
        end_pos = points[0]
        if end_pos in ("z", "Z"):
            end_pos = self.z_point
        segment = Line(start_pos, end_pos)
        segment.relative = relative
        self.append(segment)
        if len(points) > 1:
            self.line(*points[1:])
        return self

    def vertical(self, *y_points, **kwargs):
        relative = kwargs["relative"] if "relative" in kwargs else False
        start_pos = self.current_point
        if relative:
            segment = Line(start_pos, Point(start_pos.x, start_pos.y + y_points[0]))
        else:
            segment = Line(start_pos, Point(start_pos.x, y_points[0]))
        segment.relative = relative
        self.append(segment)
        if len(y_points) > 1:
            self.vertical(*y_points[1:], relative=relative)
        return self

    def horizontal(self, *x_points, **kwargs):
        relative = kwargs["relative"] if "relative" in kwargs else False
        start_pos = self.current_point
        if relative:
            segment = Line(start_pos, Point(start_pos.x + x_points[0], start_pos.y))
            segment.relative = relative
        else:
            segment = Line(start_pos, Point(x_points[0], start_pos.y))
            segment.relative = relative
        self.append(segment)
        if len(x_points) > 1:
            self.horizontal(*x_points[1:], relative=relative)
        return self

    def smooth_quad(self, *points, **kwargs):
        """Smooth curve. First control point is the "reflection" of
        the second control point in the previous path."""
        relative = kwargs["relative"] if "relative" in kwargs else False
        start_pos = self.current_point
        control1 = self.smooth_point
        end_pos = points[0]
        if end_pos in ("z", "Z"):
            end_pos = self.z_point
        segment = QuadraticBezier(start_pos, control1, end_pos)
        segment.relative = relative
        segment.smooth = True
        self.append(segment)
        if len(points) > 1:
            self.smooth_quad(*points[1:])
        return self

    def quad(self, *points, **kwargs):
        relative = kwargs["relative"] if "relative" in kwargs else False
        start_pos = self.current_point
        control = points[0]
        if control in ("z", "Z"):
            control = self.z_point
        end_pos = points[1]
        if end_pos in ("z", "Z"):
            end_pos = self.z_point
        segment = QuadraticBezier(start_pos, control, end_pos)
        segment.relative = relative
        segment.smooth = False
        self.append(segment)
        if len(points) > 2:
            self.quad(*points[2:])
        return self

    def smooth_cubic(self, *points, **kwargs):
        """Smooth curve. First control point is the "reflection" of
        the second control point in the previous path."""
        relative = kwargs["relative"] if "relative" in kwargs else False
        start_pos = self.current_point
        control1 = self.smooth_point
        control2 = points[0]

        if control2 in ("z", "Z"):
            control2 = self.z_point
        end_pos = points[1]
        if end_pos in ("z", "Z"):
            end_pos = self.z_point
        segment = CubicBezier(start_pos, control1, control2, end_pos)
        segment.relative = relative
        segment.smooth = True
        self.append(segment)
        if len(points) > 2:
            self.smooth_cubic(*points[2:])
        return self

    def cubic(self, *points, **kwargs):
        relative = kwargs["relative"] if "relative" in kwargs else False
        start_pos = self.current_point
        control1 = points[0]
        if control1 in ("z", "Z"):
            control1 = self.z_point
        control2 = points[1]
        if control2 in ("z", "Z"):
            control2 = self.z_point
        end_pos = points[2]
        if end_pos in ("z", "Z"):
            end_pos = self.z_point
        segment = CubicBezier(start_pos, control1, control2, end_pos)
        segment.relative = relative
        segment.smooth = False
        self.append(segment)
        if len(points) > 3:
            self.cubic(*points[3:])
        return self

    def arc(self, *arc_args, **kwargs):
        relative = kwargs["relative"] if "relative" in kwargs else False
        start_pos = self.current_point
        rx = arc_args[0]
        ry = arc_args[1]
        rotation = arc_args[2]
        arc = arc_args[3]
        sweep = arc_args[4]
        end_pos = arc_args[5]
        if end_pos in ("z", "Z"):
            end_pos = self.z_point
        segment = Arc(start_pos, rx, ry, rotation, arc, sweep, end_pos)
        segment.relative = relative
        self.append(segment)
        if len(arc_args) > 6:
            self.arc(*arc_args[6:])
        return self

    def closed(self, relative=False):
        start_pos = self.current_point
        end_pos = self.z_point
        segment = Close(start_pos, end_pos)
        segment.relative = relative
        self.append(segment)
        return self

    def append(self, value):
        if isinstance(value, str):
            value = Path(value)
            if len(value) == 0:
                return
            if len(value) > 1:
                self.extend(value)
                return
            value = value[0]
        self._length = None
        index = len(self._segments) - 1
        self._segments.append(value)
        self._validate_connection(index)
        if isinstance(value, Close):
            self._validate_close(index + 1)

    def insert(self, index, value):
        if isinstance(value, str):
            value = Path(value)
            if len(value) == 0:
                return
            value = value[0]
        self._length = None
        self._segments.insert(index, value)
        self._validate_connection(index - 1)
        self._validate_connection(index)
        if isinstance(value, Move):
            self._validate_move(index)
        if isinstance(value, Close):
            self._validate_close(index)

    def extend(self, iterable):
        if isinstance(iterable, str):
            iterable = Path(iterable)
        self._length = None
        index = len(self._segments) - 1
        self._segments.extend(iterable)
        self._validate_connection(index)
        self._validate_subpath(index)

    def direct_close(self):
        """Forces close operations to be zero length by introducing a direct
        line to operation just before any non-zero length close.

        This is helpful because for some operations like reverse() because the
        close must located at the very end of the path sequence. But, if it's
        in effect a line-to and close, the line-to would need to start the sequence.

        But, for some operations this won't matter since it will still result in
        a closed shape with reversed ordering. But, if the final point in the
        sequence must exactly switch with the first point in the sequence. The
        close segments must be direct and zero length.
        """
        if len(self._segments) == 0:
            return
        for i in range(len(self._segments) - 1, -1, -1):
            segment = self._segments[i]
            if isinstance(segment, Close):
                if segment.length() != 0:
                    line = Line(segment.start, segment.end)
                    segment.start = Point(segment.end)
                    self.insert(i, line)
        return self

    def reverse(self):
        if len(self._segments) == 0:
            return
        prepoint = self._segments[0].start
        self._segments[0].start = None
        p = Path()
        subpaths = list(self.as_subpaths())
        for subpath in subpaths:
            subpath.reverse()
        for subpath in reversed(subpaths):
            p += subpath
        self._segments = p._segments
        self._segments[0].start = prepoint
        return self

    def subpath(self, index):
        subpaths = list(self.as_subpaths())
        return subpaths[index]

    def count_subpaths(self):
        subpaths = list(self.as_subpaths())
        return len(subpaths)

    def as_subpaths(self):
        last = 0
        for current, seg in enumerate(self):
            if current != last and isinstance(seg, Move):
                yield Subpath(self, last, current - 1)
                last = current
        yield Subpath(self, last, len(self) - 1)

    def as_points(self):
        """Returns the list of defining points within path"""
        for seg in self:
            for p in seg:
                if not isinstance(p, Point):
                    yield Point(p)
                else:
                    yield p

    @staticmethod
    def svg_d(segments, relative=None, smooth=None):
        if len(segments) == 0:
            return ""
        parts = []
        previous_segment = None
        p = Point(0)
        if smooth is None:
            override_smooth = False
            smooth_set_value = True
        else:
            override_smooth = True
            smooth_set_value = bool(smooth)
        if relative is not None:
            for segment in segments:
                if isinstance(segment, (Move, Line, Arc, Close)):
                    parts.append(segment.d(p, relative=relative))
                elif isinstance(segment, (CubicBezier, QuadraticBezier)):
                    if (override_smooth and smooth_set_value) or (
                        not override_smooth and segment.smooth
                    ):
                        parts.append(
                            segment.d(
                                p,
                                relative=relative,
                                smooth=segment.is_smooth_from(previous_segment),
                            )
                        )
                    else:
                        parts.append(segment.d(p, relative=relative, smooth=False))
                previous_segment = segment
                p = previous_segment.end
        else:
            for segment in segments:
                if isinstance(segment, (Move, Line, Arc, Close)):
                    parts.append(segment.d(p, relative=segment.relative))
                elif isinstance(segment, (CubicBezier, QuadraticBezier)):
                    if (override_smooth and smooth_set_value) or (
                        not override_smooth and segment.smooth
                    ):
                        parts.append(
                            segment.d(
                                p,
                                relative=segment.relative,
                                smooth=segment.is_smooth_from(previous_segment),
                            )
                        )
                    else:
                        parts.append(
                            segment.d(p, relative=segment.relative, smooth=False)
                        )
                previous_segment = segment
                p = previous_segment.end
        return " ".join(parts)

    def d(self, relative=None, transformed=True, smooth=None):
        path = self
        if transformed:
            path = abs(path)
        return Path.svg_d(path._segments, relative=relative, smooth=smooth)

    def segments(self, transformed=True):
        if transformed and not self.transform.is_identity():
            return [s * self.transform for s in self._segments]
        return self._segments

    def approximate_arcs_with_cubics(self, error=0.1):
        """
        Iterates through this path and replaces any Arcs with cubic bezier curves.
        """
        sweep_limit = tau * error
        for s in range(len(self) - 1, -1, -1):
            segment = self[s]
            if isinstance(segment, Arc):
                arc_required = int(ceil(abs(segment.sweep) / sweep_limit))
                self[s : s + 1] = list(segment.as_cubic_curves(arc_required))

    def approximate_arcs_with_quads(self, error=0.1):
        """
        Iterates through this path and replaces any Arcs with quadratic bezier curves.
        """
        sweep_limit = tau * error
        for s in range(len(self) - 1, -1, -1):
            segment = self[s]
            if isinstance(segment, Arc):
                arc_required = int(ceil(abs(segment.sweep) / sweep_limit))
                self[s : s + 1] = list(segment.as_quad_curves(arc_required))


########################
#  Path Subpath Object
########################


class Subpath:
    """
    Subpath is a Path-backed window implementation. It does not store a list of segments but rather
    stores a Path, start position, end position. When a function is called on a subpath, the result of
    those events is performed on the backing Path. When the backing Path is modified the behavior is
    undefined."""

    def __init__(self, path, start, end):
        self._path = path
        self._start = start
        self._end = end

    def __copy__(self):
        return Subpath(Path(self._path), self._start, self._end)

    def __getitem__(self, index):
        return self._path[self.index_to_path_index(index)]

    def __setitem__(self, index, value):
        self._path[self.index_to_path_index(index)] = value

    def __delitem__(self, index):
        del self._path[self.index_to_path_index(index)]
        self._end -= 1

    def __iadd__(self, other):
        if isinstance(other, str):
            p = Path(other)
            self._path[self._end : self._end] = p
        elif isinstance(other, Path):
            p = copy(other)
            self._path[self._end : self._end] = p
        elif isinstance(other, PathSegment):
            self._path.insert(self._end, other)
        else:
            return NotImplemented
        return self

    def __add__(self, other):
        if isinstance(other, (str, Path, PathSegment)):
            n = copy(self)
            n += other
            return n
        return NotImplemented

    def __radd__(self, other):
        if isinstance(other, str):
            path = Path(other)
            path.extend(map(copy, self._path))
            return path
        elif isinstance(other, PathSegment):
            path = Path(self)
            path.insert(0, other)
            return path
        else:
            return NotImplemented

    def __imul__(self, other):
        if isinstance(other, str):
            other = Matrix(other)
        if isinstance(other, Matrix):
            for e in self:
                e *= other
        return self

    def __mul__(self, other):
        if isinstance(other, (Matrix, str)):
            n = copy(self)
            n *= other
            return n
        return NotImplemented

    __rmul__ = __mul__

    def __iter__(self):
        class Iterator:
            def __init__(self, subpath):
                self.n = subpath._start - 1
                self.subpath = subpath

            def __next__(self):
                self.n += 1
                try:
                    if self.n > self.subpath._end:
                        raise StopIteration
                    return self.subpath._path[self.n]
                except IndexError:
                    raise StopIteration

            next = __next__

        return Iterator(self)

    def __len__(self):
        return self._end - self._start + 1

    def __str__(self):
        return self.d()

    def __repr__(self):
        return "Path(%s)" % (", ".join(repr(x) for x in self))

    def __eq__(self, other):
        if isinstance(other, str):
            return self.__eq__(Path(other))
        if not isinstance(other, (Path, Subpath)):
            return NotImplemented
        if len(self) != len(other):
            return False
        for s, o in zip(self, other):
            if not s == o:
                return False
        return True

    def __ne__(self, other):
        if not isinstance(other, (Path, Subpath, str)):
            return NotImplemented
        return not self == other

    def segments(self, transformed=True):
        path = self._path
        if transformed:
            return [
                s * path.transform for s in path._segments[self._start : self._end + 1]
            ]
        return path._segments[self._start : self._end + 1]

    def _numeric_index(self, index):
        if index < 0:
            return self._end + index + 1
        else:
            return self._start + index

    def index_to_path_index(self, index):
        if isinstance(index, slice):
            start = index.start
            stop = index.stop
            step = index.step
            if start is None:
                start = 0
            start = self._numeric_index(start)
            if stop is None:
                stop = len(self)
            stop = self._numeric_index(stop)
            return slice(start, stop, step)
        return self._numeric_index(index)

    def bbox(self):
        """returns a bounding box for the input Path"""
        segments = self._path._segments[self._start : self._end + 1]
        bbs = [seg.bbox() for seg in segments if not isinstance(Close, Move)]
        try:
            xmins, ymins, xmaxs, ymaxs = list(zip(*bbs))
        except ValueError:
            return None  # No bounding box items existed. So no bounding box.
        xmin = min(xmins)
        xmax = max(xmaxs)
        ymin = min(ymins)
        ymax = max(ymaxs)
        return xmin, ymin, xmax, ymax

    def d(self, relative=None, smooth=None):
        segments = self._path._segments[self._start : self._end + 1]
        return Path.svg_d(segments, relative=relative, smooth=None)

    def _reverse_segments(self, start, end):
        """Reverses segments between the given indexes in the subpath space."""
        segments = self._path._segments  # must avoid path validation.
        s = self.index_to_path_index(start)
        e = self.index_to_path_index(end)
        while s <= e:
            start_segment = segments[s]
            end_segment = segments[e]
            start_segment.reverse()
            if start_segment is not end_segment:
                end_segment.reverse()
                segments[s] = end_segment
                segments[e] = start_segment
            s += 1
            e -= 1
        start = self.index_to_path_index(start)
        end = self.index_to_path_index(end)
        self._path._validate_connection(start - 1, prefer_second=True)
        self._path._validate_connection(end)

    def reverse(self):
        size = len(self)
        if size == 0:
            return
        start = 0
        end = size - 1
        if isinstance(self[-1], Close):
            end -= 1
        if isinstance(
            self[0], Move
        ):  # Move remains in place but references next element.
            start += 1
        self._reverse_segments(start, end)
        if size > 1:
            if isinstance(self[0], Move):
                self[0].end = Point(self[1].start)
        last = self[-1]
        if isinstance(last, Close):
            last.reverse()
            if last.start != self[-2].end:
                last.start = Point(self[-2].end)
            if last.end != self[0].end:
                last.end = Point(self[0].end)
        return self
