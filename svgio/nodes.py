# -*- coding: ISO-8859-1 -*-

import re

from copy import copy
from math import (
    atan2,
    cos,
    sin,
    sqrt,
    tan,
)
from xml.etree.ElementTree import iterparse

try:
    from math import tau
except ImportError:
    from math import pi

    tau = pi * 2

from .constants import *
from .basic import *

"""
The path elements are derived from regebro's svg.path project ( https://github.com/regebro/svg.path ) with
some of the math from mathandy's svgpathtools project ( https://github.com/mathandy/svgpathtools ).

The goal is to provide svg like path objects and structures. The svg standard 1.1 and elements of 2.0 will
be used to provide much of the decisions within path objects. Such that if there is a question on
implementation if the SVG documentation has a methodology it should be used.

Though not required the SVGImage class acquires new functionality if provided with PIL/Pillow as an import
and the Arc can do exact arc calculations if scipy is installed.
"""

########################
#  Core DOM Objects
########################

class Node(list):
    def __init__(self, document, parent, attributes=None):
        super().__init__()
        self.document = document
        self.parent = parent
        self.attributes = attributes if attributes is not None else dict()

    def __repr__(self):
        name = self.__class__.__name__
        values = []
        for key in self.attributes:
            value = self.attributes.get(key)
            values.append("%s=%s" % (str(key), str(value)))
        params = ", ".join(values)
        return "%s(%s)" % (name, params)

    def parse(self, source):
        node = self
        document = self.document
        parent = self.parent

        for event, elem in iterparse(source, events=("start", "end", "start-ns")):
            if event == "start":
                tag = elem.tag
                try:
                    if tag.startswith("{http://www.w3.org/2000/svg"):
                        tag = tag[28:]  # Removing namespace. http://www.w3.org/2000/svg:
                except AttributeError:
                    pass

                attributes = elem.attrib
                # Create new node.
                parent = node  # parent is now previous node context
                siblings = parent  # Parent list is my list of children.

                # define this node.
                try:
                    type_node = BOOTSTRAP[tag]
                    node = type_node(document, parent, attributes)
                except KeyError:
                    node = Node(document, parent, attributes)

                siblings.append(node)  # siblings now includes this node.

                if SVG_ATTR_ID in attributes:  # If we have an ID, we save the node.
                    document.defs[attributes[SVG_ATTR_ID]] = node  # store node value in defs.
            elif event == "end":
                # event is 'end', pop values.
                node = node.parent   # Parent is now node.

    @property
    def first_child(self):
        """
        read-only property returns the node's first child in the tree, or null if the node has no children. If the node
        is a Document, it returns the first node in the list of its direct children.
        :return:
        """
        if len(self):
            return self[0]
        else:
            return None

    @property
    def is_connected(self):
        """
        read-only property of the Node interface returns a boolean indicating whether the node is connected (directly or
         indirectly) to the context object, for example the Document object in the case of the normal DOM, or the
         ShadowRoot in the case of a shadow DOM.
        :return:
        """
        return self.parent is not None

    @property
    def last_child(self):
        """
        read-only property returns the last child of the node. If its parent is an element, then the child is generally
        an element node, a text node, or a comment node. It returns null if there are no child elements.
        :return:
        """
        if len(self):
            return self[-1]
        else:
            return None

    @property
    def next_sibling(self):
        """
        read-only property returns the node immediately following the specified one in their parent's childNodes, or
        returns null if the specified node is the last child in the parent element.
        :return:
        """
        try:
            index = self.parent.index(self)
            return self.parent[index + 1]
        except (AttributeError, IndexError):
            return None

    @property
    def node_type(self):
        """
        property is an integer that identifies what the node is. It distinguishes different kind of nodes from each
        other, such as elements, text and comments.
        :return:
        """
        return self._type

    @property
    def owner_document(self):
        """
        read-only property of the Node interface returns the top-level document object of the node.
        :return:
        """
        return self.document

    @property
    def parent_element(self):
        """read-only property returns the DOM node's parent Element, or null if the node either has no parent, or its
         parent isn't a DOM Element."""
        return self.parent

    @property
    def parent_node(self):
        """
        read-only property returns the parent of the specified node in the DOM tree.
        :return:
        """
        return self.parent

    @property
    def previous_sibling(self):
        """
        read-only property returns the node immediately preceding the specified one in its parent's childNodes list, or
        null if the specified node is the first in that list.
        :return:
        """
        try:
            index = self.parent.index(self)
            if index == 0:
                return None
            return self.parent[index - 1]
        except (AttributeError, IndexError):
            return None

    @property
    def text_content(self):
        return None

    def append_child(self, child):
        """
        Method adds a node to the end of the list of children of a specified parent node. If the given child is a
        reference to an existing node in the document, appendChild() moves it from its current position to the new
        position (there is no requirement to remove the node from its parent node before appending it to some other
        node).
        :param child:
        :return:
        """
        raise NotImplementedError

    def clone_node(self):
        """method returns a duplicate of the node on which this method was called."""
        raise NotImplementedError

    def contains(self, node):
        """
        method returns a Boolean value indicating whether a node is a descendant of a given node, i.e. the node itself,
        one of its direct children (childNodes), one of the children's direct children, and so on.
        :param node:
        :return:
        """
        raise NotImplementedError

    def get_root_node(self):
        """method of the Node interface returns the context object's root, which optionally includes the shadow root if
        it is available."""
        raise NotImplementedError

    def has_child_nodes(self):
        """method returns a Boolean value indicating whether the given Node has child nodes or not."""
        return bool(len(self))

    def insert_before(self):
        """method inserts a node before a reference node as a child of a specified parent node."""
        raise NotImplementedError

    def is_default_namespace(self):
        """
        method accepts a namespace URI as an argument and returns a Boolean with a value of true if the namespace is the
         default namespace on the given node or false if not.
        :return:
        """
        raise NotImplementedError

    def is_equal_node(self, node):
        """
        method tests whether two nodes are equal. Two nodes are equal when they have the same type, defining
        characteristics (for elements, this would be their ID, number of children, and so forth), its attributes match,
        and so on. The specific set of data points that must match varies depending on the types of the nodes.
        :return:
        """
        return self == node

    def is_same_node(self, node):
        """
        method for Node objects tests whether two nodes are the same (that is, whether they reference the same object).
        :return:
        """
        return self is node

    def lookup_prefix(self):
        """
        method returns a DOMString containing the prefix for a given namespace URI, if present, and null if not. When
        multiple prefixes are possible, the result is implementation-dependent.
        :return:
        """
        raise NotImplementedError

    def normalize(self):
        """
        method puts the specified node and all of its sub-tree into a "normalized" form. In a normalized sub-tree, no
        text nodes in the sub-tree are empty and there are no adjacent text nodes.
        :return:
        """
        raise NotImplementedError

    def remove_child(self):
        """
        method removes a child node from the DOM and returns the removed node.
        :return:
        """
        raise NotImplementedError

    def replace_child(self):
        """
        method replaces a child node within the given (parent) node.
        :return:
        """
        raise NotImplementedError


class Element(Node):
    """
    Any element within the SVG namespace.

    If additional args exist these will be passed to property_by_args
    """

    def __init__(self,document, parent, attributes):
        Node.__init__(self, document, parent, attributes)

    @property
    def id(self):
        return self.attributes.get(SVG_ATTR_ID)

    @property
    def child_element_count(self):
        return len(self)

    @property
    def children(self):
        return self

    @property
    def class_list(self):
        return None

    @property
    def class_name(self):
        return None

    @property
    def client_height(self):
        return 0

    @property
    def client_width(self):
        return 0

    @property
    def client_top(self):
        return 0

    @property
    def client_left(self):
        return 0

    @property
    def current_style(self):
        return None

    @property
    def first_element_child(self):
        return self.first_child

    @property
    def last_element_child(self):
        return self.last_child

    @property
    def next_element_sibling(self):
        return self.next_sibling

    @property
    def previous_element_sibling(self):
        return self.next_sibling

    @property
    def shadow_root(self):
        return None

    def property_by_args(self, *args):
        pass

    def property_by_object(self, obj):
        self.id = obj.id
        self.values = dict(obj.values)

    def property_by_values(self, values):
        self.id = values.get(SVG_ATTR_ID)

    def set(self, key, value):
        self.values[key] = value
        return self


########################
#  SVG Element Types
########################

class ContainerElement:
    pass


class StructuralElement:
    pass


class DescriptiveElement:
    pass


class TextContentElement:
    pass


class Transformable:
    """Any element that is transformable and has a transform property."""

    def __init__(self, *args, **kwargs):
        self._length = None
        self._lengths = None
        self.transform = None
        self.apply = None

    def property_by_object(self, s):
        self.transform = Matrix(s.transform)
        self.apply = s.apply

    def property_by_values(self, values):
        self.transform = Matrix(values.get(SVG_ATTR_TRANSFORM, ""))
        self.apply = bool(values.get("apply", True))

    def __mul__(self, other):
        if isinstance(other, (Matrix, str)):
            n = copy(self)
            n *= other
            return n
        return NotImplemented

    __rmul__ = __mul__

    def __imul__(self, other):
        if isinstance(other, str):
            other = Matrix(other)
        if isinstance(other, Matrix):
            self.transform *= other
        return self

    def __abs__(self):
        """
        The absolute value is taken to be the actual shape transformed.
        :return: transformed version of the given shape.
        """
        m = copy(self)
        m.reify()
        return m

    def reify(self):
        """
        Realizes the transform to the attributes. Such that the attributes become actualized and the transform
        simplifies towards the identity matrix. In many cases it will become the identity matrix. In other cases the
        transformed shape cannot be represented through the properties alone. And shall keep those parts of the
        transform required preserve equivalency.

        The default method will be called by submethods but will only scale properties like stroke_width which should
        scale with the transform.
        """
        self._lengths = None
        self._length = None

    def render(self, **kwargs):
        """
        Renders the transformable by performing any required length conversion operations into pixels. The element
        will be the pixel-length form.
        """
        if self.transform is not None:
            self.transform.render(**kwargs)
        return self

    def bbox(self, transformed=True):
        """
        Returns the bounding box of the given object.

        :param transformed: whether this is the transformed bounds or default.
        :return:
        """
        raise NotImplementedError

    @property
    def rotation(self):
        if not self.apply:
            return Angle.degrees(0)
        prx = Point(1, 0)
        prx *= self.transform
        origin = Point(0, 0)
        origin *= self.transform
        return origin.angle_to(prx)


class GraphicalElement:
    """Any drawn element."""

    @property
    def stroke_opacity(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_STROKE_OPACITY)

    @property
    def stroke(self):
        values = Document.get_values(self)
        stroke = values.get(SVG_ATTR_STROKE)
        stroke = Color(stroke) if stroke is not None else None
        stroke_opacity = values.get(SVG_ATTR_STROKE_OPACITY)
        if (
            stroke_opacity is not None
            and self.stroke is not None
            and self.stroke.value is not None
        ):
            try:
                self.stroke.opacity = float(stroke_opacity)
            except ValueError:
                pass
        return stroke

    @property
    def fill_opacity(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_FILL_OPACITY)

    @property
    def fill(self):
        values = Document.get_values(self)
        fill = values.get(SVG_ATTR_FILL)
        fill = Color(fill) if fill is not None else None
        fill_opacity = values.get(SVG_ATTR_FILL_OPACITY)
        if (
            fill_opacity is not None
            and fill is not None
            and fill.value is not None
        ):
            try:
                fill.opacity = float(fill_opacity)
            except ValueError:
                pass
        return fill

    @property
    def stroke_width(self):
        # TODO: Stroke width percents need to calculate the viewbox and use that.
        # if isinstance(self.stroke_width, Length):
        #     width = kwargs.get("width", kwargs.get("relative_length"))
        #     height = kwargs.get("height", kwargs.get("relative_length"))
        #     try:
        #         del kwargs["relative_length"]
        #     except KeyError:
        #         pass
        #     self.stroke_width = self.stroke_width.value(
        #         relative_length=sqrt(width * width + height * height), **kwargs
        #     )
        #     # A percentage stroke_width is always computed as a percentage of the normalized viewBox diagonal length.

        try:
            if not self.apply:
                return self.stroke_width
            if self.stroke_width is not None:
                if (
                    hasattr(self, "values")
                    and SVG_ATTR_VECTOR_EFFECT in self.values
                    and SVG_VALUE_NON_SCALING_STROKE
                    in self.values[SVG_ATTR_VECTOR_EFFECT]
                ):
                    return self.stroke_width  # we are not to scale the stroke.
                width = self.stroke_width
                det = self.transform.determinant
                return width * sqrt(abs(det))
        except AttributeError:
            return self.stroke_width

    def is_degenerate(self):
        return False


########################
#  SVG Basic Shapes
########################


class Shape(Element, GraphicalElement, Transformable):
    """
    SVG Shapes are several SVG items defined in SVG 1.1 9.1
    https://www.w3.org/TR/SVG11/shapes.html

    These shapes are circle, ellipse, line, polyline, polygon, and path.

    All shapes have methods:
    d(relative, transform): provides path_d string for the shape.
    reify(): Applies transform of the shape to modify the shape attributes.
    render(): Ensure that the shape properties have real space values.
    bbox(transformed): Provides the bounding box for the given shape.

    All shapes must implement:
    __repr__(), with a call to _repr_shape()
    __copy__()

    All shapes have attributes:
    id: SVG ID attributes. (Element)
    transform: SVG Matrix to apply to this shape. (Transformable)
    apply: Determine whether transform should be applied. (Transformable)
    fill: SVG color of the shape fill. (GraphicObject)
    stroke: SVG color of the shape stroke. (GraphicObject)
    stroke_width: Stroke width of the stroke. (GraphicObject)
    """

    def __eq__(self, other):
        if not isinstance(other, Shape):
            return NotImplemented
        if self.fill != other.fill or self.stroke != other.stroke:
            return False
        first = self
        if not isinstance(first, Path):
            first = Path(first)
        second = other
        if not isinstance(second, Path):
            second = Path(second)
        return first == second

    def __ne__(self, other):
        if not isinstance(other, Shape):
            return NotImplemented
        return not self == other

    def __iadd__(self, other):
        if isinstance(other, Shape):
            return Path(self) + Path(other)
        return NotImplemented

    __add__ = __iadd__

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
        self.reify()
        return self

    def segments(self, transformed=True):
        """
        Returns PathSegments which correctly produce this shape.

        This should be implemented by subclasses.
        """
        raise NotImplementedError

    def d(self, relative=False, transformed=True):
        """
        Returns the path_d string of the shape.

        :param relative: Returns path_d in relative form.
        :param transformed: Return path_d, with applied transform.
        :return: path_d string
        """
        return Path(self.segments(transformed=transformed)).d(relative=relative)

    def bbox(self, transformed=True):
        """
        Get the bounding box for the given shape.
        """
        bbs = [
            seg.bbox()
            for seg in self.segments(transformed=transformed)
            if not isinstance(Close, Move)
        ]
        try:
            xmins, ymins, xmaxs, ymaxs = list(zip(*bbs))
        except ValueError:
            return None  # No bounding box items existed. So no bounding box.
        return min(xmins), min(ymins), max(xmaxs), max(ymaxs)


class Rect(Shape):
    """
    SVG Rect shapes are defined in SVG2 10.2
    https://www.w3.org/TR/SVG2/shapes.html#RectElement

    These have geometric properties x, y, width, height, rx, ry
    Geometric properties can be Length values.

    Rect(x, y, width, height)
    Rect(x, y, width, height, rx, ry)
    Rect(x, y, width, height, rx, ry, matrix)
    Rect(x, y, width, height, rx, ry, matrix, stroke, fill)

    Rect(dict): dictionary values read from svg.
    """

    def property_by_values(self, values):
        Shape.property_by_values(self, values)
        self.x = Length(values.get(SVG_ATTR_X, 0)).value()
        self.y = Length(values.get(SVG_ATTR_Y, 0)).value()
        self.width = Length(values.get(SVG_ATTR_WIDTH, 1)).value()
        self.height = Length(values.get(SVG_ATTR_HEIGHT, 1)).value()
        self.rx = Length(values.get(SVG_ATTR_RADIUS_X, None)).value()
        self.ry = Length(values.get(SVG_ATTR_RADIUS_Y, None)).value()

    def _validate_rect(self):
        """None is 'auto' for values."""
        rx = self.rx
        ry = self.ry
        if rx is None and ry is None:
            rx = ry = 0
        if rx is not None and ry is None:
            rx = Length(rx).value(relative_length=self.width)
            ry = rx
        elif ry is not None and rx is None:
            ry = Length(ry).value(relative_length=self.height)
            rx = ry
        elif rx is not None and ry is not None:
            rx = Length(rx).value(relative_length=self.width)
            ry = Length(ry).value(relative_length=self.height)
        if rx == 0 or ry == 0:
            rx = ry = 0
        else:
            try:
                rx = min(rx, self.width / 2.0)
                ry = min(ry, self.height / 2.0)
            except ValueError:
                pass  # If width is in inches and rx is in units, this is unsolvable without knowing the ppi
        self.rx = rx
        self.ry = ry

    @property
    def position(self):
        values = Document.get_values(self)
        x = values.get(SVG_ATTR_X)
        y = values.get(SVG_ATTR_Y)
        return Point(x,y)

    @property
    def x(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_X)

    @property
    def y(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_Y)

    @property
    def width(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_WIDTH)

    @property
    def height(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_HEIGHT)

    @property
    def rx(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_RADIUS_X)

    @property
    def ry(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_RADIUS_Y)

    def segments(self, transformed=True):
        """
        Rect decomposition is given in SVG 2.0 10.2

        Rect:
        * perform an absolute moveto operation to location (x,y);
        * perform an absolute horizontal lineto with parameter x+width;
        * perform an absolute vertical lineto parameter y+height;
        * perform an absolute horizontal lineto parameter x;
        * ( close the path)

        Rounded Rect:
        rx and ry are used as the equivalent parameters to the elliptical arc command,
        the x-axis-rotation and large-arc-flag are set to zero, the sweep-flag is set to one

        * perform an absolute moveto operation to location (x+rx,y);
        * perform an absolute horizontal lineto with parameter x+width-rx;
        * perform an absolute elliptical arc operation to coordinate (x+width,y+ry)
        * perform an absolute vertical lineto parameter y+height-ry;
        * perform an absolute elliptical arc operation to coordinate (x+width-rx,y+height)
        * perform an absolute horizontal lineto parameter x+rx;
        * perform an absolute elliptical arc operation to coordinate (x,y+height-ry)
        * perform an absolute vertical lineto parameter y+ry
        * perform an absolute elliptical arc operation with a segment-completing close path operation

        :param transformed: provide the reified version.
        :return: path_d of shape.
        """
        scooped = False
        x = self.x
        y = self.y
        width = self.width
        height = self.height
        if self.is_degenerate():
            return ()  # a computed value of zero for either dimension disables rendering.
        rx = self.rx
        ry = self.ry
        if not self._strict:
            if rx < 0 < width and ry < 0 < height:
                scooped = True
                rx = abs(rx)
                ry = abs(ry)
        if rx < 0 < width or ry < 0 < height:
            rx = 0
            ry = 0
        if rx == ry == 0:
            segments = (
                Move(None, (x, y)),
                Line((x, y), (x + width, y)),
                Line((x + width, y), (x + width, y + height)),
                Line((x + width, y + height), (x, y + height)),
                Close((x, y + height), (x, y)),
            )
        else:
            segments = (
                Move(None, (x + rx, y)),
                Line((x + rx, y), (x + width - rx, y)),
                Arc(
                    (x + width - rx, y),
                    (x + width, y + ry),
                    rx=rx,
                    ry=ry,
                    scooped=scooped,
                ),
                Line((x + width, y + ry), (x + width, y + height - ry)),
                Arc(
                    (x + width, y + height - ry),
                    (x + width - rx, y + height),
                    rx=rx,
                    ry=ry,
                    scooped=scooped,
                ),
                Line((x + width - rx, y + height), (x + rx, y + height)),
                Arc(
                    (x + rx, y + height),
                    (x, y + height - ry),
                    rx=rx,
                    ry=ry,
                    scooped=scooped,
                ),
                Line((x, y + height - ry), (x, y + ry)),
                Arc((x, y + ry), (x + rx, y), rx=rx, ry=ry, scooped=scooped),
                Close((x + rx, y), (x + rx, y)),
            )
        if not transformed or self.transform.is_identity():
            return segments
        else:
            return [s * self.transform for s in segments]

    def reify(self):
        """
        Realizes the transform to the shape properties.

        If the realized shape can be properly represented as a rectangle with an identity matrix
        it will be, otherwise the properties will approximate the implied values.

        Skewed and Rotated rectangles cannot be reified.
        """
        scale_x = self.transform.value_scale_x()
        scale_y = self.transform.value_scale_y()
        if scale_x * scale_y < 0:
            return self  # No reification of negative values, gives negative dims.
        translate_x = self.transform.value_trans_x()
        translate_y = self.transform.value_trans_y()
        if (
            self.transform.value_skew_x() == 0
            and self.transform.value_skew_y() == 0
            and scale_x != 0
            and scale_y != 0
        ):
            GraphicalElement.reify(self)
            Transformable.reify(self)
            self.x *= scale_x
            self.y *= scale_y
            self.x += translate_x
            self.y += translate_y
            self.transform *= Matrix.translate(-translate_x, -translate_y)
            self.rx = scale_x * self.rx
            self.ry = scale_y * self.ry
            self.width = scale_x * self.width
            self.height = scale_y * self.height
            self.transform *= Matrix.scale(1.0 / scale_x, 1.0 / scale_y)
        return self

    def render(self, **kwargs):
        Shape.render(self, **kwargs)
        width = kwargs.get("width", kwargs.get("relative_length"))
        height = kwargs.get("height", kwargs.get("relative_length"))
        try:
            del kwargs["relative_length"]
        except KeyError:
            pass
        if isinstance(self.x, Length):
            self.x = self.x.value(relative_length=width, **kwargs)
        if isinstance(self.y, Length):
            self.y = self.y.value(relative_length=height, **kwargs)
        if isinstance(self.width, Length):
            self.width = self.width.value(relative_length=width, **kwargs)
        if isinstance(self.height, Length):
            self.height = self.height.value(relative_length=height, **kwargs)
        if isinstance(self.rx, Length):
            self.rx = self.rx.value(relative_length=width, **kwargs)
        if isinstance(self.ry, Length):
            self.ry = self.ry.value(relative_length=height, **kwargs)
        return self

    def is_degenerate(self):
        return self.width == 0 or self.height == 0


class _RoundShape(Shape):

    def property_by_values(self, values):
        Shape.property_by_values(self, values)
        self.cx = Length(values.get(SVG_ATTR_CENTER_X)).value()
        self.cy = Length(values.get(SVG_ATTR_CENTER_Y)).value()
        self.rx = Length(values.get(SVG_ATTR_RADIUS_X)).value()
        self.ry = Length(values.get(SVG_ATTR_RADIUS_Y)).value()
        r = Length(values.get(SVG_ATTR_RADIUS, None)).value()
        if r is not None:
            self.rx = r
            self.ry = r
        else:
            if self.rx is None:
                self.rx = 1
            if self.ry is None:
                self.ry = 1
        center = values.get("center", None)
        if center is not None:
            self.cx, self.cy = Point(center)

        if self.cx is None:
            self.cx = 0
        if self.cy is None:
            self.cy = 0

    @property
    def rx(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_RADIUS_X)

    @property
    def ry(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_RADIUS_Y)

    @property
    def cx(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_CENTER_X)

    @property
    def cy(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_CENTER_Y)

    @property
    def center(self):
        values = Document.get_values(self)
        cx = values.get(SVG_ATTR_CENTER_X)
        cy = values.get(SVG_ATTR_CENTER_Y)
        return Point(cx, cy)

    def segments(self, transformed=True):
        """
        SVG path decomposition is given in SVG 2.0 10.3, 10.4.

        A move-to command to the point cx+rx,cy;
        arc to cx,cy+ry;
        arc to cx-rx,cy;
        arc to cx,cy-ry;
        arc with a segment-completing close path operation.

        Converts the parameters from an ellipse or a circle to a string for a
        Path object d-attribute"""
        original = self.apply
        self.apply = transformed
        path = Path()
        steps = 4
        step_size = tau / steps
        if (
            transformed
            and self.transform.value_scale_x() * self.transform.value_scale_y() < 0
        ):
            step_size = -step_size
        t_start = 0
        t_end = step_size
        # zero for either dimension, or a computed value of auto for both dimensions, disables rendering of the element.
        rx = self.rx
        ry = self.ry
        if self.is_degenerate():
            return ()
        center = self.center
        path.move((self.point_at_t(0)))
        for i in range(steps):
            path += Arc(
                self.point_at_t(t_start),
                self.point_at_t(t_end),
                center,
                rx=rx,
                ry=ry,
                rotation=self.rotation,
                sweep=step_size,
            )
            t_start = t_end
            t_end += step_size
        path.closed()
        self.apply = original
        return path.segments(transformed)

    def unit_matrix(self):
        """
        return the unit matrix which could would transform the unit circle into this ellipse.

        One of the valid parameterizations for ellipses is that they are all affine transforms of the unit circle.
        This provides exactly such a matrix.

        :return: matrix
        """
        m = Matrix()
        m.post_scale(self.rx, self.ry)
        m.post_rotate(self.rotation)
        center = self.center
        m.post_translate(center.x, center.y)
        return m

    def arc_t(self, t0, t1):
        """
        return the arc found between the given values of t on the ellipse.

        :param t0: t start
        :param t1: t end
        :return: arc
        """
        return Arc(
            self.point_at_t(t0),
            self.point_at_t(t1),
            self.center,
            rx=self.rx,
            ry=self.ry,
            rotation=self.rotation,
            sweep=t1 - t0,
        )

    def arc_angle(self, a0, a1, ccw=None):
        """
        return the arc found between the given angles on the ellipse.

        :param a0: start angle
        :param a1: end angle
        :param ccw: optional flag to force clockwise or counter-clockwise arc-angles, default is smaller angle
        :return: arc
        """
        if ccw is None:
            ccw = a0 > a1
        return Arc(
            self.point_at_angle(a0),
            self.point_at_angle(a1),
            self.center,
            rx=self.rx,
            ry=self.ry,
            rotation=self.rotation,
            ccw=ccw,
        )

    def point_at_angle(self, angle):
        """
        find the point on the ellipse from the center at the given angle.
        Note: For non-circular arcs this is different than point(t).

        :param angle: angle from center to find point
        :return: point found
        """
        a = self.rx
        b = self.ry
        if a == b:
            return self.point_at_t(angle)
        angle -= self.rotation
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
        if self.apply and not self.transform.is_identity():
            return self.center.angle_to(p)
        else:
            center = Point(self.cx, self.cy)
            return center.angle_to(p)

    def t_at_point(self, p):
        """
        find the t parameter to at the point.

        :param p: point
        :return: t parameter to the given point.
        """
        angle = self.angle_at_point(p)
        angle -= self.rotation
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
        rotation = self.rotation
        a = self.rx
        b = self.ry
        center = self.center
        cx = center.x
        cy = center.y
        cos_theta = cos(rotation)
        sin_theta = sin(rotation)
        cos_t = cos(t)
        sin_t = sin(t)
        px = cx + a * cos_t * cos_theta - b * sin_t * sin_theta
        py = cy + a * cos_t * sin_theta + b * sin_t * cos_theta
        return Point(px, py)

    def point(self, position, error=ERROR):
        """
        find the point that corresponds to given value [0,1].
        Where t=0 is the first point and t=1 is the final point.

        :param position: position value between 0,1 where value equals the amount through the shape
        :param error: error permitted in determining point value (unused for this shape)
        :return: point at t
        """
        return self.point_at_t(tau * position)

    def _ramanujan_length(self):
        a = self.rx
        b = self.ry
        if b > a:
            a, b = b, a
        h = (a - b) ** 2 / (a + b) ** 2
        return pi * (a + b) * (1 + (3 * h / (10 + sqrt(4 - 3 * h))))


class Ellipse(_RoundShape):
    """
    SVG Ellipse shapes are defined in SVG2 10.4
    https://www.w3.org/TR/SVG2/shapes.html#EllipseElement

    These have geometric properties cx, cy, rx, ry
    """

    def _name(self):
        return self.__class__.__name__


class Circle(_RoundShape):
    """
    SVG Circle shapes are defined in SVG2 10.3
    https://www.w3.org/TR/SVG2/shapes.html#CircleElement

    These have geometric properties cx, cy, r
    """

    @property
    def r(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_RADIUS)

    def _name(self):
        return self.__class__.__name__


class Line(Shape):
    """
    SVG Line shapes are defined in SVG2 10.5
    https://www.w3.org/TR/SVG2/shapes.html#LineElement

    These have geometric properties x1, y1, x2, y2

    These are called Line in SVG but that name is already used for Line(PathSegment)
    """


    @property
    def x1(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_X1)

    @property
    def y1(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_Y1)

    @property
    def x2(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_X2)

    @property
    def y2(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_Y2)

    def segments(self, transformed=True):
        """
        SVG path decomposition is given in SVG 2.0 10.5.

        perform an absolute moveto operation to absolute location (x1,y1)
        perform an absolute lineto operation to absolute location (x2,y2)

        :returns Path_d path for line.
        """

        start = Point(self.x1, self.y1)
        end = Point(self.x2, self.y2)
        if transformed:
            start *= self.transform
            end *= self.transform
        return Move(None, start), Line(start, end)


class _Polyshape(Shape):
    """Base form of Polygon and Polyline since the objects are nearly the same."""


    @property
    def points(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_POINTS)

    def segments(self, transformed=True):
        """
        Polyline and Polygon decomposition is given in SVG2. 10.6 and 10.7

        * perform an absolute moveto operation to the first coordinate pair in the list of points
        * for each subsequent coordinate pair, perform an absolute lineto operation to that coordinate pair.
        * (Polygon-only) perform a closepath command

        Note:  For a polygon/polyline made from n points, the resulting path will
        be composed of n lines (even if some of these lines have length zero).
        """
        # TODO: Parse points get segments.
        if len(self.points) != 0:
            return
        if points is None:
            self.points = list()
            return
        if isinstance(points, dict):
            if SVG_ATTR_POINTS in points:
                points = points[SVG_ATTR_POINTS]
            else:
                self.points = list()
                return
        try:
            if len(points) == 1:
                points = points[0]
        except TypeError:
            pass
        if isinstance(points, str):
            findall = REGEX_COORD_PAIR.findall(points)
            self.points = [Point(float(j), float(k)) for j, k in findall]
        elif isinstance(points, (list, tuple)):
            if len(points) == 0:
                self.points = list()
            else:
                first_point = points[0]
                if isinstance(first_point, (float, int)):
                    self.points = list(map(Point, zip(*[iter(points)] * 2)))
                elif isinstance(first_point, (list, tuple, complex, str, Point)):
                    self.points = list(map(Point, points))
        else:
            self.points = list()

        if self.transform.is_identity() or not transformed:
            points = self.points
        else:
            points = list(map(self.transform.point_in_matrix_space, self.points))
        if self.is_degenerate():
            return []
        segments = [Move(None, points[0])]
        last = points[0]
        for i in range(1, len(points)):
            current = points[i]
            segments.append(Line(last, current))
            last = current
        if isinstance(self, Polygon):
            segments.append(Close(last, points[0]))
        return segments


class Polyline(_Polyshape):
    """
    SVG Polyline shapes are defined in SVG2 10.6
    https://www.w3.org/TR/SVG2/shapes.html#PolylineElement

    These have geometric properties points
    """
    def _name(self):
        return self.__class__.__name__


class Polygon(_Polyshape):
    """
    SVG Polygon shapes are defined in SVG2 10.7
    https://www.w3.org/TR/SVG2/shapes.html#PolygonElement

    These have geometric properties points
    """
    def _name(self):
        return self.__class__.__name__


class Path(Shape):

    @property
    def d(self):
        values = Document.get_values(self)
        return values.get(SVG_ATTR_DATA)


########################
#  SVG Element Nodes
########################

class Group(ContainerElement,  StructuralElement, Element, Transformable):
    """
    Group Container element can have children.
    SVG 2.0 <g> are defined in:
    5.2. Grouping: the g element
    """
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)
        # list.__init__(self)
        # if len(args) >= 1:
        #     s = args[0]
        #     if isinstance(s, Group):
        #         self.extend(list(map(copy, s)))

    def __imul__(self, other):
        if isinstance(other, str):
            other = Matrix(other)
        if isinstance(other, Matrix):
            self.transform *= other
            for e in self:
                e *= other
        return self

    def render(self, **kwargs):
        Transformable.render(self, **kwargs)

    def __copy__(self):
        return Group(self)

    def property_by_object(self, s):
        Transformable.property_by_object(self, s)
        Element.property_by_object(self, s)

    def property_by_values(self, values):
        Transformable.property_by_values(self, values)
        Element.property_by_values(self, values)

    def select(self, conditional=None):
        """
        Finds all flattened subobjects of this group for which the conditional returns
        true.

        :param conditional: function taking element and returns True to include or False if exclude
        """
        if conditional is None:
            for subitem in self:
                yield subitem
                if isinstance(subitem, Group):
                    for s in subitem.select(conditional):
                        yield s
        else:
            for subitem in self:
                if conditional(subitem):
                    yield subitem
                if isinstance(subitem, Group):
                    for s in subitem.select(conditional):
                        yield s

    def reify(self):
        Transformable.reify(self)

    @staticmethod
    def union_bbox(elements, transformed=True):
        boundary_points = []
        for e in elements:
            if not hasattr(e, "bbox"):
                continue
            box = e.bbox(False)
            if box is None:
                continue
            top_left = (box[0], box[1])
            top_right = (box[2], box[1])
            bottom_left = (box[0], box[3])
            bottom_right = (box[2], box[3])
            if transformed:
                top_left = e.transform.point_in_matrix_space(top_left)
                top_right = e.transform.point_in_matrix_space(top_right)
                bottom_left = e.transform.point_in_matrix_space(bottom_left)
                bottom_right = e.transform.point_in_matrix_space(bottom_right)
            boundary_points.append(top_left)
            boundary_points.append(top_right)
            boundary_points.append(bottom_left)
            boundary_points.append(bottom_right)
        if len(boundary_points) == 0:
            return None
        xmin = min([e[0] for e in boundary_points])
        ymin = min([e[1] for e in boundary_points])
        xmax = max([e[0] for e in boundary_points])
        ymax = max([e[1] for e in boundary_points])
        return xmin, ymin, xmax, ymax

    def bbox(self, transformed=True):
        """
        Returns the bounding box of the given object.

        In the case of groups this is the union of all the bounding boxes of all bound children.

        Setting transformed to false, may yield unexpected results if subitems are transformed in non-uniform
        ways.

        :param transformed: bounding box of the properly transformed children.
        :return:
        """
        return Group.union_bbox(self.select(), transformed)


class ClipPath(Element):
    """
    clipPath elements are defined in svg 14.3.5
    https://www.w3.org/TR/SVG11/masking.html#ClipPathElement

    Clip paths conceptually define a 1 bit mask for images these are usually defined within
    def blocks and do not render themselves but rather are attached by IRI references to the
    """
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)
        self.unit_type = SVG_UNIT_TYPE_USERSPACEONUSE

    def property_by_object(self, s):
        Element.property_by_object(self, s)
        self.unit_type = s.unit_type

    def property_by_values(self, values):
        Element.property_by_values(self, values)
        self.unit_type = self.values.get(
            SVG_ATTR_CLIP_UNIT_TYPE, SVG_UNIT_TYPE_USERSPACEONUSE
        )


class Mask(ContainerElement):
    pass


class Pattern(ContainerElement, Element):
    def __init__(self, document, parent, tag, *args, **kwargs):
        super().__init__(document, parent, tag, *args, **kwargs)
        self.viewbox = None
        self.preserve_aspect_ratio = None
        self.x = None
        self.y = None
        self.width = None
        self.height = None
        self.href = None
        self.pattern_content_units = None  # UserSpaceOnUse default
        self.pattern_transform = None
        self.pattern_units = None

    def __int__(self):
        return 0

    @property
    def viewbox_transform(self):
        if self.viewbox is None:
            return ""
        return self.viewbox.transform(self)

    def property_by_object(self, s):
        Element.property_by_object(self, s)
        self.viewbox = s.viewbox
        self.preserve_aspect_ratio = s.preserve_aspect_ratio

        self.x = s.x
        self.y = s.y
        self.width = s.width
        self.height = s.height
        self.href = s.href
        self.pattern_content_units = s.pattern_contents_units
        self.pattern_transform = (
            Matrix(s.pattern_transform) if s.pattern_transform is not None else None
        )
        self.pattern_units = s.pattern_units

    def property_by_values(self, values):
        Element.property_by_values(self, values)
        if XLINK_HREF in values:
            self.href = values[XLINK_HREF]
        elif SVG_HREF in values:
            self.href = values[SVG_HREF]
        viewbox = values.get(SVG_ATTR_VIEWBOX)
        if viewbox is not None:
            self.viewbox = Viewbox(viewbox)
        if SVG_ATTR_PRESERVEASPECTRATIO in values:
            self.preserve_aspect_ratio = values[SVG_ATTR_PRESERVEASPECTRATIO]
        self.x = Length(values.get(SVG_ATTR_X, 0)).value()
        self.y = Length(values.get(SVG_ATTR_Y, 0)).value()
        self.width = Length(values.get(SVG_ATTR_WIDTH, "100%")).value()
        self.height = Length(values.get(SVG_ATTR_HEIGHT, "100%")).value()
        if SVG_ATTR_PATTERN_CONTENT_UNITS in values:
            self.pattern_content_units = values[SVG_ATTR_PATTERN_CONTENT_UNITS]
        if SVG_ATTR_PATTERN_TRANSFORM in values:
            self.pattern_transform = Matrix(values[SVG_ATTR_PATTERN_TRANSFORM])
        if SVG_ATTR_PATTERN_UNITS in values:
            self.pattern_units = values[SVG_ATTR_PATTERN_UNITS]

    def render(self, **kwargs):
        if self.pattern_transform is not None:
            self.pattern_transform.render(**kwargs)
        width = kwargs.get("width", kwargs.get("relative_length"))
        height = kwargs.get("height", kwargs.get("relative_length"))
        try:
            del kwargs["relative_length"]
        except KeyError:
            pass
        if isinstance(self.x, Length):
            self.x = self.x.value(relative_length=width, **kwargs)
        if isinstance(self.y, Length):
            self.y = self.y.value(relative_length=height, **kwargs)
        if isinstance(self.width, Length):
            self.width = self.width.value(relative_length=width, **kwargs)
        if isinstance(self.height, Length):
            self.height = self.height.value(relative_length=height, **kwargs)
        return self


class Text(TextContentElement, Element, GraphicalElement, Transformable):
    """
    SVG Text are defined in SVG 2.0 Chapter 11

    No methods are implemented to perform a text to path conversion.

    However, if such a method exists the assumption is that the results will be
    placed in the .path attribute, and functions like bbox() will check if such
    a value exists.
    """

    def __init__(self, *args, **kwargs):
        if len(args) >= 1:
            self.text = args[0]
        else:
            self.text = ""
        self.width = 0
        self.height = 0
        self.x = 0
        self.y = 0
        self.dx = 0
        self.dy = 0
        self.anchor = "start"  # start, middle, end.
        self.font_family = "san-serif"
        self.font_size = 16.0  # 16 point font 'normal'
        self.font_weight = 400.0  # Thin=100, Normal=400, Bold=700
        self.font_face = ""

        self.path = None
        Transformable.__init__(self, *args, **kwargs)
        GraphicalElement.__init__(self, *args, **kwargs)
        Element.__init__(self, *args, **kwargs)

    def __str__(self):
        parts = list()
        parts.append("'%s'" % self.text)
        parts.append("font_family=%s" % self.font_family)
        parts.append("anchor=%s" % self.anchor)
        parts.append("font_size=%d" % self.font_size)
        parts.append("font_weight=%s" % str(self.font_weight))
        return "Text(%s)" % (", ".join(parts))

    def __repr__(self):
        parts = list()
        parts.append("%s" % self.text)
        parts.append("font_family=%s" % self.font_family)
        parts.append("anchor=%s" % self.anchor)
        parts.append("font_size=%d" % self.font_size)
        parts.append("font_weight=%s" % str(self.font_weight))
        return "Text(%s)" % (", ".join(parts))

    def property_by_object(self, s):
        Transformable.property_by_object(self, s)
        GraphicalElement.property_by_object(self, s)
        Element.property_by_object(self, s)
        self.text = s.text
        self.x = s.x
        self.y = s.y
        self.width = s.width
        self.height = s.height
        self.dx = s.dx
        self.dy = s.dy
        self.anchor = s.anchor
        self.font_family = s.font_family
        self.font_size = s.font_size
        self.font_weight = s.font_weight
        self.font_face = s.font_face

    def parse_font(self, font):
        """
        CSS Fonts 3 has a shorthand font property which serves to provide a single location to define:
        ‘font-style’, ‘font-variant’, ‘font-weight’, ‘font-stretch’, ‘font-size’, ‘line-height’, and ‘font-family’

        font-style: normal | italic | oblique
        font-variant: normal | small-caps
        font-weight: normal | bold | bolder | lighter | 100 | 200 | 300 | 400 | 500 | 600 | 700 | 800 | 900
        font-stretch: normal | ultra-condensed | extra-condensed | condensed | semi-condensed | semi-expanded | expanded | extra-expanded | ultra-expanded
        font-size: <absolute-size> | <relative-size> | <length-percentage>
        line-height: '/' <‘line-height’>
        font-family: [ <family-name> | <generic-family> ] #
        generic-family:  ‘serif’, ‘sans-serif’, ‘cursive’, ‘fantasy’, and ‘monospace’
        """
        # https://www.w3.org/TR/css-fonts-3/#font-prop
        font_elements = list(*re.findall(REGEX_CSS_FONT, font))

        font_style = font_elements[0]
        font_variant = font_elements[1]
        font_weight = font_elements[2]
        font_stretch = font_elements[3]
        font_size = font_elements[4]
        line_height = font_elements[5]
        font_face = font_elements[6]
        font_family = font_elements[7]
        if len(font_weight) > 0:
            self.font_weight = self.parse_font_weight(font_weight)
        if len(font_size) > 0:
            self.font_size = Length(font_size).value()
        if len(font_face) > 0:
            if font_face.endswith(","):
                font_face = font_face[:-1]
            self.font_face = font_face

        if len(font_family) > 0:
            self.font_family = font_family

    def parse_font_weight(self, weight):
        if weight == "bold":
            return 700
        if weight == "normal":
            return 400
        try:
            return int(weight)
        except KeyError:
            return 400

    def property_by_values(self, values):
        Transformable.property_by_values(self, values)
        GraphicalElement.property_by_values(self, values)
        self.anchor = values.get(SVG_ATTR_TEXT_ANCHOR, self.anchor)
        self.font_face = values.get(SVG_ATTR_FONT_FACE)
        self.font_family = values.get(SVG_ATTR_FONT_FAMILY, self.font_family)
        self.font_size = Length(values.get(SVG_ATTR_FONT_SIZE, self.font_size)).value()
        self.font_weight = values.get(SVG_ATTR_FONT_WEIGHT, self.font_weight)
        font = values.get(SVG_ATTR_FONT, None)
        if font is not None:
            self.parse_font(font)
        self.text = values.get(SVG_TAG_TEXT, self.text)
        self.x = Length(values.get(SVG_ATTR_X, self.x)).value()
        self.y = Length(values.get(SVG_ATTR_Y, self.y)).value()
        self.dx = Length(values.get(SVG_ATTR_DX, self.dx)).value()
        self.dy = Length(values.get(SVG_ATTR_DY, self.dy)).value()

    def reify(self):
        GraphicalElement.reify(self)
        Transformable.reify(self)

    def render(self, **kwargs):
        GraphicalElement.render(self, **kwargs)
        Transformable.render(self, **kwargs)
        width = kwargs.get("width", kwargs.get("relative_length"))
        height = kwargs.get("height", kwargs.get("relative_length"))
        try:
            del kwargs["relative_length"]
        except KeyError:
            pass
        if isinstance(self.x, Length):
            self.x = self.x.value(relative_length=width, **kwargs)
        if isinstance(self.y, Length):
            self.y = self.y.value(relative_length=height, **kwargs)
        if isinstance(self.dx, Length):
            self.dx = self.dx.value(relative_length=width, **kwargs)
        if isinstance(self.dy, Length):
            self.dy = self.dy.value(relative_length=height, **kwargs)
        return self

    def __copy__(self):
        return Text(self)

    def bbox(self, transformed=True):
        """
        Get the bounding box for the given text object.
        """
        if self.path is not None:
            return (self.path * self.transform).bbox(transformed=True)
        width = self.width
        height = self.height
        xmin = self.x
        ymin = self.y - height
        xmax = self.x + width
        ymax = self.y
        if not hasattr(self, "anchor") or self.anchor == "start":
            pass
        elif self.anchor == "middle":
            xmin -= width / 2
            xmax -= width / 2
        elif self.anchor == "end":
            xmin -= width
            xmax -= width
        if transformed:
            p0 = self.transform.transform_point([xmin, ymin])
            p1 = self.transform.transform_point([xmin, ymax])
            p2 = self.transform.transform_point([xmax, ymin])
            p3 = self.transform.transform_point([xmax, ymax])
            xmin = min(p0[0], p1[0], p2[0], p3[0])
            ymin = min(p0[1], p1[1], p2[1], p3[1])
            xmax = max(p0[0], p1[0], p2[0], p3[0])
            ymax = max(p0[1], p1[1], p2[1], p3[1])
        return xmin, ymin, xmax, ymax


class TSpan(Text):
    pass


class Image(Element, GraphicalElement, Transformable):
    """
    SVG Images are defined in SVG 2.0 12.3

    This class is called SVG Image rather than image as a guard against many Image objects
    which are quite useful and would be ideal for reading the linked or contained data.
    """

    def __init__(self, *args, **kwargs):
        self.url = None
        self.data = None
        self.viewbox = None
        self.preserve_aspect_ratio = None
        self.x = None
        self.y = None
        self.width = None
        self.height = None

        self.image = None
        self.image_width = None
        self.image_height = None
        Transformable.__init__(self, *args, **kwargs)
        GraphicalElement.__init__(self, *args, **kwargs)
        Element.__init__(
            self, *args, **kwargs
        )  # Dataurl requires this be processed first.

        if self.url is not None:
            if self.url.startswith("data:image/"):
                # Data URL
                from base64 import b64decode

                if self.url.startswith("data:image/png;base64,"):
                    self.data = b64decode(self.url[22:])
                elif self.url.startswith("data:image/jpg;base64,"):
                    self.data = b64decode(self.url[22:])
                elif self.url.startswith("data:image/jpeg;base64,"):
                    self.data = b64decode(self.url[23:])
                elif self.url.startswith("data:image/svg+xml;base64,"):
                    self.data = b64decode(self.url[26:])

    def property_by_object(self, s):
        Element.property_by_object(self, s)
        Transformable.property_by_object(self, s)
        GraphicalElement.property_by_object(self, s)
        self.url = s.url
        self.data = s.data
        self.viewbox = s.viewbox
        self.preserve_aspect_ratio = s.preserve_aspect_ratio

        self.x = s.x
        self.y = s.y
        self.width = s.width
        self.height = s.height

        self.image = s.image
        self.image_width = s.image_width
        self.image_height = s.image_height

    def property_by_values(self, values):
        Element.property_by_values(self, values)
        Transformable.property_by_values(self, values)
        GraphicalElement.property_by_values(self, values)
        if XLINK_HREF in values:
            self.url = values[XLINK_HREF]
        elif SVG_HREF in values:
            self.url = values[SVG_HREF]
        viewbox = values.get(SVG_ATTR_VIEWBOX)
        if viewbox is not None:
            self.viewbox = Viewbox(viewbox)
        if SVG_ATTR_PRESERVEASPECTRATIO in values:
            self.preserve_aspect_ratio = values[SVG_ATTR_PRESERVEASPECTRATIO]
        self.x = Length(values.get(SVG_ATTR_X, 0)).value()
        self.y = Length(values.get(SVG_ATTR_Y, 0)).value()
        self.width = Length(values.get(SVG_ATTR_WIDTH, "100%")).value()
        self.height = Length(values.get(SVG_ATTR_HEIGHT, "100%")).value()
        if "image" in values:
            self.image = values["image"]
            self.image_width, self.image_height = self.image.size

    def render(self, **kwargs):
        GraphicalElement.render(self, **kwargs)
        Transformable.render(self, **kwargs)
        width = kwargs.get("width", kwargs.get("relative_length"))
        height = kwargs.get("height", kwargs.get("relative_length"))
        try:
            del kwargs["relative_length"]
        except KeyError:
            pass
        if isinstance(self.x, Length):
            self.x = self.x.value(relative_length=width, **kwargs)
        if isinstance(self.y, Length):
            self.y = self.y.value(relative_length=height, **kwargs)
        if isinstance(self.width, Length):
            self.width = self.width.value(relative_length=width, **kwargs)
        if isinstance(self.height, Length):
            self.height = self.height.value(relative_length=height, **kwargs)
        return self

    def __copy__(self):
        """
        Copy of SVGImage. This will not copy the .image subobject in a deep manner
        since it's optional that that object will exist or not. As such if using PIL it would
        be required to either say self.image = self.image.copy() or call .load() again.
        """
        return Image(self)

    @property
    def viewbox_transform(self):
        if self.viewbox is None:
            return ""
        return self.viewbox.transform(self)

    def load(self, directory=None):
        try:
            from PIL import Image

            if self.data is not None:
                self.load_data()
            elif self.url is not None:
                self.load_file(directory)
            self.set_values_by_image()
        except ImportError:
            pass

    def load_data(self):
        try:
            # This code will not activate without PIL/Pillow installed.
            from PIL import Image

            if self.data is not None:
                from io import BytesIO

                self.image = Image.open(BytesIO(self.data))
            else:
                return
        except ImportError:
            # PIL/Pillow not found, decoding data is most we can do.
            pass

    def load_file(self, directory):
        try:
            # This code will not activate without PIL/Pillow installed.
            from PIL import Image

            if self.url is not None:
                try:
                    self.image = Image.open(self.url)
                except IOError:
                    try:
                        if directory is not None:
                            from os.path import join

                            relpath = join(directory, self.url)
                            self.image = Image.open(relpath)
                    except IOError:
                        return
        except ImportError:
            # PIL/Pillow not found, decoding data is most we can do.
            pass

    def set_values_by_image(self):
        if self.image is None:
            return  # No image to set values by.
        self.image_width = self.image.width
        self.image_height = self.image.height
        self.viewbox = Viewbox(
            "0 0 %d %d" % (self.image_width, self.image_height),
            self.preserve_aspect_ratio,
        )
        self.render(width=self.image_width, height=self.image_height)
        self.transform = Matrix(self.viewbox_transform) * self.transform

    def bbox(self, transformed=True):
        """
        Get the bounding box for the given image object
        """
        if self.image_width is None or self.image_height is None:
            p = Point(0, 0)
            p *= self.transform
            return p.x, p.y, p.x, p.y
        width = self.image_width
        height = self.image_height
        if transformed:
            p = (
                Point(0, 0) * self.transform,
                Point(width, 0) * self.transform,
                Point(width, height) * self.transform,
                Point(0, height) * self.transform,
            )
        else:
            p = (Point(0, 0), Point(width, 0), Point(width, height), Point(0, height))
        x_vals = list(s.x for s in p)
        y_vals = list(s.y for s in p)
        min_x = min(x_vals)
        min_y = min(y_vals)
        max_x = max(x_vals)
        max_y = max(y_vals)
        return min_x, min_y, max_x, max_y


class Desc(DescriptiveElement, Element):

    def __init__(self, values, desc=None):
        self.desc = desc
        Element.__init__(self, **values)

    def property_by_object(self, obj):
        Element.property_by_object(self, obj)
        self.desc = obj.desc

    def property_by_values(self, values):
        Element.property_by_values(self, values)
        if SVG_TAG_DESC in values:
            self.desc = values[SVG_TAG_DESC]


class Title(DescriptiveElement, Element):
    def __init__(self, values, title=None):
        self.title = title
        Element.__init__(self, **values)

    def property_by_object(self, obj):
        Element.property_by_object(self, obj)
        self.title = obj.title

    def property_by_values(self, values):
        Element.property_by_values(self, values)
        if SVG_TAG_TITLE in values:
            self.title = values[SVG_TAG_TITLE]


class Metadata(DescriptiveElement, Element):
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)


class SVG(Group, ContainerElement, StructuralElement):
    """
    SVG Document and Parsing.

    SVG is the SVG main object and also the embedded SVGs within it. It's a subtype of Group. The SVG has a viewbox,
    and parsing methods which can be used if given a stream, path, or svg string.
    """

    def __init__(self, document, parent, attributes):
        Group.__init__(self, document, parent, attributes)
        self.x = None
        self.y = None
        self.width = None
        self.height = None
        self.viewbox = None

    @property
    def implicit_position(self):
        if not self.apply:
            return Point(self.x, self.y)
        point = Point(self.x, self.y)
        point *= self.transform
        return point

    @property
    def implicit_x(self):
        if not self.apply:
            return self.x
        return self.implicit_position[0]

    @property
    def implicit_y(self):
        if not self.apply:
            return self.y
        return self.implicit_position[1]

    @property
    def implicit_width(self):
        if not self.apply:
            return self.width
        p = Point(self.width, 0)
        p *= self.transform
        origin = Point(0, 0)
        origin *= self.transform
        return origin.distance_to(p)

    @property
    def implicit_height(self):
        if not self.apply:
            return self.height
        p = Point(0, self.height)
        p *= self.transform
        origin = Point(0, 0)
        origin *= self.transform
        return origin.distance_to(p)

    def property_by_object(self, s):
        Group.property_by_object(self, s)
        self.x = s.x
        self.y = s.y
        self.width = s.width
        self.height = s.height
        self.viewbox = Viewbox(s.viewbox) if s.viewbox is not None else None

    def property_by_values(self, values):
        Group.property_by_values(self, values)
        self.x = Length(values.get(SVG_ATTR_X, 0)).value()
        self.y = Length(values.get(SVG_ATTR_Y, 0)).value()
        self.width = Length(values.get(SVG_ATTR_WIDTH, "100%")).value()
        self.height = Length(values.get(SVG_ATTR_HEIGHT, "100%")).value()
        viewbox = values.get(SVG_ATTR_VIEWBOX)
        par = values.get(SVG_ATTR_PRESERVEASPECTRATIO)
        self.viewbox = Viewbox(viewbox, par) if viewbox is not None else None

    def render(self, **kwargs):
        Group.render(self, **kwargs)
        width = kwargs.get("width", kwargs.get("relative_length"))
        height = kwargs.get("height", kwargs.get("relative_length"))
        try:
            del kwargs["relative_length"]
        except KeyError:
            pass
        self.width = Length(self.width).value(relative_length=width, **kwargs)
        self.height = Length(self.height).value(relative_length=height, **kwargs)
        self.x = Length(self.x).value(relative_length=width, **kwargs)
        self.y = Length(self.y).value(relative_length=height, **kwargs)

    def elements(self, conditional=None):
        if conditional is None or conditional(self):
            yield self
        for q in self.select(conditional):
            yield q

    @property
    def viewbox_transform(self):
        if self.viewbox is None:
            return ""
        return self.viewbox.transform(self)


class Defs(ContainerElement, StructuralElement, Element):
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)


class Symbol(ContainerElement, StructuralElement, Element):
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)


class Marker(ContainerElement, Element):
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)


class Use(StructuralElement, Element):
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)


class Switch(ContainerElement, Element):
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)


class Style(Element):
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)


class Script(Element):
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)


class Hyperlink(ContainerElement, Element):
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)


class LinearGradient(Element):
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)


class RadialGradient(Element):
    def __init__(self, document, parent, attributes):
        Element.__init__(self, document, parent, attributes)


########################
#  Parsing XML -> DOM
########################

def _shadow_iter(tag, elem, children):
    yield tag, "start", elem
    try:
        for t, e, c in children:
            for shadow_tag, shadow_event, shadow_elem in _shadow_iter(t, e, c):
                yield shadow_tag, shadow_event, shadow_elem
    except ValueError:
        """
        Strictly speaking it is possible to reference use from other use objects. If this is an infinite loop
        we should not block the rendering. Just say we finished. See: W3C, struct-use-12-f
        """
        pass
    yield tag, "end", elem


def parse(
    source,
    ppi=DEFAULT_PPI,
    width=1000,
    height=1000,
    color="black",
    transform=None,
    context=None,
):
    """
    Parses the SVG file. All attributes are things which the SVG document itself could not be aware of, such as
    the real size of pixels and the size of the viewport (as opposed to the viewbox).

    :param source: Source svg file or stream.
    :param ppi: How many physical pixels per inch are there in this view.
    :param width: The physical width of the viewport
    :param height: The physical height of the viewport
    :param color: the `currentColor` value from outside the current scope.
    :param transform: Any required transformations to be pre-applied to this document
    :param context: Any existing document context.
    :return:
    """
    clip = 0
    root = context
    styles = {}
    stack = []

    values = {
        SVG_ATTR_COLOR: color,
        SVG_ATTR_FILL: "black",
        SVG_ATTR_STROKE: "none",
    }

    if transform is not None:
        values[SVG_ATTR_TRANSFORM] = transform

    for tag, event, elem in SVG._use_structure_parse(source):
        """
        SVG element parsing parses the job compiling any parsed elements into their compiled object forms.
        """
        # print(event, elem)
        if event == "start":
            stack.append((context, values))
            if (
                SVG_ATTR_DISPLAY in values
                and values[SVG_ATTR_DISPLAY] == SVG_VALUE_NONE
            ):
                continue  # Values has a display=none. Do not render anything. No Shadow Dom.
            current_values = values
            values = {}
            values.update(current_values)  # copy of dictionary

            # Non-propagating values.
            if SVG_ATTR_PRESERVEASPECTRATIO in values:
                del values[SVG_ATTR_PRESERVEASPECTRATIO]
            if SVG_ATTR_VIEWBOX in values:
                del values[SVG_ATTR_VIEWBOX]
            if SVG_ATTR_ID in values:
                del values[SVG_ATTR_ID]
            if SVG_ATTR_CLIP_PATH in values:
                del values[SVG_ATTR_CLIP_PATH]

            attributes = elem.attrib  # priority; lowest
            attributes[SVG_ATTR_TAG] = tag

            # Split any Style block elements into parts; priority medium
            style = ""
            if "*" in styles:  # Select all.
                style += styles["*"]
            if tag in styles:  # selector type
                style += styles[tag]
            if SVG_ATTR_ID in attributes:  # Selector id #id
                svg_id = attributes[SVG_ATTR_ID]
                css_tag = "#%s" % svg_id
                if css_tag in styles:
                    if len(style) != 0:
                        style += ";"
                    style += styles[css_tag]
            if SVG_ATTR_CLASS in attributes:  # Selector class .class
                for svg_class in attributes[SVG_ATTR_CLASS].split(" "):
                    css_tag = ".%s" % svg_class
                    if css_tag in styles:
                        if len(style) != 0:
                            style += ";"
                        style += styles[css_tag]
                    css_tag = "%s.%s" % (
                        tag,
                        svg_class,
                    )  # Selector type/class type.class
                    if css_tag in styles:
                        if len(style) != 0:
                            style += ";"
                        style += styles[css_tag]
            # Split style element into parts; priority highest
            if SVG_ATTR_STYLE in attributes:
                style += attributes[SVG_ATTR_STYLE]

            # Process style tag left to right.
            for equate in style.split(";"):
                equal_item = equate.split(":")
                if len(equal_item) == 2:
                    key = str(equal_item[0]).strip()
                    value = str(equal_item[1]).strip()
                    attributes[key] = value
            if (
                SVG_ATTR_FILL in attributes
                and attributes[SVG_ATTR_FILL] == SVG_VALUE_CURRENT_COLOR
            ):
                if SVG_ATTR_COLOR in attributes:
                    attributes[SVG_ATTR_FILL] = attributes[SVG_ATTR_COLOR]
                else:
                    attributes[SVG_ATTR_FILL] = values[SVG_ATTR_COLOR]

            if (
                SVG_ATTR_STROKE in attributes
                and attributes[SVG_ATTR_STROKE] == SVG_VALUE_CURRENT_COLOR
            ):
                if SVG_ATTR_COLOR in attributes:
                    attributes[SVG_ATTR_STROKE] = attributes[SVG_ATTR_COLOR]
                else:
                    attributes[SVG_ATTR_STROKE] = values[SVG_ATTR_COLOR]

            if SVG_ATTR_TRANSFORM in attributes:
                # If transform is already in values, append the new value.
                if SVG_ATTR_TRANSFORM in values:
                    attributes[SVG_ATTR_TRANSFORM] = (
                        values[SVG_ATTR_TRANSFORM]
                        + " "
                        + attributes[SVG_ATTR_TRANSFORM]
                    )
                else:
                    attributes[SVG_ATTR_TRANSFORM] = attributes[SVG_ATTR_TRANSFORM]

            # All class and attribute properties are compiled.

            values.update(attributes)
            values[SVG_STRUCT_ATTRIB] = attributes
            if (
                SVG_ATTR_DISPLAY in values
                and values[SVG_ATTR_DISPLAY] == SVG_VALUE_NONE
            ):
                continue  # If the attributes flags our values to display=none, stop rendering.
            if SVG_NAME_TAG == tag:
                # The ordering for transformations on the SVG object are:
                # explicit transform, parent transforms, attribute transforms, viewport transforms
                s = SVG(values)
                s.render(ppi=ppi, width=width, height=height)
                height, width = s.width, s.height
                if s.viewbox is not None:
                    try:
                        if s.height == 0 or s.width == 0:
                            return s
                        viewport_transform = s.viewbox_transform
                    except ZeroDivisionError:
                        # The width or height was zero.
                        # https://www.w3.org/TR/SVG11/struct.html#SVGElementWidthAttribute
                        # "A value of zero disables rendering of the element."
                        return s  # No more parsing will be done.

                    if SVG_ATTR_TRANSFORM in values:
                        # transform on SVG element applied as if svg had parent with transform.
                        values[SVG_ATTR_TRANSFORM] += " " + viewport_transform
                    else:
                        values[SVG_ATTR_TRANSFORM] = viewport_transform
                    width, height = s.viewbox.width, s.viewbox.height
                if context is None:
                    stack[-1] = (context, values)
                if context is not None:
                    context.append(s)
                context = s
                if root is None:
                    root = s
            elif SVG_TAG_GROUP == tag:
                s = Group(values)
                context.append(s)
                context = s
                s.render(ppi=ppi, width=width, height=height)
            elif SVG_TAG_DEFS == tag:
                s = Group(values)
                context = s  # Non-Rendered
                s.render(ppi=ppi, width=width, height=height)
            elif SVG_TAG_CLIPPATH == tag:
                s = ClipPath(values)
                context = s  # Non-Rendered
                s.render(ppi=ppi, width=width, height=height)
                clip += 1
            elif SVG_TAG_PATTERN == tag:
                s = Pattern(values)
                context = s  # Non-rendered
                s.render(ppi=ppi, width=width, height=height)
            elif tag in (
                SVG_TAG_PATH,
                SVG_TAG_CIRCLE,
                SVG_TAG_ELLIPSE,
                SVG_TAG_LINE,  # Shapes
                SVG_TAG_POLYLINE,
                SVG_TAG_POLYGON,
                SVG_TAG_RECT,
                SVG_TAG_IMAGE,
            ):
                try:
                    if SVG_TAG_PATH == tag:
                        s = Path(values)
                    elif SVG_TAG_CIRCLE == tag:
                        s = Circle(values)
                    elif SVG_TAG_ELLIPSE == tag:
                        s = Ellipse(values)
                    elif SVG_TAG_LINE == tag:
                        s = Line(values)
                    elif SVG_TAG_POLYLINE == tag:
                        s = Polyline(values)
                    elif SVG_TAG_POLYGON == tag:
                        s = Polygon(values)
                    elif SVG_TAG_RECT == tag:
                        s = Rect(values)
                    else:  # SVG_TAG_IMAGE == tag:
                        s = Image(values)
                except ValueError:
                    continue
                s.render(ppi=ppi, width=width, height=height)
                if reify:
                    s.reify()
                if s.is_degenerate():
                    continue
                context.append(s)
            elif tag in (
                SVG_TAG_STYLE,
                SVG_TAG_TEXT,
                SVG_TAG_DESC,
                SVG_TAG_TITLE,
                SVG_TAG_TSPAN,
            ):
                # <style>, <text>, <desc>, <title>
                continue
            else:
                s = Element(values)  # SVG Unknown object return as element.
                context.append(s)

            # Assign optional linked properties.
            try:
                clip_path_url = s.values.get(SVG_ATTR_CLIP_PATH, None)
                if clip_path_url is not None:
                    clip_path = root.get_element_by_url(clip_path_url)
                    s.clip_path = clip_path
            except AttributeError:
                pass
            if clip != 0:
                try:
                    clip_rule = s.values.get(SVG_ATTR_CLIP_RULE, SVG_RULE_NONZERO)
                    if clip_rule is not None:
                        s.clip_rule = clip_rule
                except AttributeError:
                    pass
            if SVG_ATTR_ID in attributes and root is not None:
                root.objects[attributes[SVG_ATTR_ID]] = s
        elif event == "end":  # End event.
            # The iterparse spec makes it clear that internal text data is undefined except at the end.
            s = None
            if tag in (
                SVG_TAG_TEXT,
                SVG_TAG_TSPAN,
                SVG_TAG_DESC,
                SVG_TAG_TITLE,
                SVG_TAG_STYLE,
            ):
                attributes = elem.attrib
                if SVG_ATTR_ID in attributes and root is not None:
                    root.objects[attributes[SVG_ATTR_ID]] = s
            if tag in (SVG_TAG_TEXT, SVG_TAG_TSPAN):
                s = Text(values, text=elem.text)
                s.render(ppi=ppi, width=width, height=height)
                if reify:
                    s.reify()
                context.append(s)
            elif SVG_TAG_DESC == tag:
                s = Desc(values, desc=elem.text)
                context.append(s)
            elif SVG_TAG_TITLE == tag:
                s = Title(values, title=elem.text)
                context.append(s)
            elif SVG_TAG_STYLE == tag:
                assignments = list(re.findall(REGEX_CSS_STYLE, elem.text))
                for key, value in assignments:
                    key = key.strip()
                    value = value.strip()
                    for selector in key.split(","):  # Can comma select subitems.
                        sel = selector.strip()
                        if sel not in styles:
                            styles[sel] = value
                        else:
                            if not styles[sel].endswith(";"):
                                styles[sel] += ";"
                            styles[sel] += value
            elif SVG_TAG_CLIPPATH == tag:
                clip -= 1
            if s is not None:
                # Assign optional linked properties.
                try:
                    clip_path_url = s.values.get(SVG_ATTR_CLIP_PATH, None)
                    if clip_path_url is not None:
                        clip_path = root.get_element_by_url(clip_path_url)
                        s.clip_path = clip_path
                except AttributeError:
                    pass
                if clip != 0:
                    try:
                        clip_rule = s.values.get(
                            SVG_ATTR_CLIP_RULE, SVG_RULE_NONZERO
                        )
                        if clip_rule is not None:
                            s.clip_rule = clip_rule
                    except AttributeError:
                        pass

            context, values = stack.pop()
        elif event == "start-ns":
            if elem[0] != SVG_ATTR_DATA:
                values[elem[0]] = elem[1]
    return root


class Document(Node):
    def __init__(self, source=None):
        super().__init__(self, None, {})  # document, parent, attributes
        self.css_sheets = None
        self.defs = {}
        self.pixels_per_inch = DEFAULT_PPI
        self.width = 1000
        self.height = 1000
        self.color = "black"
        self.transform = None
        if source is not None:
            self.parse(source)

    @staticmethod
    def get_ancestors(node):
        yield node
        while node.parent is not None:
            yield node.parent
            node = node.parent

    @staticmethod
    def get_values(node):
        values = dict()
        for n in reversed(list(Document.get_ancestors(node))):
            values.update(n.attributes)
        return values

    def get_element_by_id(self, id):
        return self.defs.get(id)

    def get_element_by_url(self, url):
        for _id in REGEX_IRI.findall(url):
            return self.get_element_by_id(_id)

    def load(self, filename):
        raise NotImplementedError

    def save(self, pretty=False, indent=2):
        raise NotImplementedError

    def saveas(self, filename, pretty=False, indent=2):
        raise NotImplementedError

    def add_stylesheet(self, Sheet):
        raise NotImplementedError

    def render(self):
        pass

    def get_xml(self):
        return str(self)


BOOTSTRAP = {
    "svg": SVG,
    "group": Group,
    "defs": Defs,
    "symbol": Symbol,
    "marker": Marker,
    "use": Use,
    "style": Style,
    "script": Script,
    "a": Hyperlink,
    "text": Text,
    "tspan": TSpan,
    "path": Path,
    "line": Line,
    "rect": Rect,
    "circle": Circle,
    "ellipse": Ellipse,
    "polyline": Polyline,
    "polygon": Polygon,
    "image": Image,
    "title": Title,
    "desc": Desc,
    "metadata": Metadata,
    "clip-path": ClipPath,
    "mask": Mask,
    "lineargradient": LinearGradient,
    "radialgradient": RadialGradient,
    "pattern": Pattern,
}
