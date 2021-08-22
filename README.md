# svgio

`svgio` reads, writes, modifies, builds, and performs geometric rendering on SVG files. 

`svgio` does high fidelity SVG parsing into bootstrapped DOM-like nodes to facilitate reading, writing and editing of SVG files while also providing high fidelity geometric rendering from that structure.

The key difference in this project is the use of DOM Tree as a primary structure and geometric rendering as a process that is done to the DOM Tree structure. The loading process is done with node bootstrapping to allow `Rect` tags to become `Rect` nodes within the tree. This works the way javascript's access to the DOM works, as a consequence this allow operations that would be permitted on a `Rect` object as well as rendering the objects into geometric shapes, and this allows direct modifications of the DOM Tree itself to reflect on the geometric data when rendered. This means that we can build the DOM Tree itself and thus save the loaded document.

Due to the need to perform high fidelity rendering we include many standard and interacting elements from the SVG and CSS specs: Path, Matrix, Angle, Length, Color, Point and other SVG and CSS Elements. The SVG spec defines a variety of elements which generally interoperate.

# Progress
Currently nothing is done.

# License

This module is under a MIT License.
https://github.com/meerk40t/svgio/blob/master/LICENSE

# Installing
`pip install svgio`

Then in a script:

`from svgio import *`

# Requirements

None.

# Compatibility

`svgio` is compatible with Python 3+.

This module remains somewhat backwards compatible with `svg.path`, and with many parts of `svgelements` as well as intending to copy significant portions of `svgwrite`'s api.

# Philosophy

The goal of this project is to provide SVG input, output, and generation capabilities, elements and structures. This project should robustly conform to the SVG standard 1.1 and 2.0 SVG-Spec. If there is a question on implementation and the SVG documentation has a methodology, that is the preferred methodology. If the SVG spec says one thing, and `svgio` does something else, that is a bug.

`svgio` should conform to the SVG Conforming Interpreter class (2.5.4. Conforming SVG Interpreters):

>An SVG interpreter is a program which can parse and process SVG document fragments. Examples of SVG interpreters are server-side transcoding tools or optimizer (e.g., a tool which converts SVG content into modified SVG content) or analysis tools (e.g., a tool which extracts the text content from SVG content, or a validity checker).

# Acknowledgments

The Path element of this project is based in part on the `regebro/svg.path` ( https://github.com/regebro/svg.path ) project. It is also may be based, in part, on some elements of `mathandy/svgpathtools` ( https://github.com/mathandy/svgpathtools ).

