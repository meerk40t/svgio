# svgio

`svgio` reads, writes, creates, updates, and performs geometric rendering of SVG files.

The goal of the project it to provide comparable access to SVG files to Javascript executed in a browser. 

This project uses a DOM Tree as a primary structure for loading of the SVG file. If the node loaded is a node within the SVG-Namespace the node is replaced with the namespace node, and provides the same functionality as a locally instanced object.

Many practical uses for svg are to convey geometric information. However, due to the complexity of the SVG spec there are a number of edge cases that low level parsing attempts cannot overcome. The goal of `svgio` is to merge the high fidelity parsing of projects like `svgelements` with the svg writing and generating abilities of `svgwrite`, and `cairo`. But, also allow the generation of rendered geometry like with `svg.path` or `svgpathtools`.  

This isn't intended to be a jack-of-all-trades but a correct implementation of how SVG is supposed to work rather than a single minor aspect.


# Status
This is a work in progress. It currently does some loading and bootstrapping and parsing. It currently writes the data out. Most of the needed components exist and work, such as Path. There's still some questions concerning the API. And the rendering from the functionality from the `svgelements` project isn't correctly implemented yet. And there's a lot of bugs, and generally many places that have nothing hooked up.

# Reading

The loading process is done with node bootstrapping to allow `Rect` tags to become `Rect` nodes within the tree. This works the way javascript's access to the DOM works, as a consequence this allows operations that should be permitted on a `Rect` object as well as rendering the objects into geometric shapes, and this allows direct modifications of the DOM Tree itself to reflected in the geometric data when rendered. This means that loading and building the DOM Tree are equivalent operations, and the SVG and XML as a whole, can be losslessly loaded and saved.

`svgio` reads SVG files and parses them into their relevant nodes:
```python
        q = io.StringIO(u'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 80 80">
                <g stroke="blue" id="group1">
                <circle id="circle1" cx="40" cy="40" r="35"/>
                </g>
                </svg>''')
        m = Document(q)
        group = m.get_element_by_id("group1")
        group.line(0, 0, 1, 1)
        print(m[0])
```

Prints:

```xml
<svg viewBox="0 0 80 80"><g id="group1" stroke="blue"><circle cx="40" cy="40" id="circle1" r="35" /><line x1="0" x2="1" y1="0" y2="1" /></g></svg>
```
The steps done by that small section of code:
* Parse a string based svg file
* Find an element by its ID
* Add in line() object
* Print the modified code.

# Editing

`svgio` permits you to edit the DOM with similar API to javascript and what is defined in the SVG spec. This includes looking up nodes by their ID.

# Saving

`svgio` permits you to save svg files edited svg files. The files can be sent to disk as well as printed as a string.

# Creating

`svgio` allows you to create SVG files from scratch and to save them to disk. For `svgwrite` like workflows.

# Rendering

`svgio` permits you to render the geometry of the DOM tree. This breaks down the particular tree nodes according to the svg/css rules, and gives you a series of relevant Path values. The rendered path values are the result of the current state of the DOM tree. So if you alter the transformation on a group, all items in the group will be altered within the final rendering. Rendering works much as it does in a Browser or other rendering engine. But rather than filled shapes produce correct Path commands from the current state of the Document Tree. 

# Basic

Due to the need to perform high fidelity rendering we include many standard and interacting elements from the SVG and CSS specs: Path, Matrix, Angle, Length, Color, Point and other SVG and CSS Elements. The SVG spec defines a variety of elements which generally interoperate.

The basic classes provide Matrix, Angle, Length, Color, Point and other SVG and CSS Elements. These classes can be used as part of the project and are used commonly throughout the project. However, they can also be used independently to perform CSS/SVG operations.

For example, the Color() class includes not only full CSS color parsing but also several useful utility functions and conversions for Color objects.

The Angle, Length, and Point classes are generally also able to be used for other projects.

# Nodes

The Nodes module provides the needed classes for the DOM tree, bootstrapping of nodes, modifications, and saving of nodes. This works somewhat like `svgwrite` if `svgwrite` also had the ability to read and modify the DOM nodes.

# Path

The Path module is similar to `svg.path` and `svgpathtools` projects and is actually derived from those projects.

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

