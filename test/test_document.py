import unittest
from svgio import *


class TestDocument(unittest.TestCase):

    def test_document(self):
        q = io.StringIO(u'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 80 80">
        <circle stroke="blue" id="circle1" cx="40" cy="40" r="35"/>
        </svg>''')
        m = Document(q)
        circle1 = m.get_element_by_id('circle1')
        print(circle1)
        self.assertEqual(circle1.fill, None)
        self.assertEqual(circle1.stroke, "blue")

    def test_document_inherit(self):
        q = io.StringIO(u'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 80 80">
        <g stroke="blue">
        <circle id="circle1" cx="40" cy="40" r="35"/>
        </g>
        </svg>''')
        m = Document(q)
        circle1 = m.get_element_by_id('circle1')
        print(circle1)
        self.assertEqual(circle1.fill, None)
        self.assertEqual(circle1.stroke, "blue")

    def test_document_build(self):
        q = io.StringIO(u'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 80 80">
                <g stroke="blue" id="group1">
                <circle id="circle1" cx="40" cy="40" r="35"/>
                </g>
                </svg>''')
        m = Document(q)
        group = m.get_element_by_id("group1")
        group.line(0, 0, 1, 1)
        print(m[0])

    def test_document_css_class(self):
        q = io.StringIO(u'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 80 80">
        <defs>
        <style>.cls-1,.cls-2{fill:none;stroke-miterlimit:10;}.cls-1{stroke:blue;}.cls-2{stroke:red;}</style>
        </defs>
        <circle class="cls-1" id="circle1" cx="40" cy="40" r="35"/>
        </svg>''')
        m = Document(q)
        circle1 = m.get_element_by_id('circle1')
        self.assertEqual(circle1.fill, None)
        self.assertEqual(circle1.stroke, "blue")