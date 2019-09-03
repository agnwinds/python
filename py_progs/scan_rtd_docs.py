import docutils.frontend
import docutils.nodes
import docutils.parsers.rst
import docutils.utils


class LinkCheckerVisitor(docutils.nodes.GenericNodeVisitor):
    def default_visit(self, node):
        # Pass all other nodes through.
        print('NODE:', node)
        print(node.first_child_matching_class(('title')))



# Parse the file into a document with the rst parser.
fileobj = open('../docs/sphinx/source/input/parameters/central_object/Central_object.rst')
default_settings = docutils.frontend.OptionParser(
    components=(docutils.parsers.rst.Parser,)).get_default_values()
document = docutils.utils.new_document(fileobj.name, default_settings)
parser = docutils.parsers.rst.Parser()
parser.parse(fileobj.read(), document)

# Visit the parsed document with our link-checking visitor.
visitor = LinkCheckerVisitor(document)
document.walk(visitor)
