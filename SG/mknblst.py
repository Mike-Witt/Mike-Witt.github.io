#
# Make an html list of ipython notebooks
#

import glob

def mklist1():
    notebooks = glob.glob('0*.ipynb')
    notebooks += glob.glob('1*.ipynb')
    notebooks.sort()
    for nb in notebooks:
        str = '<li><b><a href="'
        str += 'http://nbviewer.ipython.org/url/mike-witt.github.io/SG/'
        str += nb + '">' + nb + '</a></b>'
        str += ' | <a href="' + nb + '">Download</a>'
        str += '\n<p>'
        print(str)
    print("</ul>")
    print('<pre>\n\n</pre>')

mklist1()
