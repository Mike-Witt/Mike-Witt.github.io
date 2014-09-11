#
# Make an html list of ipython notebooks (for SR subdirectory)
#

import glob

def mklist1():
    notebooks = glob.glob('*.ipynb')
    notebooks.sort()
    for nb in notebooks:
        str = '<li><b><a href="'
        str += 'http://nbviewer.ipython.org/url/mike-witt.github.io/SG/SR_Lectures/'
        str += nb + '">' + nb + '</a></b>'
        str += ' | <a href="' + nb + '">Download</a>'
        str += '\n<p>'
        print(str)
    print("</ul>")
    print('<pre>\n\n</pre>')

mklist1()
