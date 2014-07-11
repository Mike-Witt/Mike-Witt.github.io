#
# Make an html list of ipython notebooks
#

import glob

def mklist1():
    for nb in glob.glob('*.ipynb'):
        str = '<li><b><a href="'
        str += 'http://nbviewer.ipython.org/url/mike-witt.github.io/SG/'
        str += nb + '"> | ' + nb
        str += '<a href="' + nb + '">Download</a></b>'
        str += '\n<p>'
        print str
mklist1()
