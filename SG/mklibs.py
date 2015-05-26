#
# Make an html list of libraries
#

import glob

def mklibs():
    libs = glob.glob('sglib*.py')
    libs.sort()
    for nb in libs:
        str = '<li><b>' + nb + '</b>'
        str += ' | <a href="' + nb + '">Download</a>'
        str += '\n<p>'
        print(str)
    print("</ul>")

mklibs()
