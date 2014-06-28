#
#   s g l i b . p y
#
# Python library for the Sunday Group
#
#
# Contents
#
# 1. Basic Setup (Imports and such)
#
# 2. Utility Functions
#    Globally defined functions (and other globals variable?)
#
# 3. Linear Algebra Routines
#
# 4. Problem Generator Code
#    4.1 Vector problems
#    4.2 Matrix problems
#    4.3 Change of basis problems
#    4.4 Electron Spin (XYZ Only. Problems.)
#
# 5. Electron Spin (Arbitrary direction. Tutorial.)
#    Uses XYZ stuff from 4.4 Electron Spin
#    Includes 3D plotting code which should probably be elsewhere
#

#####################################################################
#                                                                   #
#                       1. Basic Setup                              #
#                                                                   #
#####################################################################

from random import randint
import numpy as np
import sympy as sy
from sympy import latex, Matrix, Rational, sqrt
from sympy.physics.quantum import TensorProduct as TP
sy.init_printing()
Sqrt = sqrt
matrix = Matrix
i = sy.I

#####################################################################
#                                                                   #
#                       2. Utility Functions                        #
#                                                                   #
#####################################################################

# "megasimp"
# See: https://groups.google.com/forum/#!topic/sympy/As1Nct87ieI
#
# Ideally, I guess I'd like to make a function that executes
# something like this on each part of a composite objects.
# And, possibly, put it together with (a better version of "myltx")

from sympy import simplify, expand, fu, powsimp, sqrtdenest
from sympy.strategies.tree import greedy, brute
funcs = [simplify, expand, fu, powsimp, sqrtdenest]
objective = lambda x: len(str(x))  # minimize string length
megasimp = greedy((funcs, funcs), objective)

# font_fsize is for Print()
# This is not in "points." I don't know what the "measure" is.
font_size = 4
def set_font_size(size):
    global font_size
    font_size = size
#
# Print()
#
# Notice this "Print" has a caplital "P" to distunguish it from the
# regular python "print" statements. The reason for this function is
# that I can't figure out how to use the "normal" print to do latex.
# I imagine there is some "right" way to do it, but I haven't been
# able to figure it out.
#
def Print(thing):
    global font_size
    from IPython.core.display import HTML
    from IPython.core.display import display
    display(HTML('<font size=%s>'%font_size+thing+'</font>'))

# Doesn't simpy have a better way to do this?
def canonical_complex(z):
    re, im = z.as_real_imag()
    return( re+im*i)

# frac(): Enter a fraction (OBSOLETE)
# This is no longer needed. You should be able to say:
#
#   from __future__ import division
#
# in your notebook to make expressions like 1/2 work.
#
def frac(n,d):
    if type(n) == type(sqrt(2)): return(n/d)
    elif type(n) == type(i): return(n/d)
    if type(d) == type(sqrt(2)): return(n/d)
    elif type(d) == type(i): return(n/d)
    else: return( Rational(n,d) )

#-------------------------------------------------------------------#
#                                                                   #
#           Format class and sg_print function                      #
#                                                                   #
#-------------------------------------------------------------------#

class format:
    # Format things, depending on the "exact" flag. If we are not
    # doing exact math then limit the number of digits displayed by
    # converting each element to mpmath and back. My 28 April 2014
    # question on the sympy list appears to indicate that there is
    # no "easy" way to do this.
    #
    # I SHOULD MAKE SOMETHING A BIT MORE GENERAL THAN THIS AND PUT IT
    # IN THE IPYTHON MISC SUBROUTINES IN PROJECT/IPYTHON.
    #
    def __init__(self, exact=True, ndigs=3):
        self.exact = exact
        self.ndigs=ndigs

    def fix(self, n):
        from sympy import im, sympify
        from sympy.mpmath import mpmathify

        if im(n) != 0: n = complex(n)
        else: n = float(n)
        n = mpmathify(n)
        n = sympify(n) 
        return(n)

    def do_format(self, M):

        from numpy import shape
        from sympy.mpmath import mp, mpmathify

        if self.exact: return(myltx(M))
        mp.dps = self.ndigs

        # For a scalar
        if shape(M) == (): return(myltx(self.fix(M)))

        # For a one-dim list
        the_shape = shape(M)
        shape_of_shape = shape(the_shape)
        if shape_of_shape == (1,):
            for row in range(the_shape[0]):
                M[row] = self.fix(M[row])
            return(myltx(M))

        # For a list
        if type(M) == type([1,2,3]):
            rows, cols = shape(M)
            for row in range(rows):
                for col in range(cols):
                    M[row][col] = self.fix(M[row][col])
            return(myltx(M)) 

        # For a matrix
        # Here I need to make a copy. Otherwise we'll clobber the
        # original matrix, which in some cases would be bad. Note that
        # the same problem *may* exsist for lists above, and I just
        # haven't run into it yet. But lists don't appear to have a
        # "copy" method?
        Mc = M.copy()
        rows, cols = shape(M)
        for row in range(rows):
            for col in range(cols):
                Mc[row,col] = self.fix(M[row,col])
        return(myltx(Mc)) 

def sg_print(obj, exact=True, ndigs=3):
    fm = format(exact, ndigs).do_format
    Print('$%s$'%fm(obj))

#-------------------------------------------------------------------#
#                                                                   #
#                   Creating various objects                        #
#                                                                   #
#-------------------------------------------------------------------#

def col(*l): return( Matrix(l) )

def row(*l): return( Matrix(l).transpose() )

def mat(*l):
    from sympy import sympify, Matrix, zeros
    N = sqrt(len(l))
    if type(N) != type(sympify(2)):
        print('I can\'t make a square matrix out of %s numbers'%N)
        return(0)
    M = zeros(N)
    for row in range(N):
        for col in range(N):
            M[row,col] = l[N*row+col]
    return(M)

#-------------------------------------------------------------------#
#                                                                   #
#                       find_eigenvectors                           #
#                       print_eigenvectors                          #
#                                                                   #
#-------------------------------------------------------------------#

def find_eigenvectors(O, V=False):
    # If it's one of our "basic" operators, return the eigenvectors
    # explicitly (so that we get the same "version" that we're
    # using in class.
    if O == sigma_x: return( [-1, 1], [mX, pX] )
    if O == sigma_y: return( [-1, 1], [mY, pY] )
    if O == sigma_z: return( [-1, 1], [mZ, pZ] )
    evals = []
    evecs = []
    evs = O.eigenvects()
    for ev in evs:
        evals += [ ev[0], ]
        v = ev[2][0].normalized()
        #v.simplify()
        for n in range(len(v)): v[n] = megasimp(v[n])
        evecs += [ v, ]
    return(evals, evecs)

def print_eigenvectors(O, exact=True, ndigs=3):
    fm = format(exact, ndigs).do_format
    Print(r'The operator is: $%s$' %fm(O))
    eval, evec = find_eigenvectors(O)
    for n in range(len(eval)):
        Print('eigenvalue: $%s$, eigenvector: $%s$'%(fm(eval[n]),fm(evec[n])))

#-------------------------------------------------------------------#
#                                                                   #
#   Construct a measurement operator from eigenvalues and vectors   #
#                                                                   #
#-------------------------------------------------------------------#

def construct_measurement_operator(evals, evecs, Debug=False):
    from sympy import Matrix, zeros
    dim = len(evals)
    D = zeros(dim)
    V = zeros(dim)
    for m in range(dim):
        #print('m=%s'%m)
        D[m,m] = evals[m]
        for n in range(dim):
            #print('n=%s'%n)
            V[n,m] = evecs[m][n]
    Vi = V.inv()
    M = V*D*Vi
    if Debug:
        Print('$%s %s %s = %s$'%(myltx(V), myltx(D), myltx(Vi), myltx(M)))
    return(M)

#-------------------------------------------------------------------#
#                                                                   #
#                       "myltx" code                                #
#                                                                   #
#-------------------------------------------------------------------#

#
# This is "partially successful" so I guess I'll leave it here until
# I either get "parsing" working or else convert all the problem types
# to produce their own latex code.
#
    
def myltx(obj, V=False):
    # This outer function handles matrices (which is what our 'vectors'
    # are) or anything else that can be changed into a list.
    if V: print('myltx(): obj is %s'%obj)
    try:
        gcd = sy.gcd(list(obj))
        if V: print('gcd=%s, type is: %s'%(gcd, type(gcd)))
        if type(gcd) != type(sqrt(2)): return(sy.latex(obj))
        obj = obj*gcd
        return( myltx(1/gcd) + myltx_frac(obj) )
    except Exception, e:
        if V: print str(e)
        return(myltx_frac(obj))

def myltx_frac(frac):
    # If type is not sympy "mul" return original
    if type(frac) != type(1/sqrt(2)): return(sy.latex(frac))

    # if 'args' fails then return original, but print a
    # message to look into what happened.
    try: args = frac.args
    except Exception, e:
        print('args failed with: %s'%e) 
        print('Please inform Mike')
        return(sy.latex(frac))

    # We need to have two args
    if len(args) != 2: return(sy.latex(frac))
    one = frac.args[0]; two = frac.args[1]
    # One needs to be a fraction and the other a power.
    # Hopefully, the fraction always comes first!
    # Note: Sometimes we get type 1/2 instead of rational :-)
    half = sy.Rational(1,2); third = sy.Rational(1,3)
    if type(one) != type(half) and type(one) != type(third):
        return(sy.latex(frac))
    if type(two) != type(sqrt(2)):
        return(sy.latex(frac))
    # We've passed all the tests. Create the latex we want.
    num, den = one.as_numer_denom()
    den = den/two
    return( r'\frac{%s}{%s}'%(sy.latex(num), sy.latex(den)) )

#####################################################################
#                                                                   #
#                   3. Linear Algebra Routines                      #
#                                                                   #
#####################################################################

def vector(*args):
    v = []
    for element in args: v += [ element, ]
    return Matrix(v)

def length(v):
    return( v.norm() )

norm = length
absolute_value = length

def inner_product(v1, v2):
    if len(v1) != len(v2):
        print('Vectors must be the same dimension')
        return(None)
    prod = v1.transpose().conjugate() * v2
    # Sympy returns a 1x1 matrix. We want a scalar
    prod = prod[0]
    # Try everything possible to simplify it
    prod = megasimp(prod)
    return( prod )

def projection(basis_vector, vector):
    return( inner_product(basis_vector, vector)*basis_vector )

#
# Graphics - should separate these from LA?
#

import pylab
from matplotlib import pyplot
from IPython.core.pylabtools import figsize
        
def draw_plane(n_points=100, domain=(-4, 4, -4, 4)):
    # Think of this, for now, as "creating the x axis"
    x_axis = np.linspace(domain[0], domain[1], n_points)
    y_axis = np.linspace(domain[2], domain[3], n_points)
    zeros = []
    for n in range(n_points): zeros += [0,]
    
    # Then we plot the two functions on the same axis
    pylab.figure(figsize(6,6))
    pylab.plot(x_axis, zeros, color='black')
    pylab.plot(zeros, y_axis, color='black')
    xlimits = pylab.xlim(domain[0], domain[1])
    ylimits = pylab.ylim(domain[2], domain[3])

def draw_point(x,y, color='black'):
    plot([float(x)], [float(y)], marker='o', color=color, markeredgecolor=color)

def draw_vector(v, color='black', scale=1, HS=.07):
    if np.shape(v) != (2,1):
        print('I can only plot a 2D vector')
        return
    pyplot.arrow(0, 0, float(v[0]), float(v[1]), head_width=HS*scale,
        head_length=HS*scale, color=color, length_includes_head=True)
    
def draw_projections(b1, b2, v, bcolor='black', vcolor='green'):
    v_prime = vector( ip(b1,v), ip(b2,v) )
    v_prime.simplify()
    dom = max(1, abs(v[0]), abs(v[1]), abs(v_prime[0]), abs(v_prime[1]))
    dom = float(dom)*1.2
    Print('The vector $%s$, and its projection on the basis: $%s,%s$'
        %(myltx(v), myltx(b1), myltx(b2)))
    try:
        draw_plane(domain=[-dom, dom, -dom, dom])
        draw_vector(b1, color=bcolor, scale=dom)
        draw_vector(b2, color=bcolor, scale=dom)
        draw_vector(v, color=vcolor, scale=dom)
        draw_vector( projection(b1,v), color='light'+vcolor, scale=dom)
        draw_vector( projection(b2,v), color='light'+vcolor, scale=dom)
    except:
        pylab.gcf().clf()
        Print('I couldn\'t draw the picture. Presumably because there')
        Print('were complex values involved. But in any event ...') 
    Print('The vector in this basis is: $%s$'%myltx(v_prime))

#####################################################################
#                                                                   #
#                   4. Problem Generator Code                       #
#                                                                   #
#####################################################################

#
# These are global because I'll use them in more than one class
#

# Generate a random complex number
def rcn(min_int, max_int):
    real_part = randint(min_int, max_int)
    imag_part = randint(min_int, max_int)
    return( real_part + imag_part*i)

# Generate a random number that might be either real or complex
def rnum(min_int, max_int, complex=False):
    if complex: return( rcn(min_int, max_int) )
    else: return( randint(min_int, max_int) )

# Generate a fraction. If complex is set, then the numerator
# will be a complex number.
def rfrac(min_int, max_int, complex=False):
    n = rnum(min_int, max_int, complex)
    d = rnum(min_int, max_int)
    if complex: return( (n/d).factor() )
    else: return( Rational(n,d) )

#-------------------------------------------------------------------#
#                                                                   #
#                   4.1 Vector Problems                             #
#                                                                   #
#-------------------------------------------------------------------#

class vector_problem:
    def __init__(self, max_dim=3, max_int=5, complex=False):
        self.max_dim = max_dim
        self.max_int = max_int
        self.complex = complex

    def check(self, answer):
        if self.answer == answer: print('Correct')
        else:
            print('The correct answer is:')
            Print('$%s$'%myltx(self.answer))

    def gTranspose(self):
        v = []
        for n in range(self.dim):
            v += [ rnum(0, self.max_int, self.complex), ]
        v = Matrix(v)
        # Decide whether it's a column or a row
        if randint(0,1) == 1: v = v.transpose()
        if self.complex==True: Print(r'$%s^{\dagger}$' %myltx(v))
        else: Print('$%s^T$' %myltx(v))
        self.answer = v.transpose().conjugate()
    
    def gLength(self):
        v = []
        for n in range(self.dim):
            v += [ rnum(0, self.max_int, self.complex), ]
        v = Matrix(v)
        # Decide whether it's a column or a row
        if randint(0,1) == 1: v = v.transpose()
        Print('Find the length of: $%s$' %myltx(v))
        self.answer = v.norm()
    
    def gNormalize(self):
        v = []
        for n in range(self.dim):
            v += [ rnum(0, self.max_int, self.complex), ]
        v = Matrix(v)
        # Decide whether it's a column or a row
        if randint(0,1) == 1: v = v.transpose()
        Print('Normalize: $%s$' %myltx(v))
        self.answer = v/v.norm()   
    
    def new(self, type=None):
        # Choose the dimension of the problem
        self.dim = randint(2, self.max_dim)
        # Choose the type of problem
        if type == None: ptype = randint(1, 6)
        else: ptype = type
        # Transpose problem
        if ptype == 1:
            self.gTranspose()
            return
        if ptype == 5:
            self.gLength()
            return
        if ptype == 6:
            self.gNormalize()
            return
        # Otherwise generate two vectors and choose *, +, -
        v1 = []; v2 = []
        for n in range(self.dim):
            v1 += [ rnum(0, self.max_int, self.complex), ]
            v2 += [ rnum(0, self.max_int, self.complex), ]
        v1 = Matrix(v1); v2 = Matrix(v2)
        if ptype == 2:
            # Inner product. "Bra" the first vector.
            v1 = v1.transpose().conjugate()
            Print('$%s %s$'%(myltx(v1), myltx(v2)))
            self.answer = (v1 * v2)[0].simplify()
        # For + or - we might want to switch to rows
        if randint(0,1) == 1:
            v1 = v1.transpose()
            v2 = v2.transpose()
        if ptype == 3:
            Print('$%s + %s$'%(myltx(v1), myltx(v2)))
            self.answer = v1 + v2
        if ptype == 4:
            Print('$%s - %s$'%(myltx(v1), myltx(v2)))
            self.answer = v1 - v2

#-------------------------------------------------------------------#
#                                                                   #
#                   4.2 Matrix Problems                             #
#                                                                   #
#-------------------------------------------------------------------#

class matrix_problem:
    def __init__(self, max_dim=3, max_int=5, complex=False):
        self.max_dim = max_dim
        self.max_int = max_int
        self.complex = complex

    def new(self):
        type = randint(1, 7)
        if type == 1: self.Add()
        if type == 2: self.Add()
        if type == 3: self.Mult()
        if type == 4: self.Mult()
        if type == 5: self.Mult()
        if type == 6: self.Mult()
        if type == 7: self.Trans()

    # Addition / Subtraction
    def Add(self):
        dim = randint(2, self.max_dim)

        m1 = self.gen_mat(dim, dim)
        m2 = self.gen_mat(dim, dim)
        # Choose either addition or subtraction
        if randint(0,1) == 1:
            self.answer = m1 + m2
            Print('$%s + %s$'%(myltx(m1), myltx(m2)))
        else:
            self.answer = m1 - m2
            Print('$%s - %s$'%(myltx(m1), myltx(m2)))

    # Multiplication
    def Mult(self):
        nrows = randint(1, self.max_dim)
        ncols = randint(1, self.max_dim)
        # We don't want the rows and columns both to be one
        if nrows == 1 and ncols == 1: ncols = 2
        m1 = self.gen_mat(nrows, ncols)
        nrows = ncols
        ncols = randint(1, self.max_dim)
        if nrows == 1 and ncols == 1: ncols = 2
        m2 = self.gen_mat(nrows, ncols)
        self.answer = m1 * m2
        self.answer.simplify()
        Print('$%s %s$'%(myltx(m1), myltx(m2)))

    # Transpose
    def Trans(self):
        nrows = randint(2, self.max_dim)
        ncols = randint(2, self.max_dim)
        m = self.gen_mat(nrows, ncols)
        self.answer = m.transpose().conjugate()
        if self.complex: Print(r'$%s^{\dagger}$'%(myltx(m)))
        else: Print('$%s^T$'%(myltx(m)))

    def gen_mat(self, rows, cols):
        mat = []
        for r in range(rows):
            row = []
            for c in range(cols):
                #element = randint(1, self.max_int)
                element = rnum(0, self.max_int, self.complex)
                row += [ element, ]
            mat += [ row, ]
        return( Matrix(mat) )

    def check(self, answer):
        if np.shape(self.answer) == (1,1): self.answer = self.answer[0]
        if self.answer == answer: print('Correct')
        else:
            print('The correct answer is:')
            Print('$%s$'%myltx(self.answer))

#-------------------------------------------------------------------#
#                                                                   #
#               4.3 Change of Basis Problems                        #
#                                                                   #
#-------------------------------------------------------------------#

# Change of basis - 2D only!
def ip(v1, v2): return( (v1.transpose().conjugate() * v2)[0] )
class change_basis_problem:
    def __init__(self, max_int=5, complex=False):
        self.max_int = max_int
        self.complex=complex
        self.other_bases = []
        r2 = 1/sqrt(2)
        r5 = 1/sqrt(5)
        self.other_bases += [ (Matrix([(r2),(r2)]), Matrix([(r2),(-r2)])), ]
        self.other_bases += [ (Matrix([(2*r5),(r5)]), Matrix([(-r5),(2*r5)])), ]
        # I can have only this one complex basis without changing the code.
        self.other_bases += [(Matrix([(r2),(r2*i)]),Matrix([(r2),(-r2*i)])),]

    def show_bases(self):
        Print('$%s$'%myltx(self.other_bases))

    def new(self):
        # Generate a column
        v = Matrix([ 
            (rnum(1, self.max_int, self.complex)), 
            (rnum(1, self.max_int, self.complex)) ])
        # Pick a basis
        if self.complex: n = randint(0, len(self.other_bases)-1)
        else: n = randint(0, len(self.other_bases)-2)
        b = self.other_bases[n] 
        Print(r'Transform the vector $%s$ to the basis $\left(%s,%s\right)$'
        %(myltx(v), myltx(b[0]), myltx(b[1]) ))
        self.answer = col( ip(b[0],v), ip(b[1],v) )
        self.answer.simplify()

    def check(self, answer):
        if answer == self.answer: Print('Correct')
        else:
            Print('The correct answer is: $%s$'%myltx(self.answer))

#-------------------------------------------------------------------#
#                                                                   #
#               4.4 Electron Spin                                   #
#                                                                   #
#-------------------------------------------------------------------#

# The observation matrices
sigma_x = Matrix([(0,1),(1,0)])
sigma_y = Matrix([(0,-i),(i,0)])
sigma_z = Matrix([(1,0),(0,-1)])
# The x, y, an z (plus and minus) states
pZ = col(1,0); pZ.simplify
mZ = col(0,1); mZ.simplify
pX = col(1,1)/sqrt(2); pX.simplify
mX = col(1,-1)/sqrt(2); mX.simplify
pY = col(1,i)/sqrt(2); pY.simplify
mY = col(1,-i)/sqrt(2); mY.simplify

#
# 2 bit bases for X, Y, and Z
#

x_basis_2bit = []
x_basis_2bit += [ TP(pX, pX), ]
x_basis_2bit += [ TP(pX, mX), ]
x_basis_2bit += [ TP(mX, pX), ]
x_basis_2bit += [ TP(mX, mX), ]

y_basis_2bit = []
y_basis_2bit += [ TP(pY, pY), ]
y_basis_2bit += [ TP(pY, mY), ]
y_basis_2bit += [ TP(mY, pY), ]
y_basis_2bit += [ TP(mY, mY), ]

z_basis_2bit = []
z_basis_2bit += [ TP(pZ, pZ), ]
z_basis_2bit += [ TP(pZ, mZ), ]
z_basis_2bit += [ TP(mZ, pZ), ]
z_basis_2bit += [ TP(mZ, mZ), ]

def ket(state):
    if state == 1: return(pZ)
    if state == -1: return(mZ)
    if state == 2: return(pX)
    if state == -2: return(mX)
    if state == 3: return(pY)
    if state == -3: return(mY)
    print('Error: ket() got unknown state id: %s'%state)

def pket(state):
    #Print('pket() got: $%s$'%latex(state))
    state.simplify()
    if state == pZ: var = '+z'
    elif state == mZ: var = '-z'
    elif state == pX: var = '+x'
    elif state == mX: var = '-x'
    elif state == pY: var = '+y'
    elif state == mY: var = '-y'
    else: Print('Error in pket()')
    return(r'| %s \rangle' %var)

class electron_problem:
    def __init__(self, min_int=1, max_int=5, fractions=True, use_y=True, complex=True):
        self.max_int = max_int
        self.min_int = min_int
        self.complex=complex
        self.use_y=use_y
        # This contains a (more or less) latex version of the problem
        self.latex = ''
        self.debug = False

    def set_debug(self, tf):
        self.debug = tf
    def DBG(self, msg):
        if self.debug==True: Print(msg)

    def new(self):
        # Input and output bases
        if self.use_y: ibasis = randint(1, 3)
        else: ibasis = randint(1, 2)
        if self.use_y: obasis = randint(1, 3)
        else: obasis = randint(1, 2)

        # Probability amplitudes for the input state
        # Generate arbitrary alpha and beta
        alpha = rnum(self.min_int, self.max_int, self.complex)
        # If alpha is zero, then beta can't be zero
        if alpha==0:
            beta = rnum(self.min_int+1, self.max_int, self.complex)
        else:
            beta = rnum(self.min_int, self.max_int, self.complex)
        self.DBG(r'Original: $\alpha = %s,\;\;\;\beta=%s$'%(latex(alpha),latex(beta)))

        # If alpha and beta are anything but ints, simplify them
        if type(alpha) == type(1): pass
        else:alpha = alpha.simplify().expand()
        if type(beta) == type(1): pass
        else: beta = beta.simplify().expand()
        self.DBG(r'Simplified: $\alpha = %s,\;\;\;\beta=%s$'%(latex(alpha),latex(beta)))

        k1 = ket(+ibasis); k2 = ket(-ibasis)


        # Plus or minus (in front of the first term)
        if randint(0,1)==0: 
            fpm = '' # If the first term is plus, don't print anything
            fsign = 1
        else:
            fpm = '-'
            fsign = -1

        # Plus or minus (between the two terms of the state)
        if randint(0,1)==0: 
            bpm = '+'
            bsign = 1
        else:
            bpm = '-'
            bsign = -1

        # Create the actualy (numeric) input state
        istate = fsign*alpha*k1 + bsign*beta*k2
        self.DBG('istate = $'+latex(fsign)+'*'+latex(alpha)+'*'+latex(k1)+'+'+latex(bsign)+'*'+latex(beta)+'*'+latex(k2)+'$')
        self.DBG('istate = $%s$'%latex(istate))

        # Normalize alpha and beta. Reduce the fractions.
        # But, only reduce if the gcd is an Integer.
        length = istate.norm().simplify()
        #
        alpha_den = length
        gcd = sy.gcd(alpha, alpha_den)
        if type(gcd) == type(sy.Integer(2)):
            alpha /= gcd
            alpha_den /= gcd
        alpha_num_ltx = latex(alpha)
        alpha_dltx = latex(alpha_den)
        #
        beta_den = length
        gcd = sy.gcd(beta, beta_den)
        if type(gcd) == type(sy.Integer(2)):
            beta /= gcd
            beta_den /= gcd
        beta_num_ltx = latex(beta)
        beta_dltx = latex(beta_den)
        
        self.DBG('The length was $%s$' %latex(length))
        istate = istate/length # Normalize the state
        self.DBG('normalized istate = $%s$'%latex(istate))
        self.DBG('alpha_den=$%s$' %latex(alpha_den))
        self.DBG('beta_den=$%s$' %latex(beta_den))

        # Already normalized them in the gcd part???
        #alpha /= length; beta /= length # Normalize alpha and beta

        # Create latex for alpha
        self.DBG(r'Normalized: $\alpha = %s,\;\;\;\beta=%s$'%(latex(alpha),latex(beta)))
        if alpha == 0: alpha_ltx = '0'
        elif alpha_den == 1: alpha_ltx = fpm
        else:
            alpha_ltx = r'%s\dfrac{%s}{%s}'%(fpm, alpha_num_ltx, alpha_dltx)
        # and beta
        if beta == 0: beta_ltx = '0'
        elif beta_den == 1: beta_ltx = bpm 
        else:
            beta_ltx = r'%s\dfrac{%s}{%s}'%(bpm, beta_num_ltx, beta_dltx)

        msg = 'Given the input state: $\;'
        if alpha_ltx != '0': msg += alpha_ltx + pket(k1)
        if beta_ltx != '0': msg += beta_ltx + pket(k2)
        msg += '$'
        msg += '\n<br><br>\n'

        bases = ['z', 'x', 'y']
        omats = [ sigma_z, sigma_x, sigma_y ]

        # Choose one of the two two possible output eigenvalues
        if randint(0,1)==0:
            eval = +1
            sign = '+'
        else:
            eval = -1
            sign = '-'

        # Choose form of question
        if randint(0,1)==0:
            msg += 'If you do a measurment in the '
            msg += '$%s$ basis.'%bases[obasis-1]
            msg += '\n<br><br>\n'
            msg += r'What is the probability of getting $|%s%s\rangle$'\
                %(sign, bases[obasis-1])
        else:
            msg += 'If you do the observation associated with '
            msg += 'the matrix $%s$'%latex(omats[obasis-1])
            msg += '\n<br><br>\n'
            msg += 'What is the probability of seeing '
            msg += 'the $%s$ eigenvalue'%eval
            msg += '\n<br><br>\n'

        #
        # Here's where we actually print out the problems
        #
        Print(msg)
        self.latex = msg # Save it, also.

        # Get the actual (numeric) output state
        ostate = ket(eval*obasis)
        self.DBG('Actual output state is: $%s$'%latex(ostate))
        prob = inner_product(ostate,istate)*inner_product(istate,ostate)
        self.answer = prob.simplify()
        self.DBG('With probability: $\;%s\;$ or about $\;%0.4f$'
            %(latex(self.answer), self.answer.evalf()))

    def check(self, answer):
        if type(answer)==type(.2):
            if self.answer.evalf() == answer:
                Print('Correct')
                return
            if abs(self.answer.evalf() - answer) < .002:
                Print('Close enough!')
                return
        if self.answer == answer:
                Print('Correct')
                return
        msg = 'The answer is: $\;%s$, '%latex(self.answer.together())
        msg += 'or approximately $\;%0.4f$'%self.answer.evalf()
        Print(msg)

#####################################################################
#                                                                   #
#                       5. Electron Spin                            #
#                                                                   #
#####################################################################

# ------------------------------------------------------------------
# 3D Arrow - I found this code online. I have no idea how it works.
#
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)
#
# But calling it like this will draw the vector I want
# 
def draw_vec(ax, v, color='black', label=None, **args):
    V = Arrow3D([0,float(v[0])],[0,float(v[1])],[0,float(v[2])],
                mutation_scale=20, lw=1, arrowstyle="-|>", color=color, **args)
    ax.add_artist(V)
    fudge=1.1
    ax.text3D(v[0]*fudge, v[1]*fudge, v[2]*fudge,
        label, color=color, fontsize='16')

# ------------ Start: Electron Spin Problem -------------------------------

def electron_spin_problem(s_1, s_2, exact=False, ndigs=2, draw_box=False):
    from sympy import cos, sin
    from sympy import acos as arccos

    Print('<hr size=3') #################################

    # The spin vectors must have three components
    if np.shape(s_1) != (3, 1) or np.shape(s_2) != (3, 1):
        print('Error: s_1 and s_2 must be three element column vectors.')
        return

    # Need to fix this, but it will take some work
    fm = format(exact, ndigs).do_format

    # Make sure they're normalized
    if s_1.norm() != 1:
        Print('The vector $%s$ is not normalized.'%fm(s_1))
        Print('I\'m normalizing it for you ...')
        s_1 = s_1/s_1.norm()
    if s_2.norm() != 1:
        Print('The vector $%s$ is not normalized.'%fm(s_2))
        Print('I\'m normalizing it for you ...')
        s_2 = s_2/s_2.norm()

    # Juggle the components around to do the right thing
    # based on whether we're doing exact math or not, and
    # make sure they end up as sympy single column matrices
    # (which we use to represent column vectors).

    if exact:
        s_1 = [ s_1[0], s_1[1], s_1[2] ]
        s_2 = [ s_2[0], s_2[1], s_2[2] ]
    else:
        s_1 = [ float(s_1[0]), float(s_1[1]), float(s_1[2]) ]
        s_2 = [ float(s_2[0]), float(s_2[1]), float(s_2[2]) ]
    s_1 = Matrix(s_1); s_2 = Matrix(s_2)

    Print('The spin vectors are:')
    Print(r'$s_1 = %s,\;\; s_2=%s$'%( fm(s_1), fm(s_2) ))
    Print(r'Given an electron in the state $|+s_1\rangle$')
    txt = 'We want the probabilities of measuring '
    txt += r'$|+s_2\rangle$ or $|-s_2\rangle$'
    Print(txt)

    Print('<hr size=3') #################################

    Print('<b>Step 1:</b>')
    Print('Calculate the measurement operators for both spin vectors')
    txt = r'using the formula: $O = s_x\sigma_x + s_y\sigma_y + s_z\sigma_z$'
    Print(txt)

    O_1 = s_1[0]*sigma_x + s_1[1]*sigma_y + s_1[2]*sigma_z
    O_1.simplify()
    Print(r'Calculating the measurment operator for $s_1$ ...')
    Print('$O_1 = %s %s + %s %s + %s %s = %s$'
        %( fm(s_1[0]), fm(sigma_x), fm(s_1[1]),
        fm(sigma_y), fm(s_1[2]), fm(sigma_z), fm(O_1) ))

    O_2 = s_2[0]*sigma_x + s_2[1]*sigma_y + s_2[2]*sigma_z
    O_2.simplify()
    Print(r'Calculating the measurment operator for $s_2$ ...')
    Print('$O_2 = %s %s + %s %s + %s %s = %s$'
        %( fm(s_2[0]), fm(sigma_x), fm(s_2[1]),
        fm(sigma_y), fm(s_2[2]), fm(sigma_z), fm(O_2) ))

    Print('<hr size=3') #################################

    Print('<b>Step 2:</b>')
    Print('The computer will find the eigenvalues and eigenvectors')
    evals_1, evecs_1 = find_eigenvectors(O_1)
    Print('Eigenvalues for $O_1$ are: $%s$ and $%s$'
        %( fm(evals_1[0]), fm(evals_1[1]) ))
    Print('Eigenvectors for $O_1$ are: $%s$ and $%s$'
        %( fm(evecs_1[0]), fm(evecs_1[1]) ))

    evals_2, evecs_2 = find_eigenvectors(O_2)
    Print('Eigenvalues for $O_2$ are: $%s$ and $%s$'
        %( fm(evals_2[0]), fm(evals_2[1]) ))
    Print('Eigenvectors for $O_2$ are: $%s$ and $%s$'
        %( fm(evecs_2[0]), fm(evecs_2[1]) ))

    # The basis we'll be projecting onto is { +s_2, -s_2 }.
    # We'll get back the eivenvectors in order of the lowest
    # eigenvalue. So -s will come first.

    Print('The eigenvalues tell us which vectors represent the "plus"')
    Print('states and which represent the "minus" states')

    # The state vector | -s_2 >
    state_m2 = evecs_2[0]; state_m2.simplify()
    # The state vector | +s_2 >
    state_p2 = evecs_2[1]; state_p2.simplify()

    state_1 = evecs_1[1]; state_1.simplify()
    Print(r'$|+s_1\rangle=%s$'%fm(state_1))
    state_m2 = evecs_2[0]; state_m2.simplify()
    state_p2 = evecs_2[1]; state_p2.simplify()
    Print(r'$|+s_2\rangle=%s,\;\; |-s_2\rangle=%s$'
        %( fm(state_p2), fm(state_m2) ))

    Print('<hr size=3') #################################

    Print('<b>Step 3:</b>')
    Print('Now that we know what the $s_2$ state vectors are, we want to')
    txt = r'project our given state, $|+s_1\rangle$, onto to the $s_2$ basis. '
    txt += 'This will'
    Print(txt)
    txt = 'give us an expression such as:  '
    txt += r'$\alpha |+s_2\rangle + \beta |-s_2\rangle$'
    Print(txt)
    Print('from which we can calculate the probabilities.')

    # Calculation for alpha and beta
    ip = inner_product
    alpha = ip(state_p2, state_1)
    beta  = ip(state_m2, state_1)

    Print(r'<br>The projection onto $s_2$ goes like this:')
    txt = r'$\alpha = \langle +s_2 | +s_1 \rangle$'
    txt += r'$ = %s %s$'%( fm(state_p2.adjoint()), fm(state_1))
    txt += r'$= %s$'%fm(alpha)
    Print(txt)
    txt = r'$\beta = \langle -s_2 | +s_1 \rangle$'
    txt += r'$ = %s %s$'%( fm(state_m2.adjoint()), fm(state_1))
    txt += r'$= %s$'%fm(beta)
    Print(txt)

    Print(r'<br>The state written in the $s_2$ basis is:')
    txt = '<br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'
    txt += r'$(%s)|+s_2\rangle + (%s)|-s_2\rangle$' %( fm(alpha), fm(beta) )
    Print(txt)

    Print('<hr size=3') #################################
    Print('<b>Step 4:<b> Calculate the probabilities')
    
    # Probability calculation
    prob_p2 = alpha*alpha.conjugate()
    prob_m2 = beta*beta.conjugate()
    if exact==True:
        prob_p2 = prob_p2.simplify(); prob_m2 = prob_m2.simplify()
    else:
        prob_p2 = float(prob_p2); prob_m2 = float(prob_m2)
    
    txt = r'$P(|+s_2\rangle)=|\alpha|^2=\alpha\alpha^*$'
    txt += r'$= (%s)(%s)$'%( fm(alpha), fm(alpha.conjugate()) )
    txt += '$='+fm(prob_p2)+'$'
    Print(txt)

    txt = r'$P(|-s_2\rangle)=|\beta|^2=\beta\beta^*$'
    txt += r'$= (%s)(%s)$'%( fm(beta), fm(beta.conjugate()) )
    txt += '$='+fm(prob_m2)+'$'
    Print(txt)

    Print('<hr size=3') #################################
    Print('<b>Additional Note:</b>')

    # Note: We know that the norms of both vectors are one, so ... 
    cos_theta = ip(s_1, s_2)
    if exact==False: cos_theta = float(cos_theta)
    theta = arccos(cos_theta)
    Print('The angle between $s_1$ and $s_2$ is: $%s$'%fm(theta))

    Print(r'The probabilities should be equivalent to:')
    Print(r'$\cos^2(\theta)=%s$ and $\sin^2(\theta)=%s$'
      %( fm(cos(theta/2)**2), fm(sin(theta/2)**2) ))

    # Now draw a picture

    Print('<hr size=3') #################################
    Print('<b>And finally, here\'s a picture of the spin vectors:<b>')

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Draw the axes 
    draw_vec(ax, [1,0,0], color='black', label='x')
    draw_vec(ax, [0,1,0], color='black', label='y')
    draw_vec(ax, [0,0,1], color='black', label='z')

    # Draw the spin vectors
    draw_vec(ax, s_1, label='$s_1$', color='red')
    draw_vec(ax, s_2, label='$s_2$', color='green')
    
    if draw_box: pass
    else: ax.set_axis_off()

    # If you use "equal" then the last parm on position sets size?
    ax.set_position([1,1,1,2])
    ax.set_aspect("equal")
    ax.view_init(elev=5, azim=25)

# ------------ End: Electron Spin Problem -------------------------------

