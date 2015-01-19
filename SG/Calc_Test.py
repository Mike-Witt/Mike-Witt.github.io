
from pylab import *
import numpy as np
import sympy as sy
sy.init_printing()

# This function draws a line with the specified slope
# The line is centered at "point" and has a length of "length"
def line(slope, point, length):
    # Make sure I'm using floats
    x0,y0 = [float(point[0]),float(point[1])]
    slope = float(slope); length = float(length)
    x1 = x0-length/2; x2 = x0+length/2
    dx = abs(x2-x1)
    #print('dx=%s'%dx)
    dy = slope*dx
    #print('which gives dy=%s'%dy)
    y1 = slope*(x1-x0)+y0
    y2 = slope*(x2-x0)+y0
    L = sqrt((x2-x1)**2 + (y2-y1)**2)
    #print('Actual length = %s'%L)
    #print('scaling dx by wL/aL ...')
    dx = dx*length/L
    #print('gives dx=%s'%dx)
    x1 = x0-dx/2; x2 = x0+dx/2
    y1 = slope*(x1-x0)+y0
    y2 = slope*(x2-x0)+y0
    #print('plotting: %s,%s'%([x1, x2], [y1, y2]))
    plot([x1, x2], [y1, y2], color='black')
    #axes().set_aspect('equal')
    L = sqrt((x2-x1)**2 + (y2-y1)**2)
    #print('Final length = %s'%L)
    
def myplot(objs, scale=1, V=False):
    from pylab import rcParams
    
    min_x = inf; max_x = -inf; min_y = inf; max_y = -inf
    if V: print('objs=%s'%objs)
    for o in objs:
        x = o[0]; y = o[1]
        if V: print('x=%s, y=%s'%(x,y))
        minx = min(x); miny = min(y); maxx = max(x); maxy = max(y)
        if minx < min_x: min_x = minx
        if miny < min_y: min_y = miny
        if maxx > max_x: max_x = maxx
        if maxy > max_y: max_y = maxy
        if V: print('minx=%s, maxx=%s, miny=%s, maxy=%s'
              %(minx,maxx,miny,maxy))
        if len(o) > 2: args = o[2]
        else: args = {}
        plot(x, y, **args)
    if V: print('min_x=%s, max_x=%s, min_y=%s, max_y=%s'
          %(min_x,max_x,min_y,max_y))
    xlim(min_x, max_x)
    delta_y = max_y - min_y
    y_fudge = delta_y * .2
    ylim(min_y - y_fudge, max_y + y_fudge) 

def show_derivative(f, the_range, where, figsize=[6,4]):
    import sympy as sy
    import numpy as np
    
    # This has to be done before the first plot!
    # Notes: Setting the figsize here *does* fix the problem
    #    of having to plot twice to get it right. However,
    #    setting the aspect here does *not* fix the "blurry" problem,
    #    and it actually screws things up!
    #
    width = float(figsize[0]); height = float(figsize[1])
    figure(figsize=[width,height])
    # Do *not* do the next two lines
    #ax = axes()
    #ax.set_aspect(width/height)
    
    x = sy.Symbol('x')
    df = sy.diff(f, x).simplify()
    print('f = %s, df = %s' %(f, df))
    
    xlow = float(the_range[0])
    xhigh = float(the_range[1])
    px = float(where)
    py = float(f.subs(x,px))
    X = np.linspace(xlow, xhigh, 100)
    
    # This may not play completely with myplot
    slope = float(df.subs(x,px))
    line(slope, [px,py], 1)
    
    # Mark the point of interest
    #plot(px, py, marker='.', color='black', markersize=6)
    
    Y1 = []
    for x0 in X: Y1 += [ float(f.subs(x,x0)), ]
    Y2 = []
    for x0 in X: Y2 += [ float(df.subs(x,x0)), ]
    
    # The function
    p1 = (X, Y1, {'label':r'$f=%s$'%sy.latex(f),
        'color':'blue', 'linestyle':'-', 'alpha':.5})
    # The derivative
    p2 = (X, Y2, {'label':r'$\frac{df}{dx}=%s$'%sy.latex(df),
        'color':'red', 'alpha':.5})
    # The x axis
    p3 = ([xlow,xhigh], [0,0], {'color':'gray'})
    objects = [p1, p2, p3]
    myplot(objects)
    
    # Vertical line at the point of interest
    ylim = axes().get_ylim()
    xlim = axes().get_xlim()
    plot([px,px], [ylim[0], ylim[1]], color='gray', linestyle='--')
    #title('Function and Derivative', fontsize=12)
    legend(bbox_to_anchor=(1.55,1), prop={'size':20})
    show() # Need this when not in the notebook
    
x = sy.Symbol('x')

f = sy.sin(x)

#f = x*(x-1)*(x+2)
#f = f.expand()

show_derivative(f, [-pi, pi], pi/2, figsize=[6,4])

"""
ax = axes()
lowx,highx = ax.get_xbound()
lowy,highy = ax.get_ybound()
fig = ax.get_figure()
print 'bounds:'
print highx - lowx, highy - lowy
print 'aspect:'
print ax.get_aspect()
print 'figsize:'
print fig.get_figwidth(), fig.get_figheight()
"""

