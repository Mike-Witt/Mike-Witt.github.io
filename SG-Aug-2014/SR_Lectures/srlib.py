#
#   srlib.py (MODIFIED VERSION!)
#
# This is a "cut-down" version that was meant for the Summerm 2014 SG
# group, while we were studying relativity. It doesn't have a number of
# things that are in the "real" srlib. Like the "twin" code, for example.
#

import numpy as np
import matplotlib
import matplotlib.pyplot as pyplot

c = 299792458. # See: http://en.wikipedia.org/wiki/Speed_of_light
year = 60*60*24*365 # Approximate seconds in a typical year
lightyear = c*year
g = 9.80    # Note: Different problems may use slightly different values

def gamma(v):
    if v >= c: raise(Exception('gamma: v=%s >= c'%v))
    return( 1/np.sqrt(1-(v/c)**2 ) )

# Lorentz transforms for x and t
def ltx(x,t,v): return( gamma(v)*(x-v*t) )
def ltt(x,t,v): return( gamma(v)*(t-v*x/c**2) )

#Acceleration formula from Beiser page 26
#def trans_accel(accel, velocity):
#    a_prime = accel*(1-velocity**2/c**2)**(3./2.)
#    return( a_prime)
def trans_accel(a, v):
    return( a/(gamma(v)**3) )

    #############################################################
    #                                                           #
    #           Code for space-time diagrams                    #
    #                                                           #
    #############################################################
    
class world_line:
    # Currently a world line is just a named set of x and t values.
    def __init__(self, name, x, t, color='black'):
        self.name = name
        self.x = x # The position values, in metres
        self.t = t # The time values, in seconds
        self.color = color
        
    # A world line knows how to perform the Lorentz tranform on itself,
    # given a veloctiy v, returning a new transformed world line.
    def ltrans(self, v):
        x_prime = []; t_prime = []
        for i in range(len(self.x)):
            x = self.x[i]
            t = self.t[i]
            x_prime += [ ltx(x,t,v), ]
            t_prime += [ ltt(x,t,v), ]
        return( world_line(self.name, x_prime, t_prime, self.color) )

class event:
    def __init__(self, name, x, t, color='black'):
        self.name = name
        self.x = x # The position values, in metres
        self.t = t # The time values, in seconds
        self.color = color
        
    # An event knows how to perform the Lorentz tranform on itself,
    # given a veloctiy v, returning a new transformed event.
    def ltrans(self, v):
        x = self.x
        t = self.t
        x_prime = ltx(x,t,v)
        t_prime = ltt(x,t,v)
        return( event(self.name, x_prime, t_prime, self.color) )

# A grid of time and space "iso" lines. These are made out of world
# lines, so that I can use the world lines ltrans to transform them.
#
# Notes:
#
#   - bounds: I really do need these. They are the same bounds that
#     would be use for the figure.
#
#   - x_spacing, t_spacing, nx, nt: If these are not supplied I'll default
#     to something reasonable. If they are all supplied, the "spacing"
#     parameters will take precidence.
#
#   - I don't (so far) have a strategy for "centering" the grid at
#     a point of interest. I can't assume that the "origin" will
#     always be in the picture.
#
#   --> An alternative approach would be to just specifiy a list of
#       velocities to the frame display method, and have display
#       figure out what to do.
#
class grid:
    def __init__(self, 
        bounds, x_spacing=None, t_spacing=None, nx=None, nt=None,
        color='black', thickness=.2, V=False):

        self.bounds = bounds
        self.x_spacing = x_spacing
        self.t_spacing = t_spacing
        self.nx = nx
        self.nt = nt
        self.color = color
        self.thickness = thickness

        left = bounds[0]; right = bounds[1]
        bottom = bounds[2]; top = bounds[3]

        if x_spacing == None:
            if nx == None: nx = 10
            dx = (right - left) / nx
        else:
            dx = x_spacing
            nx = int((right - left) / dx)

        if t_spacing == None:
            if nt == None: nt = 10
            dt = (top - bottom) / nt
        else:
            dt = t_spacing
            nt = int((top - bottom) / dt)

        if V: print('grid(): nx=%d, dx=%s, nt=%s, dt=%s'%(nx, dx, nt, dt))

        self.tlines = [] # Time grid lines
        self.xlines = [] # Space grid lines

        for n in range(nt):
            # Start at the bottom and draw the time lines going up
            X = [ left, right]
            T = [ bottom + n*dt, bottom + n*dt ]
            self.tlines += [ world_line('', X, T, color), ]

        for n in range(nx):
            # Start at the left and draw space lines going right
            X = [ left + n*dx, left + n*dx ]
            T = [ bottom, top ]
            self.xlines += [ world_line('', X, T, color), ]

    def ltrans(self, v):
        newgrid = grid(
            self.bounds, self.x_spacing, self.t_spacing, self.nx, self.nt,
                self.color, self.thickness)
        oldlines = newgrid.tlines; newlines = []
        for line in oldlines: newlines += [ line.ltrans(v), ]
        newgrid.tlines = newlines
        oldlines = newgrid.xlines; newlines = []
        for line in oldlines: newlines += [ line.ltrans(v), ]
        newgrid.xlines = newlines
        return(newgrid)
    
class frame:
    # A reference frame has a name and it holds a bunch of world lines.
    #
    # When you create your first frame, first you have to create all
    # the world lines you want and put them in the array "wlines" then do:
    #
    #    a_frame = frame("a name", wlines)
    #
    # If you're creating a frame that is traveling with velocity V wrt
    # some existing frame you do it like this:
    #
    #    new_frame = frame("new_name", old_frame.wlines, v=V)
    #
    def __init__(self, name, wlines=[], events=[], grids=[], v=0):
        self.name = name
        self.wlines = []
        self.events = []
        self.grids = []
        for wl in wlines:
            self.wlines += [ wl.ltrans(v), ]
        for event in events:
            self.events += [ event.ltrans(v), ]
        for grid in grids:
            self.grids += [ grid.ltrans(v), ]
            
    # Make a grid of time and space "iso" lines. There may be a more
    # efficient ways to get this to work, but right now I'll just make
    # it out of world lines.
            
    # A frame can display itself
    def display(self, figsize=[8,8], points=False, scale='si',
        bounds=None, nticks=5, V=False):

        from matplotlib.ticker import MaxNLocator
        global space_ticks, time_ticks, space_labels, time_labels

        for wl in self.wlines:
            if points: var = pyplot.plot(
                wl.x, wl.t, marker='o', linestyle='_', color=wl.color)
            else: var = pyplot.plot(
                wl.x, wl.t, label=wl.name, color=wl.color, linewidth=2)
        for event in self.events:
            var = pyplot.plot([event.x], [event.t], 
                label=event.name, marker='o', linestyle='_', color=event.color)

        for grid in self.grids: 
            for gridline in grid.tlines:
                times = np.array(gridline.t)
                var = pyplot.plot(gridline.x, times, color=grid.color,
                    linestyle='-', linewidth=grid.thickness)
            for gridline in grid.xlines:
                times = np.array(gridline.t)
                var = pyplot.plot(gridline.x, times, color=grid.color,
                    linestyle='-', linewidth=grid.thickness)
                
        matplotlib.rcParams['figure.figsize'] = figsize
        pyplot.tick_params(labelbottom='on')
        pyplot.tick_params(labelleft='on')
        pyplot.title(self.name)
        pyplot.legend()

        ax = pyplot.axes()
        ax.set_aspect(c)
        if bounds != None:
            pyplot.axis( [bounds[0], bounds[1], bounds[2], bounds[3]] )

        if V:
            print('After doing all the plots ...')
            ax = pyplot.axes()
            print('  ax.get_aspect() = %s' %ax.get_aspect())
            print('  ax.get_xbound() = %s' %[ax.get_xbound()])
            print('  ax.get_ybound() = %s' %[ax.get_ybound()])
            print('  ax.get_ybound()[1]*c = %s' %[ax.get_ybound()[1]*c])

        if bounds == None:
            bounds = []
            bounds += [ax.get_xbound()[0],]
            bounds += [ax.get_xbound()[1],]
            bounds += [ax.get_ybound()[0],]
            bounds += [ax.get_ybound()[1],]

        if V: print('bounds: %s'%bounds)

        space_ticks, junk_labels = pyplot.xticks()
        time_ticks = np.linspace(bounds[2], bounds[3], nticks)
        # We default to using the original ticks as labels.
        space_labels = space_ticks
        time_labels = time_ticks

        # Adjust the scale and axis labels based on the "scale" option

        if scale == 'si' or scale == 'SI':
            pyplot.xlabel('Metres')
            pyplot.ylabel('Seconds')

        elif scale == 'metres' or scale == 'meters':
            pyplot.xlabel('Metres')
            pyplot.ylabel('Metres')

            # For whatever reason, it didn't work to combine division by c
            # with conversion to a string. So I do the scaling of the list
            # all at once, in the following line. Then I do the formating
            # in the loop below
            time_labels = time_labels*c # Scale here
            for n in range(len(time_labels)):
                time_labels[n] = '%.1f'%(time_labels[n]) # format here

        elif scale == 'seconds' or scale == 'light-seconds':
            pyplot.xlabel('Light Seconds')
            pyplot.ylabel('Seconds')
            space_labels = space_labels/c
            for n in range(len(space_labels)):
                space_labels[n] = '%.1f'%(space_labels[n])

        elif scale == 'years' or scale == 'light-years':
            pyplot.xlabel('Light Years')
            pyplot.ylabel('Years')
            space_labels = space_labels/c
            space_labels = space_labels/year
            time_labels = time_labels/year
            for n in range(len(space_labels)):
                space_labels[n] = '%.1f'%(space_labels[n])
            for n in range(len(time_labels)):
                time_labels[n] = '%.1f'%(time_labels[n]) # format here

        else:
            raise(Exception('I don\'t know the "%s" scale options'%scale))
        ax.set_xticks(space_ticks)
        ax.set_xticklabels(space_labels)

        # The next two lines just deletes the first time label. This is to
        # prevent printing two separate "0"s at the origin.
        time_ticks = time_ticks[1:]
        time_labels = time_labels[1:]

        ax.set_yticks(time_ticks)
        ax.set_yticklabels(time_labels)

        # Return the pyplot object so that the user can do final adjustments
        return(pyplot)
