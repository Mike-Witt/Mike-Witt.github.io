#
#   srlib.py   
#
# This code is meant to be included in IPython Notebooks that want to
# create space-time diagrams, etc.
#
# Status:
#
#   May 2014 - Basic ST Diagrams are working
#   June 2014 - "Twins" code in progress
#   See the notebook(s) for further details
#
# Notes:
#
# - In this library we are assuming SI units. Mainly so that the
#   formulae look like they would in a typical textbook. There is
#   a scaling option on the method that actually displays a "frame"
#   so that you can easily display axes, consisting of either metres
#   or seconds / light-seconds.
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

    #############################################################
    #                                                           #
    #               Twins Class and related code                # 
    #                                                           #
    #############################################################

import os
from IPython.core.display import HTML
import sys
sys.path.append('/home/mike/Projects/Computer_Models/PDEs')
from model import run_stats

# This is formula (2.55) on page 36 (latex viewer page 50)
# of the new GR book.
#
# Notes:
# (1) I'm assuming that x0 = 0 in my version of the formula.
#     If I want to use the formula for the entire twin voyage,
#     I'll need to not only include x0, but I'll have to develop
#     the formula further to include v0. I *think* the book
#     gives enough information that I can figure this out.
#
# (2) 'g' is the constant acceleration. It must be constant,
#     but can be any value.
#
# (3) It calculates the postion x, based on the current time
#     (not tau). In other words, it's NOT dx.
#
def formula_255(g, t):
    x =  (c**2/g) * np.sqrt(1 + (g**2/c**2) * t**2) - c**2/g
    return(x)
    
class twins:

    def __init__(self):
        s = self
        s.req_travel_time = None
        s.req_destination = None
        s.V = False
        s.TRACE = False
        s.user = tuc_null(s)
        s.TEMP_DIR='/tmp/mike_temp_animation_frames/' # Needs ending '/'
        s.DEST_FILE='/var/www/Animations/twins.gif'  # Full path 
        s.bounds = None
        s.req_travel_time = None
        s.req_destination = None

    # run() Method
    #
    # Variables:
    #
    #   s = self # Just for convenience
    #
    #   --- These are (potentially) set by the user ---
    #
    #   s.a_t  # Constant acceleration in traveler's frame
    #          # Must be set by caller.
    #
    #   s.dt_s # Time increment (seconds).
    #          # Must be set by caller.
    #
    #   s.req_destination # Meters to requested destination.
    #                     # Must be set if requesting specific distance.
    #                     # Can always be set. Results in on line on plot.
    #
    #   s.req_travel_time # Requested travel in seconds.
    #                     # If this is set we are in "travel time" mode.
    #                     # If it is None, then we're in "distance" mode.
    #
    #   s.V # Verbose flag
    #
    #   s.user # User supplied function which provides callbacks, for
    #          # example at every time increment (dt). This is set to
    #          # null in the twin init method. Optionally supplied by user.
    #
    #   --- These are calculated in the program ---
    #
    #   # Set the all variables in both frames of reference to zero.
    #   s.t_s  # Time in stationary frame. Typically called "t."
    #   s.t_t  # Time in traveling frame. Typically called "tau."
    #   s.x_s  # Traveler's position in stationary frame.
    #   s.x_t  # Traveler's position in traveling frame.
    #   s.v_s  # Traveler's velocity in stationary frame.
    #   s.a_s  # Acceleration in stationary frame. Initialzed to a_t.
    #   s.dt_t # dtau / starts out equal to dt.
    #   s.nt = # Number of time points actually calculated.
    #   s.total_x_s # Total distance traveled, stationary frame.

    def run(s):
        # Make sure things are floats
        if s.req_destination != None:
            s.req_destination = float(s.req_destination) 
        if s.req_travel_time != None:
            s.req_travel_time = float(s.req_travel_time)
        s.a_t = float(s.a_t)
        s.dt_s = float(s.dt_s)
        s.a_s = s.a_t   # Accel starts out the same as in traveling frame.
        s.dt_t = s.dt_s # dtau starts out equal to dt
        s.t_s = 0.
        s.t_t = 0.
        s.x_s = 0.
        s.x_t = 0.
        s.v_s = 0.
        s.nt = 0 # This can be an integer
        s.total_x_s = 0.

        # If the Verbose flag is set, these DeBuG messages will be printed
        def DBG(msg):
            if s.V: print(msg)
            if s.TRACE: Trace(msg) # Copy to trace file, if tracing

        # If the TRACE flag is set, then we'll output a potentially
        # large amount of information the a file.
        if s.TRACE: s.tracefile = open('tw_tracefile.txt', 'w')
        def Trace(msg):
            if s.TRACE:
                s.tracefile.write(msg+'\n')
                s.tracefile.flush()

        # --- These are just tests to see "how far along" we are --- #

        def half_way_out():
            if s.HALF_WAY_OUT: return(True)
            if s.req_travel_time == None:
                # We're doing requested distance. If we've gone half that
                # distance, in the stationary frame, then we're half way out.
                if s.x_s >= .5 * s.req_destination:
                    s.HALF_WAY_OUT = True
                    DBG('Half way out, t_s = %0.2f years'%(s.t_s/year))
                    return(True)
            else:
                # We're doing requested travel time. If we've gone 1/4 of
                # that, in the stationary frame, then we're half way out.
                if s.t_s >= .25 * s.req_travel_time:
                    s.HALF_WAY_OUT = True
                    DBG('Half way out, t_s = %0.2f years'%(s.t_s/year))
                    return(True)
            return(False)

        def all_the_way_out():
            if s.ALL_THE_WAY_OUT: return(True)
            if s.req_travel_time == None:
                # We're doing requested distance. If we've done that,
                # in the stationary frame, then we're all the way out.
                if s.x_s >= s.req_destination:
                    s.ALL_THE_WAY_OUT = True
                    DBG('All the way out, t_s = %0.2f years'%(s.t_s/year))
                    return(True)
            else:
                # We're doing requested travel time. If we've gone 1/2 of
                # that, in the stationary frame, then we're half way out.
                if s.t_s >= .5 * s.req_travel_time:
                    s.ALL_THE_WAY_OUT = True
                    DBG('All the way out, t_s = %0.2f years'%(s.t_s/year))
                    return(True)
            return(False)

        def half_way_back():
            if s.HALF_WAY_BACK: return(True)
            if s.req_travel_time == None:
                # We're doing requested distance. If we're <= 1/2 that,
                # in the stationary frame, then we're half way back.
                if s.ALL_THE_WAY_OUT and s.x_s <= .5 * s.req_destination:
                    s.HALF_WAY_BACK = True
                    DBG('Half way back, t_s = %0.2f years'%(s.t_s/year))
                    return(True)
            else:
                # We're doing requested travel time. If we've gone 3/4 of
                # that, in the stationary frame, then we're half way back.
                if s.t_s >= .75 * s.req_travel_time:
                    s.HALF_WAY_BACK = True
                    DBG('Half way back, t_s = %0.2f years'%(s.t_s/year))
                    return(True)
            return(False)

        def all_the_way_back():
            if s.ALL_THE_WAY_BACK: return(True)
            if s.req_travel_time == None:
                # We're doing requested distance. If we're <= 0,
                # in the stationary frame, then we're back home.
                if s.HALF_WAY_BACK and s.x_s <= 0:
                    s.ALL_THE_WAY_BACK = True
                    DBG('All the way back, t_s = %0.2f years'%(s.t_s/year))
                    if s.ALL_THE_WAY_BACK: DBG('ALL_THE_WAY_BACK set OK')
                    else: DBG('ALL_THE_WAY_BACK not set!')
                    return(True)
            else:
                # We're doing requested travel time. If we've done it all,
                # in the stationary frame, then we're back home.
                if s.t_s >= s.req_travel_time:
                    s.ALL_THE_WAY_BACK = True
                    DBG('All the way back, t_s = %0.2f years'%(s.t_s/year))
                    if s.ALL_THE_WAY_BACK: DBG('ALL_THE_WAY_BACK set OK')
                    else: DBG('ALL_THE_WAY_BACK not set!')
                    return(True)
            return(False)

        # ------------ End of the "how far along" tests --------------- #

        s.HALF_WAY_OUT = False
        s.ALL_THE_WAY_OUT = False
        s.HALF_WAY_BACK = False
        s.ALL_THE_WAY_BACK = False

        # The user "time interval" function.
        s.user.ti(s) # Call it here for t = 0

        # %%% may want to change this so there is one routine to set
        # the status and then just check the constants in the code below %%%
        while not s.ALL_THE_WAY_BACK:

            if not half_way_out():
                # We are on the outbound journey, accelerating
                s.a_s = trans_accel(s.a_t, s.v_s) 
                
            elif not all_the_way_out():
                # We are on the outbound journey, decelerating
                s.a_s = trans_accel(-s.a_t, s.v_s) 

            elif not half_way_back():
                # We are on the inbound journey. We continue to "decelerate"
                # (meaning that we accelerate in the other direction).
                s.a_s = trans_accel(-s.a_t, s.v_s) 

            elif not all_the_way_back():
                # We are more than half-way home and braking
                s.a_s = trans_accel(s.a_t, s.v_s)

            # Do all the calculations for this time increment
            s.v_s += s.a_s*s.dt_s
            s.dt_t = s.dt_s/gamma(s.v_s)
            s.t_s += s.dt_s
            s.t_t += s.dt_t
            s.dx_s = s.v_s*s.dt_s
            s.x_s += s.dx_s

            s.total_x_s += abs(s.dx_s)
            s.nt += 1
            s.user.ti(s)
            Trace(\
                'nt=%s, t_s=%s, t_t=%s, x_s=%s, v_s=%s'\
                %(s.nt, s.t_s, s.t_t, s.x_s, s.v_s)\
                )

        # After we exit the while loop, do any "final" user code
        s.user.final(s)
        if s.TRACE: s.tracefile.close()

    #
    # Plot() method. Draw a static picture of the trip.
    #
    # This can be called after run() assuming that the run has used
    # the tuc_plot() user function.
    #
    def plot(s, bounds=None, figsize=[8,8]):

        if bounds == None:
            t0 = s.Ts[0]; t1 = s.Ts[len(s.Ts)-1]
            x0 = s.Es[0]; x1 = s.Ds[0]
            bnds = [ x0-.1*x1, x1+.1*x1, t0-.1*t1, t1+.01*t1 ] 
        else:
            bnds = bounds
        gbnds = [bnds[0]*4, bnds[1]*4, bnds[2]*4, bnds[3]*4]

        # If we have an "event" then create a grid at that point
        if s.event_position != None:
            """
            # This is good for the one light year example
            the_grid = grid(bounds=gbnds,
                            x_spacing=.1*lightyear,
                            t_spacing=.3*year,
                            color='green').ltrans(-s.event_velocity)
            """
            # But in general it works better to specify nx and nt
            the_grid = grid(bounds=gbnds,
                            nx=60, nt=60,
                            color='green').ltrans(-s.event_velocity)
            the_event = event('', s.event_position, s.event_time, 'black')

        # Make world lines - all in Earth's frame
        earth_wl = world_line('Earth', s.Es, s.Ts, 'black')
        dest_wl  = world_line('Dest.', s.Ds, s.Ts, 'brown')
        light_wl = world_line('light', s.Ls, s.Ts, 'yellow')
        twin_wl  = world_line('Twin',  s.Xs, s.Ts, 'green')
        earth_wls = []
        earth_wls += [ earth_wl, ]
        earth_wls += [ dest_wl, ]
        earth_wls += [ light_wl, ]
        earth_wls += [ twin_wl, ]

        # Create and display Earth frame
        earth_frame = frame('Earth Frame', earth_wls)
        if s.event_position != None: 
            earth_frame.grids = [ the_grid, ]
            earth_frame.events = [ the_event, ]

        pyplot = earth_frame.display(
            figsize=figsize, bounds=bnds, scale='light-years')

        foo=pyplot.axis('auto') # This one needs unequal scale
        foo=pyplot.xlim(bnds[0], bnds[1])
        foo=pyplot.ylim(bnds[2], bnds[3])
        pyplot.title("The Twins")

    def make_frame(s):

        if s.bounds == None:
            t0 = s.Ts[0]; t1 = s.Ts[len(s.Ts)-1]
            x0 = s.Es[0]; x1 = s.Ds[0]
            bnds = [ x0-.1*x1, x1+.1*x1, t0-.1*t1, t1+.01*t1 ]
        else:
            bnds = s.bounds
        gbnds = [bnds[0]*4, bnds[1]*4, bnds[2]*4, bnds[3]*4]

        """" 
        # This is good for the one light year example
        the_grid = grid(bounds=gbnds,
                        x_spacing=.1*lightyear,
                        t_spacing=.3*year,
                        color='green').ltrans(-s.v_s)
        """
        # But in general it works better to specify nx and nt
        the_grid = grid(bounds=gbnds,
                        nx=60, nt=60,
                        color='green').ltrans(-s.v_s)
        the_event = event('', s.x_s, s.t_s, 'black')

        # Make world lines - all in Earth's frame
        earth_wl = world_line('Earth', s.Es, s.Ts, 'black')
        dest_wl  = world_line('Dest.', s.Ds, s.Ts, 'brown')
        light_wl = world_line('light', s.Ls, s.Ts, 'yellow')
        twin_wl  = world_line('Twin',  s.Xs, s.Ts, 'green')
        earth_wls = []
        earth_wls += [ earth_wl, ]
        earth_wls += [ dest_wl, ]
        earth_wls += [ light_wl, ]
        earth_wls += [ twin_wl, ]

        # Create and display Earth frame
        earth_frame = frame('Earth Frame', earth_wls)
        earth_frame.grids = [ the_grid, ]
        earth_frame.events = [ the_event, ]

        if s.first_plot_done:
            s.pyplot.cla() # Clear the previous plot

        # %%% Ought to have a way to pass figsize!!! %%%
        s.pyplot = earth_frame.display(
            figsize=[8,8], bounds=bnds, scale='light-years')
        s.first_plot_done = True
        foo=s.pyplot.axis('auto') # This one needs unequal scale
        foo=s.pyplot.xlim(bnds[0], bnds[1])
        foo=s.pyplot.ylim(bnds[2], bnds[3])
        s.pyplot.title("The Twins")

        s.pyplot.savefig(s.TEMP_DIR + 'frame%05d'%s.frame_num + '.png')
        s.pyplot.close() # Prevent the plot from displaying in notebook
        s.frame_num += 1

    # 
    # ANIMATE THE FRAMES ALREADY IN THE SPECIFIED DIRECTORY
    #
    # This takes all the png files in the directory 'anim_dir' and creates
    # the animated gif file 'gif_name.gif' in that same directory.
    #
    # This should really be a subroutine in the IPython Project
    # directory. Not here.

    def animate_directory(self, anim_dir, name, gif_delay=20, gif_loop=1):
        import sys
        import os

        #These are the only options for now
        doGif = True 
        gif_name = 'animation'

        if doGif:
            if gif_delay==None:
                cmnd = 'cd %s; convert frame*.png '%anim_dir\
                    + '%s.gif'%gif_name
            else:
                cmnd = 'cd %s; convert '%anim_dir\
                    + '-delay %s -loop %s frame*.png '%(gif_delay, gif_loop)\
                    + '%s.gif'%gif_name
            print('Convert command:')
            print(cmnd)
            print('Working ...')
            sys.stdout.flush()
            os.system(cmnd)
            print('Convert command completed')

            # Leave the animation file in the specified place
            # "name" if the full pathname with filename and extension
            cmnd = 'mv ' + anim_dir + gif_name + '.gif ' + name
            print(cmnd)
            os.system(cmnd)

    ################ END OF ANIMATION METHOD #####################

def stats_message(s):
    if s.req_travel_time != None:
        msg = '\n    nt: %1.2e, t_s: %0.2f, x_s: %0.2f, travel time: %0.2f'\
        %(s.nt, s.t_s/year, s.x_s/lightyear, s.req_travel_time/year)
    else:
        msg = '\n    nt: %1.2e, t_s: %0.2f, x_s: %0.2f, destination: %0.2f'\
        %(s.nt, s.t_s/year, s.x_s/lightyear, s.req_destination/lightyear)
    return(msg)

#
# Twins "user code" classes
#

class tuc_null:
    def __init__(self, s): return
    def ti(self, s): return
    def final(self, s): return

class tuc_timer:
    def __init__(self, s): s.stats = run_stats()
    def ti(self, s): return
    def final(self, s): s.stats.final('Done')

class tuc_stats:
    def __init__(self, s): s.stats = run_stats()
        
    def ti(self, s):
        rc = s.stats.update(stats_message(s))

    def final(self, s): s.stats.final('Done')
    
class tuc_one:
    # Note that we don't save anything in "self" for this class. The twins
    # object will hand us its "self" called "s" here. We store everything
    # in the twins object. That's where the caller will look afterward.

    def __init__(self, s):
        s.X = []; s.T = []; s.E = []; s.D = [];
        s.Tau = []; s.DT = []; s.DTau = []
        s.AS = []; s.AT = []; s.V = [];
        s.maxv = 0
        s.found_999 = False

        # Since I don't know how many time points there will be, I can
        # only specify what fraction I want to save.
        s.ti_mod_factor = 1000 # Save 1 out of this many points

        s.stats = run_stats()

    # This function is called at every dt interval, incuding t=0
    def ti(self, s):

        rc = s.stats.update('')
    
        # Save max velocity
        if s.v_s > s.maxv: s.maxv = s.v_s

        # Look for .999c
        if (not s.found_999) and s.v_s >= .999*c:
            print('--> At %s of light speed'%(s.v_s/c))
            print('--> t = %0.2f, tau = %0.2f' %(s.t_s/year, s.t_t/year))
            print('--> pos = %0.2f light years' %(s.x_s/lightyear))
            s.found_999 = True
    
        # Only save the specified number of points
        if np.mod(s.nt, s.ti_mod_factor) != 0: return
        
        s.X += [ s.x_s, ]
        s.E += [ 0, ]
        s.D += [ s.req_destination, ]
        s.V += [ s.v_s/c, ]
        s.T += [ s.t_s, ]
        s.Tau += [ s.t_t, ]
        s.DT += [ s.dt_s, ]
        s.DTau += [ s.dt_t, ]
        s.AS += [ abs(s.a_s), ]
        s.AT += [ s.a_t, ]

    def final(self, s):
        s.stats.final('Done')

class tuc_plot:

    def __init__(self, s):

        # The following variables are set by the user, to control
        # aspects of the plot or animation. This means that they
        # must be set *after* the tuc_plot() class is instantiated.
        #
        s.grid_at_nt = None # For plot() - grid at specific time slice
        s.grid_at_ts = None # For plot() - grid at specific time
        s.ti_mod_factor = 1000
        s.animate = False
        s.animate_mod_factor = 1000 # (Animate only)

        # ------- End of user set variables ------------#

        s.stats = run_stats()
        s.event_position = None # We might create one in the ti() routine
        s.grid_done = False # For plot(). Has a grid been made yet.

        # The lower case "s" in each name is a reminder that each of
        # these variables is as seen from the stationary frame.
        s.Ts = [] # All the calculated times in the stationary frame
        s.Xs = [] # All the calculated positions in the stationary frame
        s.Ds = [] # The position of the destination as given by the user
        s.Es = [] # The position of the starting point (Earth), likely zero
        s.Ls = [] # The position of a light "ray" at east time
        s.NT = [] # Number of the current time slice
        #
        s.XX = [] # Position by the "2.55" formula

        # These variables are for a potential animation
        s.frame_num = 0
        s.bounds = None
        s.first_plot_done = False

        # Delete and re-create the temp directory.
        # Again, this is "just in case" we're animating.
        # (May want to make a separate tuc_animation class?)
        cmnd = 'rm -rf ' + s.TEMP_DIR
        print(cmnd)
        os.system(cmnd)
        cmnd = 'mkdir ' + s.TEMP_DIR
        print(cmnd)
        os.system(cmnd)

    def ti(self, s):
        rc = s.stats.update(stats_message(s))

        # We will save one out of ti_mod_factor time frames + the last frame
        if  np.mod(s.nt, s.ti_mod_factor) == 0 or s.ALL_THE_WAY_BACK:
            s.Ts += [ s.t_s, ]
            s.Xs += [ s.x_s, ]

            # I've got some work to do before I can calculate the whole
            # path using 255.
            if not s.HALF_WAY_OUT:
                xx = formula_255(s.a_t, s.t_s)
                #print('t = %1.2e, xs = %s, xx = %s'%(s.t_s, s.x_s, xx))
                s.XX += [ xx, ]
            else:
                s.XX += [ 0, ]

            if s.req_destination != None: s.Ds += [ s.req_destination, ]
            s.Es += [ 0, ]
            s.Ls += [ s.t_s*c, ]
            # Saving the time slice number so that we can identify exactly
            # what time is was at a given time slice, etc.
            s.NT += [ s.nt, ]

        if not s.animate and not s.grid_done and (\
            (s.grid_at_nt != None and s.grid_at_nt == s.nt)\
            or (s.grid_at_ts != None and s.t_s >= s.grid_at_ts)\
            ):
            print('Making grid: nt=%s, x_s=%0.2f, v_s=%0.3f, t_s=%0.2f'
                %(s.nt, s.x_s/lightyear, s.v_s/c, s.t_s/year))
            s.event_position = s.x_s
            s.event_velocity = s.v_s
            s.event_time = s.t_s
            s.grid_done = True

        # NOTE: animate_mod_factor needs to be an even multiple
        # of ti_mod_factor!
        if s.animate and np.mod(s.nt, s.animate_mod_factor) == 0:
            s.make_frame()
        if s.animate and s.ALL_THE_WAY_BACK:
            s.make_frame()
            print('plot.ti(): We are all the way back. s.x_s = %s'%s.x_s)

    def final(self, s):
        if s.animate:
            s.animate_directory(
                s.TEMP_DIR, s.DEST_FILE, gif_delay=s.frame_delay)
        s.stats.final('Done')

#---------------------------------------------------------------------------#
#----------------- End of ALL code related to twins class ------------------#
#---------------------------------------------------------------------------#

