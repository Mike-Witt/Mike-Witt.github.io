#
#   one_bit_visualzations.py
#
# Status
#
#   28 May 2015: I can animate a bit moving around (two dots) on complex plane
#   29 May 2015: test_complex() shows that P does change in the X basis.

SHOWP = True

class vb_complex:

    def __init__(self, alpha, beta, limits=[-1.2, 1.2, -1.2, 1.2]):
        from matplotlib import pyplot as pl
        pl.ion()
        fig = pl.figure()
        ax = fig.gca()
        ax.set_aspect(1.0)
        self.disp_alpha, = ax.plot(
            alpha.real, alpha.imag, marker='o', linestyle='_', color='blue')
        self.disp_beta, = ax.plot(
            beta.real, beta.imag, marker='o', linestyle='_', color='red')
        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])
        ax.plot([limits[0], limits[1]], [0, 0], color='black')
        ax.plot([0, 0], [limits[2], limits[3]], color='black')
        fig.canvas.set_window_title('Time = %.2f' %0)
        pl.draw()
        pl.show()
        self.fig = fig

    def update(self, t, alpha, beta):
        from matplotlib import pyplot as pl
        fig = self.fig
        ax = fig.gca()
        self.disp_alpha.set_xdata(alpha.real)
        self.disp_alpha.set_ydata(alpha.imag)
        self.disp_beta.set_xdata(beta.real)
        self.disp_beta.set_ydata(beta.imag)
        fig.canvas.set_window_title('Time = %.2f' %t)
        pl.draw()
        pl.show()

def test_complex(alpha=None, beta=None, T=50, num_frames=25, sleep_time=.05):
    global SHOWP
    from time import sleep
    from numpy import exp, sqrt

    dt = float(T)/float(num_frames)
    U_alpha = exp(.3*1j*dt)
    U_beta = exp(.2*1j*dt)
    if alpha == None: alpha = 1/sqrt(2)
    if beta == None: beta = 1/sqrt(2)
    length = alpha*alpha.conjugate() + beta*beta.conjugate()
    if abs(length - 1.0) > 10**-10:
        print('WARNING: Length of vector is %s' %length)
    if SHOWP: bit = vb_prob(alpha, beta, x_basis=True)
    else: bit = vb_complex(alpha, beta)

    for n in range(num_frames):
        alpha = U_alpha*alpha 
        beta = U_beta*beta 
        length = alpha*alpha.conjugate() + beta*beta.conjugate()
        if abs(length - 1.0) > 10**-10:
            print('WARNING: Length of vector is %s' %length)
        bit.update(n*dt, alpha, beta)
        sleep(sleep_time)

class vb_real:

    def __init__(self, alpha, beta, limits=[-1.2, 1.2, -1.2, 1.2]):
        from matplotlib import pyplot as pl
        pl.ion()
        fig = pl.figure()
        ax = fig.gca()
        ax.set_aspect(1.0)

        scale=1; HS=.07; color='black'
        pl.arrow(0, 0, float(alpha), float(beta),  head_width=HS*scale,
            head_length=HS*scale, color=color, length_includes_head=True)

        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])
        ax.plot([limits[0], limits[1]], [0, 0], color='black')
        ax.plot([0, 0], [limits[2], limits[3]], color='black')
        fig.canvas.set_window_title('Time = %.2f' %0)
        pl.draw()
        pl.show()
        self.fig = fig
        self.limits = limits

    def update(self, t, alpha, beta):
        from matplotlib import pyplot as pl
        limits = self.limits
        fig = self.fig
        ax = fig.gca()
        pl.cla()
        scale=1; HS=.07; color='black'
        pl.arrow(0, 0, float(alpha), float(beta),  head_width=HS*scale,
            head_length=HS*scale, color=color, length_includes_head=True)

        ax.set_xlim(limits[0], limits[1])
        ax.set_ylim(limits[2], limits[3])
        ax.plot([limits[0], limits[1]], [0, 0], color='black')
        ax.plot([0, 0], [limits[2], limits[3]], color='black')
        fig.canvas.set_window_title('Time = %.2f' %t)
        pl.draw()

def test_real(alpha=None, beta=None, T=None, num_frames=100, sleep_time=.05):
    global SHOWP
    from time import sleep
    from numpy import exp, sqrt, pi, sin, cos

    if T == None: T = 2*pi
    dt = float(T)/float(num_frames)
    theta = 0
    alpha = cos(theta); beta = sin(theta)
    length = alpha*alpha.conjugate() + beta*beta.conjugate()
    if abs(length - 1.0) > 10**-10:
        print('WARNING: Length of vector is %s' %length)
    if alpha.imag != 0: print('WARNING: alpha has a complex part')
    if beta.imag != 0: print('WARNING: beta has a complex part')
    if SHOWP: r_bit = vb_prob(alpha, beta)
    else: r_bit = vb_real(alpha, beta)

    for n in range(num_frames):
        theta += dt
        alpha = cos(theta); beta = sin(theta)
        length = alpha*alpha.conjugate() + beta*beta.conjugate()
        if abs(length - 1.0) > 10**-10:
            print('WARNING: Length of vector is %s' %length)
        if alpha.imag != 0: print('WARNING: alpha has a complex part')
        if beta.imag != 0: print('WARNING: beta has a complex part')
        r_bit.update(n*dt, alpha, beta)
        sleep(sleep_time)

class vb_prob:

    def __init__(self, alpha, beta, x_basis=False):
        from matplotlib import pyplot as pl
        pl.ion()
        fig = pl.figure()
        self.x_basis = x_basis

        if x_basis:
            from sglib import ip, col
            from numpy import sqrt
            alpha_x = ip(col(1,1)/sqrt(2), col(alpha,beta))
            beta_x = ip(col(1,-1)/sqrt(2), col(alpha,beta))
            alpha = alpha_x; beta = beta_x

        prob_plus = alpha*alpha.conjugate()
        prob_minus = beta*beta.conjugate()
        pos_plus = 1
        pos_minus = 2
        width = .5 
        p1 = pl.bar([pos_plus], prob_plus, width, color='blue')
        p2 = pl.bar([pos_minus], prob_minus, width, color='red')
        pl.ylabel('Probability')
        pl.title('Time = %.2f' %0)
        xticks = [0, pos_plus+width/2., pos_minus+width/2., 3.5]
        pl.xticks(xticks, ('', '+Z', '-Z', ''))
        if x_basis: pl.xticks(xticks, ('', '+X', '-X', ''))
        pl.yticks([0, .25, .5, .75, 1.00])
        pl.ylim(0, 1.1)

        fig.canvas.set_window_title('Time = %.2f' %0)
        pl.draw()
        pl.show()
        self.fig = fig

    def update(self, t, alpha, beta, x_basis=False):
        from matplotlib import pyplot as pl
        fig = self.fig
        ax = fig.gca()
        pl.cla()
        x_basis = self.x_basis

        if x_basis:
            from sglib import ip, col
            from numpy import sqrt
            alpha_x = ip(col(1,1)/sqrt(2), col(alpha,beta))
            beta_x = ip(col(1,-1)/sqrt(2), col(alpha,beta))
            alpha = alpha_x; beta = beta_x

        prob_plus = alpha*alpha.conjugate()
        prob_minus = beta*beta.conjugate()
        pos_plus = 1
        pos_minus = 2
        width = .5
        p1 = pl.bar([pos_plus], prob_plus, width, color='blue')
        p2 = pl.bar([pos_minus], prob_minus, width, color='red')
        pl.ylabel('Probability')
        pl.title('Time = %.2f' %0)
        xticks = [0, pos_plus+width/2., pos_minus+width/2., 3.5]
        pl.xticks(xticks, ('', '+Z', '-Z', ''))
        if x_basis: pl.xticks(xticks, ('', '+X', '-X', ''))
        pl.yticks([0, .25, .5, .75, 1.00])
        pl.ylim(0, 1.1)

        fig.canvas.set_window_title('Time = %.2f' %t)
        pl.draw()
        pl.show()

def draw_xyz(ax, v, color='green', draw_box=True):
            import matplotlib.pyplot as pl
            from sglib import draw_vec

            # Draw the axes
            draw_vec(ax, [1,0,0], color='black', label='x')
            draw_vec(ax, [0,1,0], color='black', label='y')
            draw_vec(ax, [0,0,1], color='black', label='z')

            # Draw the spin vectors
            draw_vec(ax, v, label=r'$\vec{v}$', color=color)

            if draw_box: pass
            else: ax.set_axis_off()

            # If you use "equal" then the last parm on position sets size?
            # But this causes it to fail outsize a notebook!
            #ax.set_position([1,1,1,2])

            ax.set_aspect("equal")
            ax.view_init(elev=5, azim=25)
            pl.draw()
            pl.show()

class vb_xyz:
        def __init__(self, v, color='green', draw_box=True, size=[12,5]):
            import matplotlib.pyplot as pl
            from mpl_toolkits.mplot3d import Axes3D
            from sglib import draw_vec

            pl.ion()
            fig = pl.figure(figsize=size)
            ax = fig.add_subplot(121, projection='3d')
            ax2 = fig.add_subplot(122, projection='3d')
            draw_xyz(ax, v, color=color, draw_box=draw_box)
            draw_xyz(ax2, v, color=color, draw_box=draw_box)

def draw_zx(ax, alpha, beta, limits=[-1, 1, -1, 1]):
    from matplotlib import pyplot as pl
    ax.set_aspect(1.0)
    scale=1; HS=.07; color='black'
    pl.arrow(0, 0, float(alpha), float(beta),  head_width=HS*scale,
        head_length=HS*scale, color=color, length_includes_head=True)
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    ax.plot([limits[0], limits[1]], [0, 0], color='black')
    ax.plot([0, 0], [limits[2], limits[3]], color='black')
    pl.draw()
    pl.show()

class vb_xyz_zx:
        def __init__(self, v, color='green', draw_box=True, size=[12,5]):
            import matplotlib.pyplot as pl
            from mpl_toolkits.mplot3d import Axes3D
            from sglib import draw_vec, ip, col
            from numpy import arccos, cos, sin

            pl.ion()
            fig = pl.figure(figsize=size)
            ax = fig.add_subplot(121, projection='3d')
            ax2 = fig.add_subplot(122)
            draw_xyz(ax, v, color=color, draw_box=draw_box)

            # Now I want to find the corresponding state vector (z basis)
            # NOTE: THIS IS DEFINITELY NOT RIGHT. ALL I'M REALLY DOING HERE
            # IS FINDING A STATE (ON THE ZX PLANE) THAT WOULD GIVE THE SAME
            # PROBABILITIES OF MEASUREING Z AS THE GIVEN SPIN VECTOR.
            z = col(0,0,1)
            v = col(v[0],v[1],v[2]) # Make v is a "vector"
            # Both vectors have a length of 1, so their inner product will
            # give the cosine of the angle between them 
            cos_theta = float(ip(z, v))
            theta = arccos(cos_theta)
            alpha = cos(theta/2)
            beta = sin(theta/2)
            draw_zx(ax2, alpha, beta) 

