from pylab import plot, show, subplots
from scipy.special import hankel2
from numpy import arange, pi, exp, ma, log
from matplotlib import animation

rho_0 = 1  # density of the equilibrium state
u_0 = 1e-5      # velocity amplitude
a = 1./8      # initial cylinder radius

fig, ax = subplots()

r = arange(a, 20*a, a/10) # in units of wavelengths

t_ini = 0        # initial time in units of wave periods

f = 1        # frequency = 1 cycle per period
w = 2*pi*f
c = 1        # one wavelength per period
k  = w/c
ld = 2*pi/k  # = 1

# analytic solutions
# blackstock p. 409

def p_cmpx(r,t):
    return ( 1j*rho_0*c*u_0 * hankel2(0,r*k) / hankel2(1,a*k) ) * exp(1j*w*t)

def p(r,t):
    return p_cmpx(r,t).real

def u_cmpx(r,t):
    return ( u_0 * hankel2(1,r*k) / hankel2(1,a*k) ) * exp(1j*w*t)

def u(r,t):
    return u_cmpx(r,t).real

def Z(r):
    return 1j*w*rho_0*r*log(1/(k*r))

def p_approx(t):
    return (u_cmpx(a,t)*Z(a)).real

# -------------------------

dt = 1
t_change = dt

def animate(t):
    global t_change,r,a
    if t>t_change:
        a = a/2
        t_change += dt
        r = arange(a, 20*a, a/10)
        ax.set_ylim( (-abs(Z(a))*u_0*2,abs(Z(a))*u_0*2) )
        ax.set_xlim( 0,r[-1] )
        fig.canvas.draw()
        lines[0].set_xdata(r)
        lines[1].set_xdata(r)
    timetext.set_text('')
    atext.set_text("tube radius in wavelengths = " + "{:.2e}".format(a))
    lines[0].set_ydata(p(r,t))
    lines[1].set_ydata(u(r,t))
    scatter.set_offsets([[a],[p_approx(t)]])
    return tuple(lines) + (timetext,) + (atext,) + (scatter,)

def init():
    timetext.set_text('')
    atext.set_text('')
    scatter.set_offsets([[],[]])
    for line in lines:
        line.set_ydata(ma.array(line.get_xdata(), mask=True))
    return tuple(lines) + (timetext,) + (atext,) + (scatter,)

lines = []
line = ax.plot(r,p(r,t_ini),label='p')[0]
lines.append(line)
line = ax.plot(r,u(r,t_ini),label='u')[0]
lines.append(line)
scatter = ax.scatter([a],[p_approx(t_ini)],label='p as Z*u',s=10)
ax.set_ylim( (-abs(Z(a))*u_0*2,abs(Z(a))*u_0*2) )
ax.set_xlim( 0,r[-1] )
ax.set_xlabel("r (wavelengths)")
ax.legend(handles=tuple(lines) + (scatter,))
timetext = ax.text(0.1, 0.95, '', transform = ax.transAxes, va='center')
atext = ax.text(0.1, 0.91, '', transform = ax.transAxes, va='center')

ani = animation.FuncAnimation(fig, animate, arange(0, 20, 0.003),
                              repeat=False, init_func=init,
                              interval=25, blit=True)

ani.save("pulsating-tube.mp4", writer="avconv")

show()
