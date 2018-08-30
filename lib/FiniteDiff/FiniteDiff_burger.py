import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def IniCon(x):
	return np.where (x % 1. < 0.5, np.power(np.sin(2*x*np.pi),2), 0)

def FiniteDiff_burger(nt, nx):
	"""
	Returns the velocity field and distance for 1D non-linear burger's equation 
	input: no. of time points nt
		   no. of space points nx 
		   courant number C
		
	"""
	# Increments
	dt = 1./(nt-1);  dx = 1./(nx-1)
	c = dt/dx
	# Initialise data structures
	u = np.zeros((nx,nt))
	x = np.linspace(0., 1., nx)
	
	# Initial condition
	u[:,0] = IniCon(x)
	# boundary condition
	#Boundary conditions
	u[0,:] = u[nx-1,:] = 0
	# Loop over time
	for n in range(0,nt-1):
		# print dt\dx
		# loop over space
		for i in range(1,nx-1):
			u[i,n+1] = u[i,n]-u[i,n]*c*(u[i,n]-u[i-1,n])
	c = round(c,2)
	return u, x, c


def FiniteDiff_burgelax_windroff(nt, nx):
	def f ( u ):
		"""
		F evaluates the conservation quantity.
		"""
		value = 0.5 * np.power(u,2);
		return value
	def df ( u ):
		"""
		DF evaluates the derivative of the conservation quantity.
		"""
		value = u*1.0;
		return value

	dt = 1./(nt-1);  dx = 1./(nx-1)
	c = dt/dx
	u = np.zeros((nx,nt))
	x = np.linspace(0., 1., nx)
	u[:,0] = IniCon(x)
	u[0,:] = u[nx-1,:] = 0


	for i in range(0,nt-2):
		u[1,i+1] = u[1,i]-\
		0.5*c * ( f(u[2,i]) - f(u[nx-1,i]) )+\
		0.5*c**2 *\
		( df(0.5*(u[2,i] + u[1,i])) * (f(u[2,i]) - f(u[1,i])) -\
		df(0.5*(u[1,i] + u[nx-1,i])) * (f(u[1,i]) - f(u[nx-1,i])) );

		u[1:nx-2,i+1] = u[1:nx-2,i]-\
		0.5*c * ( f(u[2:nx-1,i]) - f(u[0:nx-3,i]) )+\
		0.5*c**2 *\
		( df(0.5*(u[2:nx-1,i] + u[1:nx-2,i])) * (f(u[2:nx-1,i]) - f(u[1:nx-2,i])) -\
		df(0.5*(u[1:nx-2,i] + u[0:nx-3,i])) * (f(u[1:nx-2,i]) - f(u[0:nx-3,i])) );

		u[nx-1,i+1] = u[nx-1,i]- \
		0.5*c * ( f(u[0,i]) - f(u[nx-2,i]) )+ \
		0.5*c**2 *\
		( df(0.5*(u[1,i] + u[nx-1,i])) * (f(u[1,i]) - f(u[nx-1,i])) -\
		df(0.5*(u[nx-1,i] + u[nx-2,i])) * (f(u[nx-1,i]) - f(u[nx-2,i])) );

	c = round(c,2)
	return u, x, c


def plotburger(ax,u,x,nt,title):
	"""
	Plots the 1D velocity field
	"""

	color=iter(cm.jet(np.linspace(0,1,nt/1)))
	for i in range(0,nt,1):
		ax.plot(x,u[:,i],c=next(color))
		plt.xlabel('x (m)')
		plt.ylabel('u (m\s)')
		plt.ylim([0,1])
		plt.title(title)

fig = plt.figure(facecolor='White',figsize=(16, 10)); 

ax = plt.subplot(2,2,1);
u,x,c = FiniteDiff_burger(160, 20)
plotburger(ax,u,x,160,'(a): FTBS c='+str(c)+'m\s,  nt=160,  nx=20 ')

ax = plt.subplot(2,2,2);
u,x,c= FiniteDiff_burger(160, 40)
plotburger(ax,u,x,160,'(b): FTBS c='+str(c)+'m\s,  nt=160,  nx=40 ')

ax = plt.subplot(2,2,3);
u,x,c = FiniteDiff_burgelax_windroff(160, 20)
plotburger(ax,u,x,160,'(c): LXWF c='+str(c)+'m\s,  nt=160,  nx=20 ')

ax = plt.subplot(2,2,4);
u,x,c= FiniteDiff_burgelax_windroff(160, 160)
plotburger(ax,u,x,160,'(d): LXWF c='+str(c)+'m\s,  nt=160,  nx=160 ')


plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95, wspace=0.1, hspace=0.2);
plt.savefig('Burger_test.png', format='png', dpi=600)
plt.show()
