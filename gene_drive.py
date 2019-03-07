from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

K = 0.1
mu = 0.05
def breeding(state,t):
	Mw = state[0]
	Fw = state[1]
	Mg = state[2]
	dMw = 0.5 * K * Mw * Fw / (Mw + Mg) - mu * Mw
	dFw = 0.5 * K * Mw * Fw / (Mw + Mg) - mu * Fw
	dMg = K * Mg * Fw / (Mw + Mg) - mu * Mg
	return np.array((dMw,dFw,dMg))

y0 = [99,100,1]
t = np.arange(0,400)
y = odeint(breeding,y0,t)

plt.plot(t,y[:,0],label='Wild-type male population')
plt.plot(t,y[:,1],label='Wild-type female population')
plt.plot(t,y[:,2],label='Gene-drive male population')
plt.legend()
plt.show()