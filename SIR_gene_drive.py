from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

Bhh = 0.01/100 # Human-to-human transmission rate
Gh = 0.1 # Human recovery rate
Pi = 0.6 # dengue protection factor
Mv = 0.05
Bvh = 0.5
Bhv = 0.4
A = 0.5
K = 0.5
cc = 0.0005
def func(state,t):
	Shu = state[0]
	Shd = state[1]
	Ih = state[2]
	Rh = state[3]
	Nh = Shu + Shd + Ih + Rh
	Sv = state[4]
	Iv = state[5]
	Mw = state[6]
	Mg = state[7]
	Nv = Sv + Iv + Mw + Mg
	dShu = - Bhv * A * Iv * Shu / Nh \
		- Bhh * Ih * Shu / Nh
	dShd = - Bhv * A * Iv * Shd * Pi / Nh \
		- Bhh * Ih * Shd * Pi / Nh
	dIh = Bhv * A * Iv * (Shu + Pi * Shd) / Nh \
		+ Bhh * Ih * (Shu + Pi * Shd) / Nh \
		- Gh * Ih
	dRh = Gh * Ih
	dSv = - Bvh * A * Sv * Ih / Nh \
		+ 0.5 * K * Mw * (Sv + Iv) / (Mw + Mg) \
		- Mv * Sv - cc * Nv* Sv
	dIv = Bvh * A * Sv * Ih / Nh \
		- Mv * Iv - cc * Nv* Iv
	dMw = 0.5 * K * Mw * (Sv + Iv) / (Mw + Mg) \
		- Mv * Mw - cc * Nv* Mw
	dMg = K * Mg * (Sv + Iv) / (Mw + Mg) \
		- Mv * Mg - cc * Nv* Mg
	return (dShu,dShd,dIh,dRh,dSv,dIv,dMw,dMg)

y0 = [30,70,10,0,0,0,0,0]
t = np.arange(0,400)
y = odeint(func,y0,t)

plt.plot(t,y[:,0],label='Susceptible population (non-dengue)')
plt.plot(t,y[:,1],label='Susceptible population (dengue)')
plt.plot(t,y[:,2],label='Infected population')
plt.plot(t,y[:,3],label='Recovered population')
plt.plot(t,y[:,4],label='Susceptible mosquitos')
plt.plot(t,y[:,5],label='Infected mosquitos')

plt.legend()
plt.show()