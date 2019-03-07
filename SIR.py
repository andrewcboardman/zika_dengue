from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

Bhh = 0.01/100 # Human-to-human transmission rate
Gh = 0.1 # Human recovery rate
Pi = 0.6 # dengue protection factor
Mv = 0.01
Bvh = 0.5
Bhv = 0.4
A = 0.5
cc = 0.0005

def func(state,t):
	Shu = state[0]
	Shd = state[1]
	Ih = state[2]
	Rh = state[3]
	Nh = Shu + Shd + Ih + Rh
	Sv = state[4]
	Iv = state[5]
	Nv = Sv + Iv

	dShu = - Bhv * A * Iv * Shu / Nh \
		- Bhh * Ih * Shu / Nh
	dShd = - Bhv * A * Iv * Shd * Pi/ Nh \
		- Bhh * Ih * Shd * Pi / Nh \
	dIh = Bhv * A * Iv * (Shu + Pi * Shd) / Nh \
		+ Bhh * Ih * (Shu + Pi * Shd) / Nh \
		- Gh * Ih
	dRh = Gh * Ih
	dSv = - Bvh * A * Sv * Ih / Nh \
		+ Mv * Iv
	dIv = Bvh * A * Sv * Ih / Nh \
		- Mv * Iv
	return (dShu,dShd,dIh,dRh,dSv,dIv)

y0 = [30,70,10,0,100,0]
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





