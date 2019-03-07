from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

Bhh = 0.01/100 # Human-to-human transmission rate
Gh = 0.1 # Human recovery rate
Pi = 0.6 # dengue protection factor
Muv = 0.01
Bvh = 0.5
Bhv = 0.4
A = 0.5
cc = 0.0005
K = 0.02
insecticide = 1

def func(state,t):
	Shu = state[0]
	Shd = state[1]
	Ih = state[2]
	Rh = state[3]

	Svaa = state[4]
	Ivaa = state[5]
	Mvaa = state[6]
	SvAa = state[7]
	IvAa = state[8]
	MvAa = state[9]
	SvAA = state[10]
	IvAA = state[11]
	MvAA = state[12]

	Nh = np.sum(state[:4])
	Nv = np.sum(state[4:])
	Fvaa = Svaa + Ivaa
	FvAa = SvAa + IvAa
	FvAA = SvAA + IvAA
	Iv = Ivaa + IvAa + IvAA
	NvM = Mvaa + MvAa + MvAA

	dShu = - Bhv * A * Iv * Shu / Nh \
		- Bhh * Ih * Shu / Nh
	dShd = - Bhv * A * Iv * Shd * Pi/ Nh \
		- Bhh * Ih * Shd * Pi / Nh 
	dIh = Bhv * A * Iv * (Shu + Pi * Shd) / Nh \
		+ Bhh * Ih * (Shu + Pi * Shd) / Nh \
		- Gh * Ih
	dRh = Gh * Ih
	dSvaa = - Bvh * A * Svaa * Ih / Nh \
		+ 0.5 * K * Mvaa * Fvaa / NvM \
		+ 0.25 * K * MvAa * Fvaa / NvM \
		+ 0.25 * K * Mvaa * FvAa / NvM \
		+ 0.125 * K * MvAa * FvAa / NvM \
		- insecticide * Muv * Svaa - cc * Nv * Svaa
	dIvaa = Bvh * A * Svaa * Ih / Nh \
		- insecticide * Muv * Ivaa - cc * Nv * Ivaa
	dMvaa = 0.5 * K * Mvaa * Fvaa / NvM \
		+ 0.25 * K * MvAa * Fvaa / NvM \
		+ 0.25 * K * Mvaa * FvAa / NvM \
		+ 0.125 * K * MvAa * FvAa / NvM \
		- insecticide * Muv * Svaa - cc * Nv * Svaa
	dSvAa = - Bvh * A * Svaa * Ih / Nh \
		+ 0.5 * K * Mvaa * FvAA / NvM \
		+ 0.25 * K * Mvaa * FvAa / NvM \
		+ 0.25 * K * MvAa * Fvaa / NvM \
		+ 0.25 * K * MvAa * FvAa / NvM \
		+ 0.25 * K * MvAa * FvAA / NvM \
		+ 0.25 * K * MvAA * FvAa / NvM \
		+ 0.5 * K * MvAA * Fvaa / NvM \
		- Muv * SvAa -  cc * Nv * SvAa
	dIvAa = Bvh * A * Svaa * Ih / Nh \
		- Muv * IvAa -  cc * Nv * IvAa
	dMvAa = 0.5 * K * Mvaa * FvAA / NvM \
		+ 0.25 * K * Mvaa * FvAa / NvM \
		+ 0.25 * K * MvAa * Fvaa / NvM \
		+ 0.25 * K * MvAa * FvAa / NvM \
		+ 0.25 * K * MvAa * FvAA / NvM \
		+ 0.25 * K * MvAA * FvAa / NvM \
		+ 0.5 * K * MvAA * Fvaa / NvM \
		- Muv * SvAa -  cc * Nv * SvAa
	dSvAA = - Bvh * A * Svaa * Ih / Nh \
		+ 0.125 * K * MvAa * FvAa / NvM \
		+ 0.25 * K * MvAA * FvAa / NvM \
		+ 0.25 * K * MvAa * FvAA / NvM \
		+ 0.5 * K * MvAA * FvAA / NvM \
		- Muv * MvAA - cc * Nv * MvAA
	dIvAA = Bvh * A * Svaa * Ih / Nh \
		- Muv * IvAA -  cc * Nv * IvAA
	dMvAA = - Bvh * A * Svaa * Ih / Nh \
		+ 0.125 * K * MvAa * FvAa / NvM \
		+ 0.25 * K * MvAA * FvAa / NvM \
		+ 0.25 * K * MvAa * FvAA / NvM \
		+ 0.5 * K * MvAA * FvAA / NvM \
		- Muv * MvAA - cc * Nv * MvAA
	return (dShu,dShd,dIh,dRh,dSvaa,dIvaa,dMvaa,dSvAa,dIvAa,dMvAa,dSvAA,dIvAA,dMvAA)

y0 = [30,70,10,0,0,0,0,0,0,0,0,0,0]
t = np.arange(0,10)
y = odeint(func,y0,t)

plt.plot(t,y[:,0],label='Susceptible population (non-dengue)')
plt.plot(t,y[:,1],label='Susceptible population (dengue)')
plt.plot(t,y[:,2],label='Infected population')
plt.plot(t,y[:,3],label='Recovered population')
plt.plot(t,y[:,4],label='Susceptible mosquitos')
plt.plot(t,y[:,5],label='Infected mosquitos')

plt.legend()
plt.show()
