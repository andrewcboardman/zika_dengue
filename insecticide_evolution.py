from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

K = 0.2
mu = 0.05
insecticide = 2
resist = 1
cc = 0.00025
def breeding(state,t):
	Maa = state[0]
	Faa = state[1]
	MAa = state[2]
	FAa = state[3]
	MAA = state[4]
	FAA = state[5]
	Mtot = Maa + MAa + MAA
	Ftot = Faa + FAa + FAA
	Ntot = Mtot + Ftot
	dMaa = 0.5 * K * Maa * Faa / Mtot \
		+ 0.25 * K * MAa * Faa / Mtot \
		+ 0.25 * K * Maa * FAa / Mtot \
		+ 0.125 * K * MAa * FAa / Mtot \
		- insecticide * mu * Maa - cc*Ntot * Maa
	dFaa = 0.5 * K * Maa * Faa / Mtot \
		+ 0.25 * K * MAa * Faa / Mtot \
		+ 0.25 * K * Maa * FAa / Mtot \
		+ 0.125 * K * MAa * FAa / Mtot \
		- insecticide * mu * Faa -  cc*Ntot*Faa
	dMAa = 0.5 * K * Maa * FAA / Mtot \
		+ 0.25 * K * Maa * FAa / Mtot \
		+ 0.25 * K * MAa * Faa / Mtot \
		+ 0.25 * K * MAa * FAa / Mtot \
		+ 0.25 * K * MAa * FAA / Mtot \
		+ 0.25 * K * MAA * FAa / Mtot \
		+ 0.5 * K * MAA * Faa / Mtot \
		- mu * MAa -  cc*Ntot*MAa
	dFAa = 0.5 * K * Maa * FAA / Mtot \
		+ 0.25 * K * Maa * FAa / Mtot \
		+ 0.25 * K * MAa * Faa / Mtot \
		+ 0.25 * K * MAa * FAa / Mtot \
		+ 0.25 * K * MAa * FAA / Mtot \
		+ 0.25 * K * MAA * FAa / Mtot \
		+ 0.5 * K * MAA * Faa / Mtot \
		- mu * FAa-  cc*Ntot*FAa
	dMAA = 0.125 * K * MAa * FAa / Mtot \
		+ 0.25 * K * MAA * FAa / Mtot \
		+ 0.25 * K * MAa * FAA / Mtot \
		+ 0.5 * K * MAA * FAA / Mtot \
		- mu * MAA-  cc*Ntot*MAA
	dFAA = 0.125 * K * MAa * FAa / Mtot \
		+ 0.25 * K * MAA * FAa / Mtot \
		+ 0.25 * K * MAa * FAA / Mtot \
		+ 0.5 * K * MAA * FAA / Mtot \
		- mu * FAA-  cc*Ntot*FAA
	return np.array((dMaa,dFaa,dMAa,dFAa,dMAA,dFAA))

y0 = [100,100,1,1,0,0]
t = np.arange(0,1000)
y = odeint(breeding,y0,t)

plt.plot(t,y[:,0],label='aa male population')
plt.plot(t,y[:,1],label='aa female population')
plt.plot(t,y[:,2],label='Aa male population')
plt.plot(t,y[:,3],label='Aa female population')
plt.plot(t,y[:,4],label='AA male population')
plt.plot(t,y[:,5],label='AA female population')
plt.plot(t,np.sum(y,axis=1),label='total mosquito population')
plt.legend()
plt.show()