from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
BIRTHRATEM = 1/14 # Rate of mosquito births
DEATHRATEM = 1/14 # Rate of mosquito deaths
MINFH = 0.4 # Prob. biting mosquito will infect human
HINFM = 0.5 # Prob. biting mosquito will be infected
HINFH = 0.01 # Rate of sexual transmission (includes sex rate)
BITERATE = 0.5 # Rate of mosquito bites
DPROTECT = 0.55 # Protection factor offered by dengue
HRECOVER = 1/10 # Human recovery rate


def deriv(t,y):
	Shu,Shd,Ih,Rh,Sm,Im = y 
	dShu = -MINFH*Shu*Im -HINFH*Shu*Ih
	dShd = DPROTECT*(-MINFH*Shd*Im -HINFH*Shd*Ih)
	dIh = MINFH*Shu*Im+HINFH*Shu*Ih + DPROTECT*(MINFH*Shd*Im+HINFH*Shd*Ih) - HRECOVER*Ih
	dRh = HRECOVER*Ih
	dSm = BIRTHRATEM*(Sm+Im)-DEATHRATEM*Sm-HINFM*Ih*Sm
	dIm = HINFM*Ih*Sm - DEATHRATEM*Im
	return np.array((dShu,dShd,dIh,dRh,dSm,dIm))

t, y = solve_ivp(deriv,(0,100),(100,50,10,0,200,0))

plt.plot(t,y)
plt.show()
