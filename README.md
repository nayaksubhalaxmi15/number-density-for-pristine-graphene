# number-density-for-pristine-graphene
Here we have calculated the number density of grapehe without any disorder and impurities.
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
plt.rc('mathtext', fontset='cm')
plt.rc('font', family='serif')




mut = np.arange(0.1, 5+ 0.001, 0.001)  # Corrected range, including 5

# Function for fermi integral
def integrand(x, t, mut_val):
    return (x**(t-1)) / (np.exp(x - mut_val) + 1)

# Compute fermi dirac integral
def compute_fermi_integral(t, mut, integrand):
    f_integral = []
    for mut_val in mut:
        result, error = integrate.quad(integrand, 0, np.inf, args=(t, mut_val))
        f_integral.append(result)
    return np.array(f_integral)

# Fermi integrals
f_1 = compute_fermi_integral(1, mut, integrand)
f_2 = compute_fermi_integral(2, mut, integrand)
f_3 = compute_fermi_integral(3, mut, integrand)
# Function for fermi integral
def integrand1(x, t, mut_val):
    return (x**(t-1)) / (np.exp(x + mut_val) + 1)

def compute_fermi_integral(t, mut, integrand1):
    f_integral = []
    for mut_val in mut:
        result, error = integrate.quad(integrand1, 0, np.inf, args=(t, mut_val))
        f_integral.append(result)
    return np.array(f_integral)

f1_1 = compute_fermi_integral(1, mut, integrand1)
f1_2 = compute_fermi_integral(2, mut, integrand1)
f1_3 = compute_fermi_integral(3, mut, integrand1)
#calculate number density
k_B = 1.38e-23  # Boltzmann constant in Joules per Kelvin
h = 6.63e-34    # Planck's constant in Joules seconds
e=1.6*10**(-19)
T=60 # K
v_f=1.11*10**6
n1=(8*np.pi * k_B**2 * T**2)/(h**2 * v_f**2)*(f_2 )*10**(-4)
#- f1_2)
n2=(8*np.pi * k_B**2 * T**2)/(h**2 * v_f**2)*(f1_2)*10**(-4)
nt=(8*np.pi * k_B**2 * T**2)/(h**2 * v_f**2)*(f_2 - f1_2)*10**(-4)

#plotting
plt.figure(figsize=(6.5, 4.5))
fig, ax = plt.subplots()
plt.plot(mut, n1,color='blue', label='$\ n_e $',linewidth='3')
plt.plot(mut, n2 ,color='red',label='$\ n_h $',linewidth='3')
plt.plot(mut, nt ,color='green',label='$n=\ n_e $-$ \ n_h $',linewidth='3')
plt.ylim(9*10**8, 10**11)
plt.xlim(0.1, 5)
ax. yaxis.set_ticks_position('both')
ax.tick_params(axis="y" , direction="in" , pad=1)
ax.tick_params(axis="x" , direction="in" , pad=1)
ax.tick_params(axis="y", which='minor', direction="in", pad=1)
ax.tick_params(axis="x", which='minor', direction="in", pad=1)
plt.tick_params(axis="both", which='major', length=12)
plt.tick_params(axis="both", which='minor', length=6)
plt.tick_params(width=3)
plt.tick_params(which='minor', width=2.5) 
plt.xlabel('X-axis', fontsize=16)
plt.ylabel('Y-axis', fontsize=16)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax.spines['bottom'].set_linewidth(3)  # Thicken x-axis
ax.spines['top'].set_linewidth(3) 
ax.spines['right'].set_linewidth(3) 
ax.spines['left'].set_linewidth(3)    # Thicken y-axis
#ax.text(4,80000000000,'(a)',color='black',ha='center', va='center',fontsize=14)
#plt.plot(mut,n)
plt.xlabel(' Normalized chemical potential $\mu/k_BT$')
plt.ylabel('Net number density $ \\tilde {n} (cm^{-2})$')
plt.yscale('log')
plt.xscale('log')
plt.legend(fontsize=18, edgecolor='Black',frameon= False)
#plt.legend(linewidth=2)
np.savetxt('nt.dat', np.vstack((mut, nt)).T, delimiter='\t' )
plt.show()
