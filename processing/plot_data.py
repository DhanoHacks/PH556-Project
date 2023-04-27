import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import LombScargle
target = "V* 471 Vul"
data = np.genfromtxt("results copy.txt",delimiter=",")
data = data[1:,:]
data = data[data[:, 0].argsort()]
# print(data)

minJD = data[0,0]
meanmag = np.mean(data[:,1])
stdmag = np.std(data[:,1])
print(f"Mean PSF Mag: {meanmag:.2f}")

period = np.linspace(3,40,3801)
frequency = 1/period
ls = LombScargle(data[:,0], data[:,1], data[:,2])
power = ls.power(frequency)
best_frequency = frequency[np.argmax(power)]
best_period = 1/best_frequency

t_fit = np.linspace(np.min(data[:,0])-1,np.max(data[:,0])+1,100)
y_fit = ls.model(t_fit, best_frequency)

print(f"LS best period: {best_period:.1f} days")
plt.plot(frequency, power) 
plt.xlabel(f"Frequency")
plt.ylabel(f"Power")
plt.title(f"Lomb-Scargle Periodogram")
plt.savefig("ls_power_map.png")
# plt.show()
plt.close()

# plt.plot(plotrange,np.real(fit_data))
plt.plot(t_fit-minJD,y_fit)
plt.scatter(data[:,0]-minJD,data[:,1])
plt.errorbar(data[:,0]-minJD,data[:,1],yerr=data[:,2],fmt="o")
plt.xlabel(f"JD since {minJD}")
plt.ylabel(f"PSF Magnitude")
plt.title("Light curve + Lomb-Scargle model fit")
plt.savefig("lightcurve.png")
# plt.show()
plt.close()


# FU mode cepheid, GDR2+external parallaxes, V band, N=374 samples
a = -1.813
a_err = 0.051
b = -2.490
b_err = 0.057
# Mv = a + b log P
Mv = a + b*np.log10(best_period)
print(f"Absolute magnitude of {target}: {Mv:.2f} (V-band)")
color_index = 1.81
print(f"Color index (B-V) of {target} is {color_index:.2f}")
#  g = V + 0.634*(B-V) - 0.108 (Bilir et. al. 2004) - https://ui.adsabs.harvard.edu/abs/2005AN....326..321B/abstract
Mg = Mv + 0.634*color_index - 0.108
print(f"Absolute magnitude of {target}: {Mg:.2f} (g-band)")
mg = meanmag
print(f"Apparent magnitude of {target}: {mg:.2f} (g'-band, close enough to g-band)")
d = 10**((mg-Mg+5)/5)
print(f"Distance modulus of {target}: {d:.0f} parsec")