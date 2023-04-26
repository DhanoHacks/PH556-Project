import matplotlib.pyplot as plt
import numpy as np
from nfft import nfft,ndft,ndft_adjoint
from scipy.signal import argrelextrema
N=80
data = np.genfromtxt("results.txt",delimiter=",")
data = data[1:,:]
data = data[data[:, 0].argsort()]
# print(data)

minJD = data[0,0]
meanmag = np.mean(data[:,1])
stdmag = np.std(data[:,1])
days_shifted = data[:,0]-minJD
mag_shifted = (data[:,1]-meanmag)/stdmag
freqs = ndft_adjoint(days_shifted,mag_shifted,N)

plotrange = np.arange(days_shifted[-1]+1)
# plotrange = days_shifted
print(np.arange(-N//2,N//2)[argrelextrema(np.abs(np.real(freqs)),np.greater)])
plt.plot(np.arange(-N//2,N//2),np.real(freqs))
plt.axhline(0)
plt.savefig("lol.png")
plt.close()
# print(np.real(f),np.real(f/np.abs(f)))

fit_data = np.zeros(shape=np.shape(plotrange))
for j in range(np.shape(plotrange)[0]):
    for k in range(np.shape(freqs)[0]):
        fit_data[j] += freqs[k]*np.exp(-2*np.pi*1j*(-N//2+k)*plotrange[j])
fit_data = fit_data*stdmag+meanmag

# plt.plot(plotrange,np.real(fit_data))
plt.scatter(days_shifted,mag_shifted*stdmag+meanmag)
plt.errorbar(days_shifted,mag_shifted*stdmag+meanmag,yerr=data[:,2],fmt="o")
plt.xlabel(f"JD since {minJD}")
plt.ylabel(f"PSF Magnitude")
plt.title("Light curve of V0471 Vul")
plt.savefig("lightcurve.png")
# plt.show()
plt.close()
print(meanmag)