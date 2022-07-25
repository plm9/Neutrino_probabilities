import numpy as np
import matplotlib.pyplot as plt

R = 6371. #km

s12 = np.sin(33.82/180. * np.pi)
c12 = np.cos(33.82/180. * np.pi)

s23 = np.sin(48.3/180. * np.pi)
c23 = np.cos(48.3/180. * np.pi)

s13 = np.sin(8.61/180. * np.pi)
c13 = np.cos(8.61/180. * np.pi)


dm21 = 7.39e-5 # eV^-2
dm32 = 2.449e-3
#dm31 = 2.55e-3 # eV^-2
dm31 = dm32+dm21

d_cp=90/180.  * np.pi

#phi = 

U = np.array( [
                  [c12 * c13,                    s12 * c13,             s13*np.exp(-d_cp*1j)],
                  [-s12*c23-c12*s23*s13*np.exp(d_cp*1j), c12*c23-s12*s23*s13*np.exp(d_cp*1j), s23*c13],
                  [s12*s23-c12*c23*s13*np.exp(d_cp*1j), -c12*s23-s12*c23*s13*np.exp(d_cp*1j), c23*c13]
                ] )

Uc = np.conjugate(U)


#phi = 1.27 * dm21 * l / E

phis = np.empty((3,3), dtype=float)

def P_Nu(cosZenith, E):
  global phis
  l = -2*R*cosZenith
  p = 0

  phis[1, 0] = 2*1.27 * dm21 * l / E
  phis[2, 0] = 2*1.27 * dm31 * l / E
  phis[2, 1] = 2*1.27 * dm32 * l / E

  for i in range(1, 3):
    for j in range(0, i):
      p -= 4 * np.real( U[1, i] * Uc[0, i] * Uc[1, j] * U[0, j] ) * np.sin(phis[i, j] / 2)**2
      p -= 2 * np.imag( U[1, i] * Uc[0, i] * Uc[1, j] * U[0, j] ) * np.sin(phis[i, j])

  return p

def P_antiNu(cosZenith, E):
  global phis
  l = -2*R*cosZenith
  p = 0

  phis[1, 0] = 2*1.27 * dm21 * l / E
  phis[2, 0] = 2*1.27 * dm31 * l / E
  phis[2, 1] = 2*1.27 * dm32 * l / E

  for i in range(1, 3):
    for j in range(0, i):
      p -= 4 * np.real( Uc[1, i] * U[0, i] * U[1, j] * Uc[0, j] ) * np.sin(phis[i, j] / 2)**2
      p -= 2 * np.imag( Uc[1, i] * U[0, i] * U[1, j] * Uc[0, j] ) * np.sin(phis[i, j])

  return p

def main():
  N = 200
  En  = np.logspace(-1, 2, N+1)
  cz = np.linspace(-1, 0.0, N+1)

  Ps = np.empty((N, N), dtype=float)
  for i in range(N):
    for j in range(N):
      E = 0.5*(En[j] + En[j+1])
      c = 0.5*(cz[i] + cz[i+1])
      #E = En[j]
      #c = cz[i]
      Ps[i, j] = P_Nu(c, E) - P_antiNu(c, E)

  plt.figure()
  plt.pcolormesh(En, cz, Ps, cmap="bwr")
  plt.xscale("log")
  plt.colorbar()
  plt.show()


if __name__=="__main__":
  main()