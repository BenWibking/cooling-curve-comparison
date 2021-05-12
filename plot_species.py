import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import colors

if __name__ == '__main__':
    dens, temp, P, tcool, Ne, HI, HII, HeI, HeII, HeIII, converged = np.loadtxt("grackle_curve.txt", unpack=True, skiprows=3)
    dens_FIRE, temp_FIRE, P_FIRE, tcool_FIRE, Ne_FIRE, HI_FIRE, HII_FIRE, HeI_FIRE, HeII_FIRE, HeIII_FIRE, converged_FIRE = np.loadtxt("gizmo_spcool_FG2009_curve.txt", unpack=True, skiprows=3)

    lw = 1.5
    plt.figure(figsize=(4,4))

    plt.plot(dens, HI, linewidth=lw, color='black', label=r"HI")
    plt.plot(dens_FIRE, HI_FIRE, linewidth=lw, color='black', linestyle='dashed')

    plt.plot(dens, HII, linewidth=lw, color='blue', label=r"HII")
    plt.plot(dens_FIRE, HII_FIRE, linewidth=lw, color='blue', linestyle='dashed')

    plt.plot(dens, HeI, linewidth=lw, color='red', label=r"HeI")
    plt.plot(dens_FIRE, HeI_FIRE, linewidth=lw, color='red', linestyle='dashed')

    plt.plot(dens, HeII, linewidth=lw, color='orange', label=r"HeII")
    plt.plot(dens_FIRE, HeII_FIRE, linewidth=lw, color='orange', linestyle='dashed')

    plt.plot(dens, HeIII, linewidth=lw, color='green', label=r"HeIII")
    plt.plot(dens_FIRE, HeIII_FIRE, linewidth=lw, color='green', linestyle='dashed')

    plt.xlabel(r"gas density ($m_H$ cm$^{-3}$)")
    plt.ylabel(r"species number per H")
    plt.xlim(1e-7, 1e5)
    plt.ylim(1e-6, 2.0)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig("species_fraction_comparison.pdf")
