from gizmo_analysis import compute_pressure, compute_temperature, load_hydro_data, compute_mesh, m_H, unitdensity_cgs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import colors

if __name__ == '__main__':
    dens, temp, P, tcool, Ne, converged = np.loadtxt("grackle_curve.txt", unpack=True, skiprows=3)
    dens_FIRE, temp_FIRE, P_FIRE, tcool_FIRE, Ne_FIRE, converged_FIRE = np.loadtxt("gizmo_spcool_FG2009_curve.txt", unpack=True, skiprows=3)

    grackle_label = r"Grackle (constant $\Gamma_{pe}$)" # (UVB: Haardt \& Madau 2012)
    FIRE_label = r"FIRE-2 (with $G_0 = 1.7$)" #(UVB: Faucher-Giguere et al. 2009)

    lognHmin = -7
    lognHmax = 5
    plt.figure(figsize=(4,4))
    plt.plot((dens), (temp), linewidth=2, color='black', label=grackle_label)
    plt.plot((dens_FIRE), (temp_FIRE), linewidth=2, color='red', label=FIRE_label)
    plt.xlabel(r"gas density ($m_H$ cm$^{-3}$)")
    plt.ylabel("temperature (K)")
    plt.xlim(10.**lognHmin, 10.**lognHmax)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.savefig("temperature_curve_comparison.pdf")
