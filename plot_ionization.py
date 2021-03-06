from gizmo_analysis import compute_pressure, compute_temperature, load_hydro_data, compute_mesh, m_H, unitdensity_cgs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import colors

kboltz = 1.3806504e-16

if __name__ == '__main__':
    dens, temp, P, tcool, Ne, HI, HII, HeI, HeII, HeIII, converged = np.loadtxt("grackle_curve.txt", unpack=True, skiprows=3)
    dens_FIRE, temp_FIRE, P_FIRE, tcool_FIRE, Ne_FIRE, HI_FIRE, HII_FIRE, HeI_FIRE, HeII_FIRE, HeIII_FIRE, converged_FIRE = np.loadtxt("gizmo_spcool_FG2009_curve.txt", unpack=True, skiprows=3)

    ion_threshold = [0.4, 0.5, 0.9]
    for threshold in ion_threshold:
        ion_temp = np.interp(threshold, Ne, temp)
        print(f"ionization temperature: {ion_temp}")

    grackle_label = r"Grackle" # (UVB: Haardt \& Madau 2012)
    FIRE_label = r"FIRE-2" #(UVB: Faucher-Giguere et al. 2009)
    thickness = 2

    plt.figure(figsize=(4,4))
    plt.plot(dens, Ne, linewidth=thickness, color='black', label=grackle_label)
    plt.plot(dens_FIRE, Ne_FIRE, linewidth=thickness, color='red', label=FIRE_label)
    ne_obs = 0.047 # cm^-3
    nH_obs = 0.5 # n(H_tot) = 0.5 cm^-3 assumed in Jenkins (2013)
    e_per_H_obs = ne_obs / nH_obs
    print(f"free electrons observationally-inferred per H nucleon: {e_per_H_obs}")
    plt.scatter([nH_obs], [e_per_H_obs], marker='x', color='black', label='observation (Jenkins 2013)')
    plt.xlabel("gas density ($m_H$ cm$^{-3}$)")
    plt.ylabel("free electrons per H nucleon")
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1.0e-6, 1.3)
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.savefig("ionization_curve.pdf")

    plt.figure(figsize=(4,4))
    plt.plot(Ne, temp, linewidth=thickness, color='black', label=grackle_label)
    plt.plot(Ne_FIRE, temp_FIRE, linewidth=thickness, color='red', label=FIRE_label)
    plt.xlabel("equilibrium ionisation fraction (dimensionless)")
    plt.ylabel("temperature (K)")
    plt.yscale('log')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig("ionization_threshold_curve.pdf")