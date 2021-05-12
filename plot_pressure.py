from gizmo_analysis import compute_pressure, compute_temperature, load_hydro_data, compute_mesh, m_H, unitdensity_cgs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import colors

kboltz = 1.3806504e-16

def find_unstable_phase(dens, temp, P):
    dP_dn = np.diff(P/kboltz) / np.diff(dens)

    # detect sign change of dP/dn to find unstable phase
    signchange = np.sign(dP_dn[:-1]*dP_dn[1:]) # negative iff sign change between i, i+1
    signchange_at = np.where(signchange == -1)
    density_mid = 0.5*(dens[:-1]+dens[1:])
    density_signchange = 0.5*(density_mid[:-1]+density_mid[1:])
    density_signchange_at = density_signchange[signchange_at]
    temp_signchange_at = np.interp(density_signchange_at, dens, temp)
    pressure_signchange_at = np.interp(density_signchange_at, dens, P) / kboltz

    # restrict to phases below 10^4 K (there are sometimes numerically-spurious sign changes above this)
    neutral_mask = (temp_signchange_at < 1.0e4)
    temp_unm_bracket = temp_signchange_at[neutral_mask]
    dens_unm_bracket = density_signchange_at[neutral_mask]
    pres_unm_bracket = pressure_signchange_at[neutral_mask]

    print(f"unstable neutral phase:")
    print(f"\ttemperature: [{temp_unm_bracket[0]}, {temp_unm_bracket[1]}] K")
    print(f"\tdensity: [{dens_unm_bracket[0]}, {dens_unm_bracket[1]}] H cm^-3")
    print(f"\tpressure: [{pres_unm_bracket[0]}, {pres_unm_bracket[1]}] cm^-3 K")
    return temp_unm_bracket, dens_unm_bracket, pres_unm_bracket

if __name__ == '__main__':
    dens, temp, P, tcool, Ne, HI, HII, HeI, HeII, HeIII, converged = np.loadtxt("grackle_curve.txt", unpack=True, skiprows=3)
    dens_FIRE, temp_FIRE, P_FIRE, tcool_FIRE, Ne_FIRE, HI_FIRE, HII_FIRE, HeI_FIRE, HeII_FIRE, HeIII_FIRE, converged_FIRE = np.loadtxt("gizmo_spcool_FG2009_curve.txt", unpack=True, skiprows=3)

    temp_unm_bracket, dens_unm_bracket, pres_unm_bracket = find_unstable_phase(dens, temp, P)
    find_unstable_phase(dens_FIRE, temp_FIRE, P_FIRE)

    grackle_label = r"Grackle (constant $\Gamma_{pe}$)" # (UVB: Haardt \& Madau 2012)
    FIRE_label = r"FIRE-2 (with $G_0 = 1.7$)" #(UVB: Faucher-Giguere et al. 2009)

    plt.figure(figsize=(4,4))
    plt.plot((dens), (P/kboltz), linewidth=2, color='black', label=grackle_label)
    plt.plot((dens_FIRE), (P_FIRE/kboltz), linewidth=2, color='red', label=FIRE_label)
    plt.xlabel(r"gas density ($m_H$ cm$^{-3}$)")
    plt.ylabel(r"pressure (K cm$^{-3}$)")
    plt.xlim(1e-7, 1e5)
    plt.ylim(1e0, 1e4) # to show differences
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig("pressure_curve_comparison.pdf")
