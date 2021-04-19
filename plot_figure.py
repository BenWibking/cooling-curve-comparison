from gizmo_analysis import compute_pressure, compute_temperature, load_hydro_data, compute_mesh, m_H, unitdensity_cgs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import colors

kboltz = 1.3806504e-16

if __name__ == '__main__':
    dens, temp, P, tcool, Ne, converged = np.loadtxt(
        "grackle_curve.txt", unpack=True, skiprows=3)

    pdata = load_hydro_data("../outputs_mhd/snapshot_1000.hdf5")
    ptemp = compute_temperature(pdata)
    pmesh = compute_mesh(pdata)
    pdensity_code = pmesh.Density()
    pdensity_cgs = pdensity_code * unitdensity_cgs
    pdensity_nHcgs = pdensity_cgs / m_H
    ppressure_cgs = compute_pressure(pdata, pdensity_code)

    # (UVB: Haardt \& Madau 2012)
    grackle_label = r"Grackle (constant $\Gamma_{pe}$)"
    # (UVB: Faucher-Giguere et al. 2009)
    FIRE_label = r"FIRE-2 (with $G_0 = 1.7$)"

    fig = plt.figure(figsize=(8, 4), constrained_layout=True)
    axes = fig.subplot_mosaic([['temp', 'pressure']])

    # temperature panel
    lognHmin = -7
    lognHmax = 5
    logTmin = np.log10(6.0)
    logTmax = np.log10(2.0e6)
    range = [(lognHmin, lognHmax), (logTmin, logTmax)]
    axes['temp'].hist2d(np.log10(pdensity_nHcgs), np.log10(ptemp), weights=pdata['Masses'],
                        bins=100, density=True, norm=colors.LogNorm(), range=range)
    axes['temp'].plot(np.log10(dens), np.log10(temp), linewidth=2,
                      color='black', label="thermal equilibrium")
    axes['temp'].set_xlabel(r"$\log_{10}$ density ($m_H$ cm$^{-3}$)")
    axes['temp'].set_ylabel(r"$\log_{10}$ temperature (K)")
    axes['temp'].set_xlim(lognHmin, lognHmax)
    axes['temp'].set_ylim(logTmin, logTmax)
    axes['temp'].legend(loc='lower left')

    # pressure panel
    logPmin = -2
    logPmax = 7
    lognHmin = -7
    lognHmax = 5
    (counts, xedges, yedges, hist_im) = axes['pressure'].hist2d(np.log10(pdensity_nHcgs),
                        np.log10(ppressure_cgs/kboltz), weights=pdata['Masses'],
                        bins=100, density=True, norm=colors.LogNorm(),
                        range=[(lognHmin, lognHmax), (logPmin, logPmax)])
    axes['pressure'].plot(np.log10(dens), np.log10(
        P/kboltz), linewidth=2, color='black', label="thermal equilibrium")
    axes['pressure'].set_xlabel(r"$\log_{10}$ density ($m_H$ cm$^{-3}$)")
    axes['pressure'].set_ylabel(r"$\log_{10}$ pressure (K cm$^{-3}$)")
    axes['pressure'].set_xlim(-7, 5)
    axes['pressure'].set_ylim(logPmin, logPmax)
    axes['pressure'].legend(loc='lower right')

    cbar = fig.colorbar(hist_im, ax=list(axes.values()), label='probability density by mass',
                        location='right', fraction=0.05, shrink=0.9)
    cbar.ax.yaxis.set_ticks_position('right')

    #plt.tight_layout()
    plt.savefig("temperature_pressure_figure.pdf")
