from Meshoid import Meshoid
from load_from_snapshot import load_from_snapshot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import colors


# physical constants
G_cgs = 6.6725985e-8  # cgs
boltzmann_cgs = 1.380658e-16  # erg K^{-1}
m_H = 1.6733e-24    # g
gamma = 5./3.   # assumed constant

# conversion factor from specific energy in code units to specific energy in cgs
seconds_per_year = 3.154e7   # s
solarmass_cgs = 1.989e33     # g  (1 Msun in g)
kpc_in_cm = 3.085678e21      # cm (1 kpc in cm)
pc_in_cm = kpc_in_cm / 1.0e3  # cm (1 pc in cm)
solarluminosity_cgs = 3.83e33  # ergs s^{-1}

unittime_cgs = 3.08568e16    # s  (0.976 Gyr in s)
unitlength_cgs = 3.085678e21  # cm (1 kpc in cm)
unitmass_cgs = 1.989e43      # g  (1e10 Msun in g)
unitvelocity_cgs = unitlength_cgs / unittime_cgs
unitenergypermass_cgs = unitlength_cgs**(2) * unittime_cgs**(-2)  # erg g^{-1}
unitdensity_cgs = unitmass_cgs * unitlength_cgs**(-3)  # g cm^{-3}
unitdensity_per_H = unitdensity_cgs / m_H
unitenergydensity_cgs = unitmass_cgs * \
    unitvelocity_cgs**(2) * unitlength_cgs**(-3)
unitbfield_cgs = 1.0  # gauss


# plot resolution
fig_dpi = 300


def load_hydro_data(filename):
    pdata = {}
    for field in "Masses", "Coordinates", "SmoothingLength", "Velocities", "InternalEnergy", "GrackleHI", "GrackleHII", "GrackleHeII", "GrackleHeIII", "MagneticField":
        pdata[field] = load_from_snapshot(field, 0, filename)

    pdata["Metallicity"] = load_from_snapshot("Metallicity", 0, filename)
    return pdata


def load_dm(filename):
    dmdata = {}
    for field in "Masses", "Coordinates":
        dmdata[field] = load_from_snapshot(field, 1, filename)
    if dmdata is None:
        return None
    else:
        return dmdata


def load_disk(filename):
    dmdata = {}  # old disk stars
    for field in "Masses", "Coordinates":
        dmdata[field] = load_from_snapshot(field, 2, filename)
    if dmdata is None:
        return None
    else:
        return dmdata


def load_bulge(filename):
    dmdata = {}  # old bulge stars
    for field in "Masses", "Coordinates":
        dmdata[field] = load_from_snapshot(field, 3, filename)
    if dmdata is None:
        return None
    else:
        return dmdata


def load_stars(filename):
    stardata = {}
    for field in "Masses", "Coordinates", "StellarFormationTime":
        stardata[field] = load_from_snapshot(field, 4, filename)
    if stardata is None:
        return None
    else:
        return stardata


def compute_temperature(pdata):
    # assume solar metallicity (Grackle assumes this)
    y_Helium = 0.23
    z_metals = 0.02
    x_H = 1.0 - y_Helium - z_metals

    H_term = x_H
    He_term = y_Helium/4.0
    # z_term = z_metals/2.0   # this term is ignored

    #e_term_fullyionized = x_H + y_Helium*(2)/4.0
    #mu_fullyionized = 1.0 / (H_term + He_term + e_term_fullyionized)

    HI = pdata['GrackleHI']
    HII = pdata['GrackleHII']
    nH = (HI + HII)
    HeII = pdata['GrackleHeII']
    HeIII = pdata['GrackleHeIII']
    ne = HII + HeII + 2.0*HeIII
    e_term = ne / nH  # free electrons per hydrogen

    mu = 1.0 / (H_term + He_term + e_term)
    mean_molecular_weight = mu * m_H
    InternalEnergy = pdata["InternalEnergy"] * unitenergypermass_cgs
    T = (mean_molecular_weight / boltzmann_cgs) * (gamma-1) * InternalEnergy
    return T


def compute_pressure(pdata, density):
    InternalEnergy = pdata["InternalEnergy"] * unitenergypermass_cgs
    dens = density * unitdensity_cgs
    P = (gamma-1) * dens * InternalEnergy  # erg cm^-3
    return P


def apply_radius_cut(pdata, T, rmax=40., zmax=1000.):
    coords = pdata["Coordinates"]
    x = coords[:, 0]
    y = coords[:, 1]
    z = coords[:, 2]
    R = np.sqrt(x*x + y*y)
    radius_cut = (R < rmax)
    height_cut = (z < zmax)
    mask = np.logical_and(radius_cut, height_cut)

    ndata = {}
    ndata['Coordinates'] = coords[mask]
    ndata['Masses'] = pdata['Masses'][mask]
    if 'SmoothingLength' in ndata.keys():
        ndata['SmoothingLength'] = pdata['SmoothingLength'][mask]
    ndata['Velocities'] = pdata['Velocities'][mask]
    ndata['MagneticField'] = pdata['MagneticField'][mask]
    ndata['InternalEnergy'] = pdata['InternalEnergy'][mask]
    if T is not None:
        temp = T[mask]
    else:
        temp = None

    return ndata, temp


def save_phase_plot(input_dens, temp, filename):
    # make phase plot! (temperature vs n_H)
    dens = input_dens * unitdensity_per_H

    lognHmin = -7.5
    lognHmax = 5.5
    logTmin = np.log10(3.0)
    logTmax = 8.5

    h = plt.hist2d(np.log10(dens), np.log10(temp),
                   bins=100, norm=colors.LogNorm(),
                   range=[[lognHmin, lognHmax], [logTmin, logTmax]])

    # unclear what units of this are
    plt.colorbar(label=r'proportional to mass')
    plt.xlabel(r'$\log_{10}$ density ($n_{H}$ [cm$^{-3}$])')
    plt.ylabel(r'$\log_{10}$ temperature (K)')
    plt.savefig(filename, dpi=fig_dpi)
    plt.close()


def plot_stars_on_axis(ax, star_coords, rmax=10., s=0.05):
    if star_coords is not None:
        x = star_coords[:, 0]
        y = star_coords[:, 1]
        r = np.sqrt(x**2 + y**2)
        ax.scatter(x[r < rmax], y[r < rmax], s=s, color='black')


def save_slice_plot(mesh, field, filename, colorbar_label="",
                    bfield=None,
                    star_coords=None, rmax=10., plane='z', vmin=10., vmax=1.0e7):
    res = 800
    x = y = np.linspace(-rmax, rmax, res)
    X, Y = np.meshgrid(x, y)

    slice = mesh.Slice(field, center=np.array(
        [0, 0, 0]), size=2*rmax, res=res, plane=plane)

    fig, ax = plt.subplots(figsize=(6, 6))
    p = ax.imshow(slice, cmap='viridis',
                  extent=[x.min(), x.max(), y.min(), y.max()],
                  interpolation='nearest',
                  origin='lower',
                  aspect='equal',
                  norm=colors.LogNorm(vmin=vmin, vmax=vmax))

    if bfield is not None:
        bfield_res = 40
        bfield_slice = mesh.Slice(bfield,
                                  center=np.array([0, 0, 0]),
                                  size=2*rmax,
                                  res=bfield_res,
                                  plane=plane)

        x = y = np.linspace(-rmax, rmax, bfield_res)
        bfield_x = bfield_slice[:, :, 0]
        bfield_y = bfield_slice[:, :, 1]
        bfield_z = bfield_slice[:, :, 2]
        bfield_norm = np.sqrt(bfield_x**2 + bfield_y**2 + bfield_z**2)
        bfield_x /= bfield_norm
        bfield_y /= bfield_norm
        # plot B-fields
        ax.quiver(x, y, bfield_x, bfield_y, angles='xy', pivot='middle')

    if star_coords is not None:
        plot_stars_on_axis(ax, star_coords, rmax=rmax)

    ax.set_aspect('equal')
    fig.colorbar(p, label=colorbar_label)
    ax.set_xlabel("X (kpc)")
    ax.set_ylabel("Y (kpc)")
    plt.savefig(filename, dpi=fig_dpi)
    plt.close()


def compute_density_projection(mesh,
                               rmax=10.,
                               plane='z',
                               projection_length=1.0,  # code units [kpc]
                               res=1000):
    unitlength_pc = 1.0e3  # == 1 kpc
    unitmass_msun = 1.0e10
    unitsurfacedensity_msun_persq_pc = unitmass_msun * unitlength_pc**(-2)

    proj_dens = unitsurfacedensity_msun_persq_pc * mesh.SurfaceDensity(mesh.m,
                                                                       center=np.array(
                                                                           [0, 0, 0]),
                                                                       size=2*rmax, res=res,
                                                                       plane=plane,
                                                                       zmax=0.5*projection_length)

    proj_dens = np.swapaxes(proj_dens, 0, 1)
    x = y = np.linspace(-rmax, rmax, res)
    return proj_dens, x, y


def compute_bfield_slice(mesh, bfield,
                         rmax=10.,
                         plane='z',
                         res=40):
    x = y = np.linspace(-rmax, rmax, res)
    bfield_res = res
    bfield_slice = mesh.Slice(bfield,
                              center=np.array([0, 0, 0]),
                              size=2*rmax,
                              res=bfield_res,
                              plane=plane)

    x = y = np.linspace(-rmax, rmax, bfield_res)
    return bfield_slice, x, y


def save_density_projection_plot(proj_dens, x, y, filename,
                                 vmin=1.0,
                                 vmax=1.0e3,
                                 bfield_slice=None,
                                 bfield_xcoord=None, bfield_ycoord=None,
                                 colorbar_label=r"$\Sigma_{gas}$ $(\rm M_\odot\,pc^{-2})$"):
    proj_dens[proj_dens < vmin] = vmin  # set zeros to vmin

    fig, ax = plt.subplots(figsize=(5, 5))
    im = ax.imshow(proj_dens, cmap='viridis',
                   extent=[x.min(), x.max(), y.min(), y.max()],
                   interpolation='nearest',
                   origin='lower',
                   aspect='equal',
                   norm=colors.LogNorm(vmin=vmin, vmax=vmax))

    if bfield_slice is not None:
        bfield_x = bfield_slice[:, :, 0]
        bfield_y = bfield_slice[:, :, 1]
        bfield_z = bfield_slice[:, :, 2]
        bfield_norm = np.sqrt(bfield_x**2 + bfield_y**2 + bfield_z**2)
        bfield_x /= bfield_norm
        bfield_y /= bfield_norm
        # plot B-fields
        ax.quiver(bfield_xcoord, bfield_ycoord, bfield_x, bfield_y,
                  angles='xy', pivot='middle', color='white')
        # headwidth=1, headlength=2)

    ax.set_aspect('equal')
    #fig.colorbar(p, label=colorbar_label)

    location = 'right'
    cbar = fig.colorbar(im, ax=[ax], label=colorbar_label,
                        location=location, pad=0.09, fraction=0.08, shrink=0.9)
    cbar.ax.yaxis.set_ticks_position(location)

    ax.set_xlabel("X (kpc)")
    ax.set_ylabel("Y (kpc)")
    plt.savefig(filename, dpi=fig_dpi)
    plt.close()


def save_zdensity_projection_plot(mesh, filename, star_coords=None,
                                  vmin=1.0,
                                  vmax=1.0e3,
                                  bfield=None,
                                  rmax=10.,
                                  plane='y',
                                  projection_length=1.0,  # code units [kpc]
                                  colorbar_label=r"y-average density $(\rm H\,cm^{-3})$"):
    res = 1000
    x = y = np.linspace(-rmax, rmax, res)
    X, Y = np.meshgrid(x, y)

    unitlength = unitlength_cgs
    unitmass = unitmass_cgs / m_H
    unitsurfacedensity = unitmass * unitlength**(-2)

    proj_dens = unitsurfacedensity * mesh.SurfaceDensity(mesh.m,
                                                         center=np.array(
                                                             [0, 0, 0]),
                                                         size=2*rmax, res=res,
                                                         plane=plane,
                                                         zmax=0.5*projection_length)
    proj_dens = np.swapaxes(proj_dens, 0, 1)

    mean_dens = proj_dens / (projection_length * unitlength)
    mean_dens[mean_dens < vmin] = vmin  # set zeros to vmin

    fig, ax = plt.subplots(figsize=(6, 6))
    p = ax.imshow(mean_dens, cmap='viridis',
                  extent=[x.min(), x.max(), y.min(), y.max()],
                  interpolation='nearest',
                  origin='lower',
                  aspect='equal',
                  norm=colors.LogNorm(vmin=vmin, vmax=vmax))

    ax.set_aspect('equal')
    fig.colorbar(p, label=colorbar_label)
    ax.set_xlabel("X (kpc)")
    ax.set_ylabel("Z (kpc)")
    plt.title(r"average density for $|\Delta y| \, < \, 500$ pc")
    plt.savefig(filename, dpi=fig_dpi)
    plt.close()


def compute_mesh(pdata):
    if 'SmoothingLength' in pdata.keys():
        return Meshoid(pdata["Coordinates"], pdata["Masses"], pdata["SmoothingLength"])
    else:
        return Meshoid(pdata["Coordinates"], pdata["Masses"])


def compute_rzavg(mesh, field, weights=None, rmax=17.5, nangles=128, res=200, transform_field=None):
    """compute the azimuthal (axisymmetric) average of scalar field 'field'."""

    angles = np.linspace(0., 2.0*np.pi, nangles, endpoint=False)
    r = np.linspace(0., rmax, res)
    z = np.linspace(-rmax/2.0, rmax/2.0, res)
    rzavg = []

    # for mass weighting, weights=mesh.Density()
    if weights is None:
        weights = np.ones_like(field)

    if field is not None:
        slices = np.empty((nangles, res, res))
        weight_arr = np.empty((nangles, res, res))

        for k, angle in enumerate(angles):
            this_field = field
            if transform_field is not None:
                this_field = transform_field(angle)

            slices[k, :, :] = mesh.Slice(weights*this_field, center=np.array(
                [0.5*rmax, 0, 0]), size=rmax, res=res, plane='y', rot=angle)
            weight_arr[k, :, :] = mesh.Slice(weights, center=np.array(
                [0.5*rmax, 0, 0]), size=rmax, res=res, plane='y', rot=angle)

        rzavg = np.nanmean(slices, axis=0)
        rzavg *= 1.0/np.nanmean(weight_arr, axis=0)

    return rzavg, r, z


def save_rzavg_plot(filename, avg_slices, r, z,
                    colorbar_label="",
                    vmin=[1.0e-6], vmax=[1.0],
                    norm=[colors.LogNorm, colors.LogNorm],
                    cmaps=['viridis', 'viridis'],
                    do_colorbar=True,
                    contour_levels=[None, None]):

    figsize = (9, 4)
    if len(avg_slices) == 1:
        figsize = (5, 5)
    fig = plt.figure(figsize=figsize)
    ncols = len(avg_slices)
    if ncols == 1:
        width_ratios = [1]
    elif ncols == 2:
        width_ratios = [1, 1]
    gs = fig.add_gridspec(
        1, ncols, width_ratios=width_ratios, hspace=0, wspace=0)
    axes = gs.subplots(sharex='col', sharey='row')
    if len(avg_slices) == 1:
        axes = [axes]

    locations = ['left', 'right']
    pad = [0.09, 0.09]
    for i, (ax, avg_slice, location) in enumerate(zip(axes, avg_slices, locations)):
        this_r = r.copy()

        if (i % 2 == 0) and (len(avg_slices) > 1):
            avg_slice = avg_slice[:, ::-1]
            this_r *= -1.0

        im = ax.imshow(avg_slice,
                       extent=[this_r.min(), this_r.max(), z.min(), z.max()],
                       interpolation='bilinear',
                       origin='lower',
                       aspect='equal',
                       norm=norm[i](vmin=vmin[i], vmax=vmax[i]),
                       cmap=cmaps[i])

        if contour_levels[i] is not None:
            ax.contour(this_r[::-1], z, avg_slice, contour_levels[i])

        ax.set_xlabel("R (kpc)")

        if do_colorbar:
            cbar = fig.colorbar(im, ax=[ax], label=colorbar_label[i],
                                location=location, pad=pad[i], fraction=0.08, shrink=0.9)
            cbar.ax.yaxis.set_ticks_position(location)
        else:
            ax.set_ylabel("z (kpc)")

    if len(axes) > 1:
        for ax in axes.flat:
            ax.label_outer()
            # ax.set_yticks([])

    plt.savefig(filename, dpi=fig_dpi)
    plt.close()


def compute_poloidal_toroidal(pdata):
    # decompose magnetic field into poloidal, toroidal components
    # [https://en.wikipedia.org/wiki/Toroidal_and_poloidal_coordinates]

    coords = pdata['Coordinates']
    x = coords[:, 0]
    y = coords[:, 1]
    z = coords[:, 2]
    r = np.sqrt(x*x + y*y + z*z)
    theta = np.arctan2(y, x)
    phi = np.arcsin(z/r)

    bfield = pdata['MagneticField'] * unitbfield_cgs  # gauss
    bx = bfield[:, 0]
    by = bfield[:, 1]
    bz = bfield[:, 2]

    thetahat_x = -np.sin(theta)
    thetahat_y = np.cos(theta)
    btheta = bx * thetahat_x + by * thetahat_y

    phihat_x = -np.sin(phi) * np.cos(theta)
    phihat_y = -np.sin(phi) * np.sin(theta)
    phihat_z = np.cos(phi)
    bphi = bx * phihat_x + by * phihat_y + bz * phihat_z

    poloidal = np.abs(bphi)
    toroidal = np.abs(btheta)
    total = np.sqrt(bx**2 + by**2 + bz**2)

    return poloidal, toroidal, total
