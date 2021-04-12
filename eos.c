#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"
//#include "../kernel.h"

/*! Routines for gas equation-of-state terms (collects things like calculation of gas pressure)
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */


/* return an estimate of the Hydrogen molecular fraction of gas, intended for simulations of e.g. molecular clouds, galaxies, and star formation */
double Get_Gas_Molecular_Mass_Fraction(struct species_data *data, double temperature, double neutral_fraction, double free_electron_ratio, double urad_from_uvb_in_G0)
{
    /* if tracking chemistry explicitly, return the explicitly-evolved H2 fraction */
    
#if (COOL_GRACKLE_CHEMISTRY >= 2) // Use GRACKLE explicitly-tracked H2 [using the molecular network if this is valid]
    return DMIN(1,DMAX(0, data->grH2I + data->grH2II)); // include both states of H2 tracked
#endif
    
#if defined(COOL_MOLECFRAC_NONEQM) // use our simple 1-species network for explicitly-evolved H2 fraction
    return DMIN(1, DMAX(0, data->MolecularMassFraction));
#endif


    
#if defined(COOL_MOLECFRAC_LOCALEQM) || defined(COOL_MOLECFRAC_KMT) || defined(COOL_MOLECFRAC_GD) // here are some of the 'fancy' molecular fraction estimators which need various additional properties
    double T=1, nH_cgs=1, Z_Zsol=1, urad_G0=1, xH0=1, x_e=0; // initialize definitions of some variables used below to prevent compiler warnings
    if(temperature > 3.e5) {return 0;} else {T=temperature;} // approximations below not designed for high temperatures, should simply give null
    xH0 = DMIN(DMAX(neutral_fraction,0.),1.); // get neutral fraction [given by call to this program]
    x_e = DMIN(DMAX(free_electron_ratio,0.),2.); // get free electron ratio [number per H nucleon]
    nH_cgs = data->Density*All.cf_a3inv * UNIT_DENSITY_IN_NHCGS; // get nH defined as number of nucleons per cm^3
    Z_Zsol=1; urad_G0=1; // initialize metal and radiation fields. will assume solar-Z and spatially-uniform Habing field for incident FUV radiation unless reset below.
#ifdef METALS
    Z_Zsol = data->Metallicity[0]/All.SolarAbundances[0]; // metallicity in solar units [scale to total Z, since this mixes dust and C opacity], and enforce a low-Z floor to prevent totally unphysical behaviors at super-low Z [where there is still finite opacity in reality; e.g. Kramer's type and other opacities enforce floor around ~1e-3]
#endif
    /* get incident radiation field from whatever module we are using to track it */
#if defined(RT_PHOTOELECTRIC) || defined(RT_LYMAN_WERNER)
    int whichbin = RT_FREQ_BIN_LYMAN_WERNER;
#if !defined(RT_LYMAN_WERNER)
    whichbin = RT_FREQ_BIN_PHOTOELECTRIC; // use photo-electric bin as proxy (very close) if don't evolve LW explicitly
#endif
    urad_G0 = data->Rad_E_gamma[whichbin] * (data->Density*All.cf_a3inv/data->Mass) * UNIT_PRESSURE_IN_CGS / 3.9e-14; // convert to Habing field //
#endif
    urad_G0 += urad_from_uvb_in_G0; // include whatever is contributed from the meta-galactic background, fed into this routine
    urad_G0 = DMIN(DMAX( urad_G0 , 1.e-10 ) , 1.e10 ); // limit values, because otherwise exponential self-shielding approximation easily artificially gives 0 incident field
#endif
        
    
#if defined(COOL_MOLECFRAC_LOCALEQM)
    /* estimate local equilibrium molecular fraction actually using the real formation and destruction rates. expressions for the different rate terms
        as used here are collected in Nickerson, Teyssier, & Rosdahl et al. 2018. Expression for the line self-shielding here
        including turbulent and cell line blanketing terms comes from Gnedin & Draine 2014. below solves this all exactly, using the temperature, metallicity,
        density, ionization states, FUV incident radiation field, and column densities in the simulations. */
    /* take eqm of dot[nH2] = a_H2*rho_dust*nHI [dust formation] + a_GP*nHI*ne [gas-phase formation] + b_3B*nHI*nHI*(nHI+nH2/8) [3-body collisional form] - b_H2HI*nHI*nH2 [collisional dissociation]
        - b_H2H2*nH2*nH2 [collisional mol-mol dissociation] - Gamma_H2^LW * nH2 [photodissociation] - Gamma_H2^+ [photoionization] - xi_H2*nH2 [CR ionization/dissociation] */
    double fH2=0, sqrt_T=sqrt(T), nH0=xH0*nH_cgs, n_e=x_e*nH_cgs, EXPmax=40.; // define some variables for below, including neutral H number density, free electron number, etc.
    double a_Z  = (9.e-19 * T / (1. + 0.04*sqrt_T + 0.002*T + 8.e-6*T*T)) * (0.5*Z_Zsol) * nH_cgs * nH0; // dust formation
    double a_GP = (1.833e-21 * pow(T,0.88)) * nH0 * n_e; // gas-phase formation
    double b_3B = (6.0e-32/sqrt(sqrt_T) + 2.0e-31/sqrt_T) * nH0 * nH0 * nH0; // 3-body collisional formation
    double b_H2HI = (7.073e-19 * pow(T,2.012) * exp(-DMIN(5.179e4/T,EXPmax)) / pow(1. + 2.130e-5*T , 3.512)) * nH0 * (nH0/2.); // collisional dissociation
    double b_H2H2 = (5.996e-30 * pow(T,4.1881) * exp(-DMIN(5.466e4/T,EXPmax)) / pow(1. + 6.761e-6*T , 5.6881)) * (nH0/2.) * (nH0/2.); // collisional mol-mol dissociation
    double G_LW = 3.3e-11 * urad_G0 * (nH0/2.); // photo-dissociation (+ionization); note we're assuming a spectral shape identical to the MW background mean, scaling by G0
    double xi_cr_H2 = (7.525e-16) * (nH0/2.); // CR dissociation (+ionization)
    // can write this as a quadtratic: 0 = x_a*f^2 - x_b*f + x_c, with f = molec mass fraction
    double x_a = (b_3B + b_H2HI - b_H2H2); // terms quadratic in f -- this term can in principle be positive or negative, usually positive
    double x_b = (a_GP + a_Z + 2.*b_3B + b_H2HI + G_LW + xi_cr_H2); // terms linear in f [note sign, pulling the -sign out here] -- positive-definite
    double x_c = (a_GP + a_Z + b_3B); // terms independent of f -- positive-definite
    double y_a = x_a / (x_c + MIN_REAL_NUMBER), y_b = x_b / (x_c + MIN_REAL_NUMBER), z_a = 4. * y_a / (y_b*y_b + MIN_REAL_NUMBER); // convenient to convert to dimensionless variable needed for checking definite-ness
    if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // checking limits of terms for accuracy

    /* now comes the tricky bit -- need to account for the -molecular- self-shielding [depends on fH2, not just the dust external shielding already accounted for */
    double xb0 = a_GP + a_Z + 2.*b_3B + b_H2HI + xi_cr_H2;
    if(fH2 > 1.e-10 && fH2 < 0.99 && G_LW > 0.1*xb0) // fH2 is non-trivial, and the radiation term is significant, so we need to think about molecular self-shielding
    {
        double fH2_min = fH2; // we have just calculated fH2 with -no- molecular self-shielding, so this number can only go up from here
        // calculate a bundle of variables we will need below, to account for the velocity-gradient Sobolev approximation and slab attenuation of G0 //
        double dx_cell = Get_Particle_Size(i) * All.cf_atime; // cell size
        double surface_density_H2_0 = 5.e14 * PROTONMASS, x_exp_fac=0.00085, w0=0.2; // characteristic cgs column for -molecular line- self-shielding
        double surface_density_local = xH0 * data->Density * All.cf_a3inv * dx_cell * UNIT_SURFDEN_IN_CGS; // this is -just- the [neutral] depth through the local cell/slab. that's closer to what we want here, since G0 is -already- attenuated in the pre-processing step!
        double v_thermal_rms = 0.111*sqrt(T); // sqrt(3*kB*T/2*mp), since want rms thermal speed of -molecular H2- in kms
        double dv2=0; int j,k; for(j=0;j<3;j++) {for(k=0;k<3;k++) {double vt = data->Gradients.Velocity[j][k]*All.cf_a2inv; /* physical velocity gradient */
            if(All.ComovingIntegrationOn) {if(j==k) {vt += All.cf_hubble_a;}} /* add hubble-flow correction */
            dv2 += vt*vt;}} // calculate magnitude of the velocity shear across cell from || grad -otimes- v ||^(1/2)
        double dv_turb=sqrt(dv2)*dx_cell*UNIT_VEL_IN_KMS; // delta-velocity across cell
        double x00 = surface_density_local / surface_density_H2_0, x01 = x00 / (sqrt(1. + 3.*dv_turb*dv_turb/(v_thermal_rms*v_thermal_rms)) * sqrt(2.)*v_thermal_rms), y_ss, x_ss_1, x_ss_sqrt, fH2_tmp, fH2_max, Qmax, Qmin; // variable needed below. note the x01 term corrects following Gnedin+Draine 2014 for the velocity gradient at the sonic scale, assuming a Burgers-type spectrum [their Eq. 3]

        fH2_tmp = 1.; // now consider the maximally shielded case, if you had fmol = 1 in the shielding terms
        x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=xb0+y_ss*G_LW; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
        z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant
        fH2_max = DMAX(0,DMIN(1,fH2)); // this serves as an upper-limit for f
        
        if(fH2_max > 1.1*fH2_min)
        {
            fH2_tmp = fH2_max; // re-calculate the maximally-shielded case
            x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=xb0+y_ss*G_LW; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
            z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant
            fH2_tmp=fH2; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=xb0+y_ss*G_LW; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
            fH2_max = fH2; Qmax = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // set the new max fH2, from this, and set the corresponding value of the function we are trying to root-find for

            fH2_tmp = fH2_min; // re-calculate the minimally-shielded case
            x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=xb0+y_ss*G_LW; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
            z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant
            fH2_tmp=fH2; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=xb0+y_ss*G_LW; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
            fH2_min = fH2; Qmin = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // set the new min fH2, from this, and set the corresponding value of the function we are trying to root-find for

            fH2 = exp( (log(fH2_min)*Qmax - log(fH2_max)*Qmin) / (Qmax-Qmin) ); // do a Newton-Raphson step in log[f_H2] space now that we have good initial brackets
            if((fH2_max > 1.5*fH2_min) && (Qmax*Qmin < 0) && (fH2_max > 1.1*fH2)) // have a big enough dynamic range, and bracketing Qmin/max, to make further iteration meaningful
            {
                double f_p=fH2_min, Q_p=Qmin, Q, fH2_new; int iter=0; // define variables for iteration below
                while(1)
                {
                    x_ss_1=1.+fH2*x01; x_ss_sqrt=sqrt(1.+fH2*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=xb0+y_ss*G_LW; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
                    Q = 1 + y_a*fH2*fH2 - y_b*fH2; // update the value of the function we are trying to zero
                    if(iter==0) {if(Q*Q_p>=0) {f_p=fH2_max; Q_p=Qmax;}} // check in case we attempted to bracket from the 'wrong side'
                    if(Q*Q_p >= 0) {break;} // no longer bracketing, end while loop
                    fH2_new = exp( (log(f_p)*Q - log(fH2)*Q_p) / (Q-Q_p) ); f_p=fH2; fH2=fH2_new; Q_p=Q; // update guess and previous values //
                    iter++; // count iterations
                    if(fabs(fH2-f_p) < 0.1*0.5*(f_p+fH2)) {break;} // converged well enough for our purposes!
                    if((y_ss > 0.85) || (y_ss*G_LW < xb0)) {break;} // negligible shielding, or converged to point where external LW is not dominant dissociator so no further iteration needed
                    if((fH2 > 0.95*fH2_max) || (fH2 > 0.99) || (fH2 < 1.e-10) || (fH2 < 1.05*fH2_min) || (iter > 10)) {break;} // approached physical limits or bounds of validity, or end of iteration cycle
                } // end of convergence iteration to find solution for fmol
            } // opening condition for iteration requiring large enough dynamic range, valid bracketing
        } // opening condition for even checking iteration with fmax > 1.5*fmin
    } // opening condition for considering any molecular self-shielding terms at all
    if(!isfinite(fH2)) {fH2=0;} else {if(fH2>1) {fH2=1;} else if(fH2<0) {fH2=0;}} // check vs nans, valid values
    return xH0 * fH2; // return answer
#endif

    
#if defined(COOL_MOLECFRAC_KMT)
    /* use the simpler Kumholz, McKee, & Tumlinson 2009 sub-grid model for molecular fractions in equilibrium, which is a function modeling spherical clouds
        of internally uniform properties exposed to incident radiation. Depends on column density, metallicity, and incident FUV field. */
    /* get estimate of mass column density integrated away from this location for self-shielding */
    double surface_density_Msun_pc2_infty = 0.05 * evaluate_NH_from_GradRho(data->GradRho,data->Hsml,data->Density,data->NumNgb,1,i) * UNIT_SURFDEN_IN_CGS / 0.000208854; // approximate column density with Sobolev or Treecol methods as appropriate; converts to M_solar/pc^2
    /* 0.05 above is in testing, based on calculations by Laura Keating: represents a plausible re-scaling of the shielding length for sub-grid clumping */
    double surface_density_Msun_pc2_local = data->Density * Get_Particle_Size(i) * All.cf_a2inv * UNIT_SURFDEN_IN_CGS / 0.000208854; // this is -just- the depth through the local cell/slab. that's closer to what we want here, since G0 is -already- attenuated in the pre-processing step!
    double surface_density_Msun_pc2 = DMIN( surface_density_Msun_pc2_local, surface_density_Msun_pc2_infty);
    //double surface_density_Msun_pc2 = surface_density_Msun_pc2_local;
    /* now actually do the relevant calculation with the KMT fitting functions */
    double clumping_factor_for_unresolved_densities = 1; // Gnedin et al. add a large clumping factor to account for inability to resolve high-densities, here go with what is resolved
    double chi = 0.766 * (1. + 3.1*pow(Z_Zsol, 0.365)); // KMT estimate of chi, useful if we do -not- know anything about actual radiation field
    if(urad_G0 >= 0) {chi = 71. * urad_G0 / (clumping_factor_for_unresolved_densities * nH_cgs);} // their actual fiducial value including radiation information
    double psi = chi * (1.+0.4*chi)/(1.+1.08731*chi); // slightly-transformed chi variable
    double s = (Z_Zsol + 1.e-3) * surface_density_Msun_pc2 / (MIN_REAL_NUMBER + psi); // key variable controlling shielding in the KMT approximaton
    double q = s * (125. + s) / (11. * (96. + s)); // convert to more useful form from their Eq. 37
    double fH2 = 1. - pow(1.+q*q*q , -1./3.); // full KMT expression [unlike log-approximation, this extrapolates physically at low-q]
    if(q<0.2) {fH2 = q*q*q * (1. - 2.*q*q*q/3.)/3.;} // catch low-q limit more accurately [prevent roundoff error problems]
    if(q>10.) {fH2 = 1. - 1./q;} // catch high-q limit more accurately [prevent roundoff error problems]
    fH2 = DMIN(1,DMAX(0, fH2)); // multiple by neutral fraction, as this is ultimately the fraction of the -neutral- gas in H2
    return xH0 * fH2;
#endif
    

#if defined(COOL_MOLECFRAC_GD)
    /* use the sub-grid final expression calibrated to ~60pc resolution simulations with equilibrium molecular chemistry and post-processing radiative
        transfer from Gnedin & Draine 2014 (Eqs. 5-7) */
    double S_slab = Get_Particle_Size(i) * All.cf_atime * UNIT_LENGTH_IN_PC / 100.; // slab size in units of 100 pc
    double D_star = 0.17 * (2. + S_slab*S_slab*S_slab*S_slab*S_slab) / (1. + S_slab*S_slab*S_slab*S_slab*S_slab); // intermediate variable
    double U_star = 9. * D_star / S_slab, n_star = 14. * sqrt(D_star) / S_slab; // intermediate variables
    double g_eff = sqrt(D_star*D_star + Z_Zsol*Z_Zsol); // intermediate variable parameterizing the dust-to-gas ratio here [assuming the dust-to-gas ratio relative to solar scales linearly with metallicity, giving Z_Zsol = their D_MW parameter]
    double Lambda_incident = log(1. + pow(0.05/g_eff + urad_G0, 2./3.) * pow(g_eff, 1./3.) / U_star); // intermediate variable parameterizing the incident radiation, takes input UV radiation field relative to MW
    double nHalf = n_star * Lambda_incident / g_eff; // intermediate variable
    double w_x = 0.8 + sqrt(Lambda_incident) / pow(S_slab, 1./3.); // intermediate variable
    double x_f = w_x * log(nH_cgs / nHalf); // intermediate variable
    fH2 = 1./(1. + exp(-x_f*(1.-0.02*x_f+0.001*x_f*x_f)));
    return xH0 * fH2;
#endif
    
    
#if (SINGLE_STAR_SINK_FORMATION & 256) || defined(GALSF_SFR_MOLECULAR_CRITERION) || defined(COOL_MOLECFRAC_KG) /* estimate f_H2 with Krumholz & Gnedin 2010 fitting function, assuming simple scalings of radiation field, clumping, and other factors with basic gas properties so function only of surface density and metallicity, truncated at low values (or else it gives non-sensical answers) */
    double clumping_factor=1, fH2_kg=0, tau_fmol = (0.1 + data->Metallicity[0]/All.SolarAbundances[0]) * evaluate_NH_from_GradRho(data->GradRho,data->Hsml,data->Density,data->NumNgb,1,i) * 434.78 * UNIT_SURFDEN_IN_CGS; // convert units for surface density. also limit to Z>=0.1, where their fits were actually good, or else get unphysically low molecular fractions
    if(tau_fmol>0) {double y = 0.756 * (1 + 3.1*pow(data->Metallicity[0]/All.SolarAbundances[0],0.365)) / clumping_factor; // this assumes all the equilibrium scalings of radiation field, density, SFR, etc, to get a trivial expression
        y = log(1 + 0.6*y + 0.01*y*y) / (0.6*tau_fmol); y = 1 - 0.75*y/(1 + 0.25*y); fH2_kg=DMIN(1,DMAX(0,y));}
    return fH2_kg * neutral_fraction;
#endif
    
    
#if defined(COOLING) || defined(COOL_MOLECFRAC_GC) /* if none of the above is set, default to a wildly-oversimplified scaling set by fits to the temperature below which gas at a given density becomes molecular from cloud simulations in Glover+Clark 2012 */
    double T_mol = DMAX(1.,DMIN(8000., data->Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS));
    return neutral_fraction / (1. + temperature*temperature/(T_mol*T_mol));
#endif
    
    return 0; // catch //
}


/* return helium -number- fraction, not mass fraction */
double yhelium(struct species_data* data)
{
#ifdef COOL_METAL_LINES_BY_SPECIES
    {double ytmp=DMIN(0.5,data->Metallicity[1]); return 0.25*ytmp/(1.-ytmp);}
#else
    return ((1.-HYDROGEN_MASSFRAC)/(4.*HYDROGEN_MASSFRAC)); // assume uniform H-He gas
#endif
}


/* return mean molecular weight, appropriate for the approximations of the user-selected chemical network[s] */
double Get_Gas_Mean_Molecular_Weight_mu(double T_guess, double rho, double *xH0, double *ne_guess, double urad_from_uvb_in_G0, struct species_data* data)
{
#if defined(COOLING)
    double X=HYDROGEN_MASSFRAC, Y=1.-X, Z=0, fmol;
#ifdef METALS
    {
        Z = DMIN(0.25,data->Metallicity[0]); if(NUM_METAL_SPECIES>=10) {Y = DMIN(0.35,data->Metallicity[1]);}
        X = 1. - (Y+Z);
    }
#endif
    fmol = Get_Gas_Molecular_Mass_Fraction(data, T_guess, *xH0, *ne_guess, urad_from_uvb_in_G0); /* use our simple subroutine to estimate this, ignoring UVB and with clumping factor=1 */
    return 1. / ( X*(1-0.5*fmol) + Y/4. + *ne_guess*HYDROGEN_MASSFRAC + Z/(16.+12.*fmol) ); // since our ne is defined in some routines with He, should multiply by universal
#else
    return 4./(3.+5.*HYDROGEN_MASSFRAC); // fully-ionized H-He plasma
#endif
}


