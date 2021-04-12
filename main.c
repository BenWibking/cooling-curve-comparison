#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "allvars.h"
#include "cooling.h"

#ifdef COOL_GRACKLE
#include <grackle.h>
#endif // COOL_GRACKLE

#define mh 1.67262171e-24
#define kboltz 1.3806504e-16

struct species_data grackle_species;

int main(int argc, char *argv[]) {
  // set output settings
  const int verbose = 0;

  // initialize abundances (copied from GIZMO init.c)
#ifdef METALS
  for (int j = 0; j < NUM_METAL_SPECIES; j++) {
    All.SolarAbundances[j] = 0;
  } // initialize all to zero
  All.SolarAbundances[0] =
      0.02; // all metals (by mass); present photospheric abundances from
            // Asplund et al. 2009 (Z=0.0134, proto-solar=0.0142) in notes;
            //   also Anders+Grevesse 1989 (older, but hugely-cited compilation;
            //   their Z=0.0201, proto-solar=0.0213)
#ifdef COOL_METAL_LINES_BY_SPECIES
  All.SolarAbundances[1] =
      0.28; // He  (10.93 in units where log[H]=12, so photospheric mass
            // fraction -> Y=0.2485 [Hydrogen X=0.7381]; Anders+Grevesse
            // Y=0.2485, X=0.7314), with proto-solar Y=0.27
  All.SolarAbundances[2] = 3.26e-3; // C   (8.43 -> 2.38e-3, AG=3.18e-3);
                                    // proto-solar from Asplund=8.47 -> 2.53e-3
  All.SolarAbundances[3] =
      1.32e-3; // N   (7.83 -> 0.70e-3, AG=1.15e-3); PS=7.87->7.41e-4
  All.SolarAbundances[4] =
      8.65e-3; // O   (8.69 -> 5.79e-3, AG=9.97e-3); PS=8.73->6.13e-3
  All.SolarAbundances[5] =
      2.22e-3; // Ne  (7.93 -> 1.26e-3, AG=1.72e-3); PS=7.97->1.34e-3
  All.SolarAbundances[6] =
      9.31e-4; // Mg  (7.60 -> 7.14e-4, AG=6.75e-4); PS=7.64->7.57e-4
  All.SolarAbundances[7] =
      1.08e-3; // Si  (7.51 -> 6.71e-4, AG=7.30e-4); PS=7.55->7.12e-4
  All.SolarAbundances[8] =
      6.44e-4; // S   (7.12 -> 3.12e-4, AG=3.80e-4); PS=7.16->3.31e-4
  All.SolarAbundances[9] =
      1.01e-4; // Ca  (6.34 -> 0.65e-4, AG=0.67e-4); PS=6.38->6.87e-5
  All.SolarAbundances[10] =
      1.73e-3; // Fe (7.50 -> 1.31e-3, AG=1.92e-3); PS=7.54->1.38e-3
#endif
#endif

  // initialize units
  All.UnitMass_in_g = 1.67e-24;
  All.UnitLength_in_cm = 1.0;
  double time_units = 1.0e12;
  All.UnitVelocity_in_cm_per_s = All.UnitLength_in_cm / time_units;

  // set cosmological parameters
  All.cf_atime = 1.;    // non-cosmological
  All.HubbleParam = 1.; // non-cosmological

  // set temperature units (approximate)
  double temperature_units =
      mh * pow(All.UnitLength_in_cm / time_units, 2) / kboltz;

  const int ngas = 600;           // bins of density
  const double logrho_min = -6.0; // 1e-6 H cm^-3
  const double logrho_max = 4.0;  // 1e4 H cm^-3
  double *rho_bins = malloc(sizeof(double) * (ngas + 1));

  for (int i = 0; i <= ngas; i++) {
    const double logrho =
        (logrho_max - logrho_min) * ((double)i / (double)ngas) + logrho_min;
    rho_bins[i] = pow(10., logrho);
  }

  fprintf(stdout, "# [H cm^-3] [K] [ergs cm^-3] [s] [dimensionless]\n");
  fprintf(stdout, "rho temperature pressure cooling_time HI_fraction\n");

#ifdef COOL_GRACKLE
  // set photoionizing background
  strcpy(All.GrackleDataFile, "../grackle/input/CloudyData_UVB=HM2012.h5");
  // initialize Grackle
  InitGrackle();
#else
  // initialize GIZMO cooling
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, 42); /* start-up seed */
  InitCool();
#endif

  for (int i = 0; i <= ngas; i++) {
    // set gas properties
    const double rho = rho_bins[i];
    const double T0 = 1.0e6; // K [starting temperature is arbitrary, unless
                             // there are multiple equilibria]
#ifdef COOL_GRACKLE
    const double dt_initial = 1.0e2; // yr
#else
    const double dt_initial = 1.0e6; // yr
#endif

    if (verbose) {
      fprintf(stderr, "rho = %g H cm^-3\n", rho);
    }

    double eint = (T0 / temperature_units);
    double dt = 3.15e7 * dt_initial / time_units;

    // initialize species abundances
    double ne_guess = 1.0e-2;

#ifdef COOL_GRACKLE
    double tiny_number = 1.0e-20;
    grackle_species.grHI = HYDROGEN_MASSFRAC;
    grackle_species.grHII = tiny_number;
    grackle_species.grHM = tiny_number;

    grackle_species.grHeI = (1.0 - HYDROGEN_MASSFRAC);
    grackle_species.grHeII = tiny_number;
    grackle_species.grHeIII = tiny_number;
#else
    grackle_species.Density = rho;
    grackle_species.Ne = ne_guess;
    grackle_species.Metallicity[0] = 0.02; // assume solar
    // copied from init.c:
    /* initialize abundance ratios. for now, assume [scaled] solar */
    for (int j = 0; j < NUM_METAL_SPECIES; j++) {
      grackle_species.Metallicity[j] =
          All.SolarAbundances[j] *
          (grackle_species.Metallicity[0] / All.SolarAbundances[0]);
    }
    /* need to allow for a primordial He abundance */
    if (NUM_LIVE_SPECIES_FOR_COOLTABLES >= 10)
      grackle_species.Metallicity[1] =
          (1. - HYDROGEN_MASSFRAC) +
          (All.SolarAbundances[1] - (1. - HYDROGEN_MASSFRAC)) *
              grackle_species.Metallicity[0] / All.SolarAbundances[0];
#endif

    // iterate until temperature converges
    const int maxiter = 1e5;
#ifdef COOL_GRACKLE
    const double HI_rel_tol = 1.0e-4;
    const double temp_rel_tol = 1.0e-4;
#else
    const double HI_rel_tol = 1.0e-3;
    const double temp_rel_tol = 1.0e-3;
#endif
    double safety_fac = 0.1; // relative to the cooling time
    int is_converged = 0;
    double T_prev = T0;
    double HI_prev = HYDROGEN_MASSFRAC;
    double HI_converged = NAN;
    double T_converged = NAN;
    double P_converged = NAN;
    double tcool_converged = NAN;

    for (int j = 0; j < maxiter; j++) {
#ifdef COOL_GRACKLE
      // 0 == compute equilibrium ionization state
      eint = CallGrackle(eint, rho, dt, ne_guess, &grackle_species, 0);
      const double HI_fraction = grackle_species.grHI / HYDROGEN_MASSFRAC;
      if (verbose) {
        fprintf(stderr, "Ionization state:\n");
        fprintf(stderr, "\tHI = %g\n", grackle_species.grHI);
        fprintf(stderr, "\tHII = %g\n", grackle_species.grHII);
        fprintf(stderr, "\tHeI = %g\n", grackle_species.grHeI);
        fprintf(stderr, "\tHeII = %g\n", grackle_species.grHeII);
        fprintf(stderr, "\tHeIII = %g\n", grackle_species.grHeIII);
      }

      // 1 == calculate cooling time
      const double cooling_time =
          CallGrackle(eint, rho, dt, ne_guess, &grackle_species, 1);

      // 2 == calculate temperature
      const double temperature =
          CallGrackle(eint, rho, dt, ne_guess, &grackle_species, 2);

      // 3 == calculate pressure
      const double pressure =
          CallGrackle(eint, rho, dt, ne_guess, &grackle_species, 3);
      const double pressure_cgs = pressure * UNIT_PRESSURE_IN_CGS;
#else // GIZMO cooling
      // set global variables
      All.Timebase_interval = dt;

      // compute cooling
      eint = DoCooling(eint, rho, dt, ne_guess, &grackle_species);
      if (verbose) {
        double nHcgs = HYDROGEN_MASSFRAC * rho;
        double photoheating_per_H = nHcgs * grackle_species.PhotoheatingRate;
        fprintf(stderr, "Photoheating rate (per H) = %g ergs/s\n",
                photoheating_per_H);
        fprintf(stderr, "Photoheating efficiency = %g\n",
                grackle_species.PhotoheatingEfficiency);
        fprintf(stderr, "Electron number per H = %g\n", ne_guess);
        fprintf(stderr, "Heating rate (per H) = %g\n",
                nHcgs * grackle_species.HeatingRate);
        fprintf(stderr, "Cooling rate (per H) = %g\n",
                nHcgs * grackle_species.CoolingRate);
      }

      // compute temperature and HI fraction
      const double rho_cgs = rho * UNIT_DENSITY_IN_CGS;
      const double eint_cgs = eint * UNIT_SPECEGY_IN_CGS;
      double nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu;
      const double temperature = convert_u_to_temp(
          eint_cgs, rho_cgs, &grackle_species, &ne_guess, &nH0_guess,
          &nHp_guess, &nHe0_guess, &nHep_guess, &nHepp_guess, &mu);
      const double HI_fraction = nH0_guess;

      // compute cooling time
      const double cooling_time =
          GetCoolingTime(eint, rho, ne_guess, &grackle_species);

      // compute pressure (assuming constant EOS_GAMMA)
      const double pressure = (EOS_GAMMA - 1.) * eint * rho;
      const double pressure_cgs = pressure * UNIT_PRESSURE_IN_CGS;
#endif
      if (verbose) {
        fprintf(stderr, "Cooling time = %g s.\n", cooling_time * time_units);
        fprintf(stderr, "Temperature = %g K.\n", temperature);
        fprintf(stderr, "P/k_B = %g.\n\n", pressure_cgs / kboltz);
      }

      // check if converged
      if ((fabs((T_prev - temperature) / temperature) < temp_rel_tol) &&
          (fabs((HI_prev - HI_fraction) / HI_fraction) < HI_rel_tol)) {
        is_converged = 1;
        break;
      }

      // reduce safety factor if it is taking a really long time to converge
      if ((j > 0) && ((j % 1000) == 0)) {
        safety_fac *= 0.5;
      }

      // calculate new dt based on cooling time
      if (cooling_time != 0.) {
        dt = safety_fac * fabs(cooling_time);
      }

      // save temperature, HI fraction
      T_prev = temperature;
      HI_prev = HI_fraction;

      // save variables
      HI_converged = HI_fraction;
      T_converged = temperature;
      P_converged = pressure_cgs;
      tcool_converged = cooling_time;
    }

    // output data
    fprintf(stdout, "%.10e %.10e %.10e %.10e %.10e %d\n", rho, T_converged,
            P_converged, tcool_converged * time_units, HI_converged,
            is_converged);
  }

  return 0;
}
