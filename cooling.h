#ifndef COOLING_H_
#define COOLING_H_
/*
 * This file contains the definitions for the cooling.c routines
 *
 * This file was originally part of the GADGET3 code developed by
 *   Volker Springel. The code has been modified by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#include <gsl/gsl_rng.h>

extern struct species_data
{
    double Ne; // electron number per hydrogen nucleon
    double DtInternalEnergy;
    double Density;
    double Metallicity[NUM_METAL_SPECIES]; /*!< metallicity (species-by-species) of gas or star particle */
#ifdef OUTPUT_COOLRATE_DETAIL
    double CoolingRate;
    double MetalCoolingRate;
    double HeatingRate;
    double NetHeatingRateQ;
    double PhotoheatingRate;
    double PhotoheatingEfficiency;
    double HydroHeatingRate;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 1)
    gr_float grHI;
    gr_float grHII;
    gr_float grHM;
    gr_float grHeI;
    gr_float grHeII;
    gr_float grHeIII;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 2)
    gr_float grH2I;
    gr_float grH2II;
#endif
#if (COOL_GRACKLE_CHEMISTRY >= 3)
    gr_float grDI;
    gr_float grDII;
    gr_float grHDI;
#endif
} grackle_species;

#ifdef COOL_GRACKLE
void InitGrackle(void);
double CallGrackle(double u_old, double rho, double dt, double ne_guess, struct species_data *data, int mode);
#endif

void   InitCool(void);
void   InitCoolMemory(void);
void   IonizeParams(void);
void   IonizeParamsFunction(void);
void   IonizeParamsTable(void);
void   MakeCoolingTable(void);
void   ReadIonizeParams(char *fname);
void   SetZeroIonization(void);
void   TestCool(void);

double find_abundances_and_rates(double logT, double rho, struct species_data *data, double shieldfac, int return_cooling_mode,
                                 double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess, double *mu_guess);
double convert_u_to_temp(double u, double rho, struct species_data *data, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess, double *mu_guess);
double CoolingRate(double logT, double rho, double nelec, struct species_data *data);
double CoolingRateFromU(double u, double rho, double ne_guess, struct species_data *data);
double DoCooling(double u_old, double rho, double dt, double ne_guess, struct species_data *data);
double GetCoolingTime(double u_old, double rho,  double ne_guess, struct species_data *data);
double ThermalProperties(double u, double rho, struct species_data *data, double *mu_guess, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess);


#endif // COOLING_H_