/**********************************************************************
UDF for defining the heat and mass transport for n-dodecane droplet
Has been tested with ANSYS Fluent 17.2

Physical properties are from
Abramzon, B. and W.A. Sirignano, Droplet Vaporization Model for Spray Combustion Calculations. International Journal of Heat and Mass Transfer, 1989. 32(9): p. 1605-1618.


If you use our software, please cite the following:
1. O. Rybdylova, M. Al Qubeissi, M. Braun, C. Crua, J. Manin, L.M. Pickett, G. de Sercey, E.M. Sazhina, S.S. Sazhin, M. Heikal, A model for droplet heating and its implementation into ANSYS Fluent, International Communications in Heat and Mass Transfer, Volume 76, 2016, Pages 265-270, ISSN 0735-1933,
https://doi.org/10.1016/j.icheatmasstransfer.2016.05.032.

2. Timur S. Zaripov, Oyuna Rybdylova, Sergei S. Sazhin, A model for heating and evaporation of a droplet cloud and its implementation into ANSYS Fluent, International Communications in Heat and Mass Transfer, Volume 97, 2018, Pages 85-91, ISSN 0735-1933,
https://doi.org/10.1016/j.icheatmasstransfer.2018.06.007.

Copyright (C) 2018 Oyuna Rybdylova, Timur Zaripov - All Rights Reserved
You may use, distribute and modify this code under the terms of the MIT license

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
***********************************************************************/
#include "udf.h"
#include <time.h>

#define BM_MAX 1.E20
#define BM_MIN -0.99999
#define ACCURACY 1.e-6
#define PI 3.1415926535897932384626433832795
#define N_Lambda 44 // number of terms in the series
#define N_INT 100 // number of layers inside a droplet
#define Delta_R 0.01 // = 1/N_INT

// 116 DPM_USER_REALs have to be enabled in ANSYS Fluent if 100 layers inside a droplet is used
#define NCOMPONENTS 1 // Number of components, there is a check in Heat and Mass transfer on the number of components
#define VAP_END (116)

// BEGIN VAP functions 
int Lambda(real h_0, real lambda[])
{
    FILE * fout;
    int i;
    double lambda_left, lambda_right, f_left, f_right, lambda_mid, f_mid;
    double conv_crit = 1.e-8;
    double step = 1.e-7;

    for (i = 0; i < N_Lambda; i++) lambda[i] = -1.0;

    for (i = 0; i < N_Lambda; i++)
    {
        lambda_left = ((double)(i))*PI + step;
        lambda_right = (((double)(i + 1)) - 0.5)*PI - step;

        if (h_0 > 0.0)
        {
            lambda_left += 0.5*PI;
            lambda_right += PI*0.5;
        }

        f_left = lambda_left*cos(lambda_left) + h_0*sin(lambda_left);
        f_right = lambda_right*cos(lambda_right) + h_0*sin(lambda_right);
        if (f_left*f_right < 0.0)
        {
            while (lambda_right - lambda_left > conv_crit)
            {
                lambda_mid = (lambda_left + lambda_right)*0.5;
                f_mid = lambda_mid*cos(lambda_mid) + h_0*sin(lambda_mid);
                if (f_left*f_mid < 0.0)
                {
                    lambda_right = lambda_mid;
                    f_right = lambda_right*cos(lambda_right) + h_0*sin(lambda_right);
                }
                else
                {
                    lambda_left = lambda_mid;
                    f_left = lambda_left*cos(lambda_left) + h_0*sin(lambda_left);
                }
                //Message("Lambda left = %f\nLambda right = %f\n f left = %e\nf right = %e\n------------\n", lambda_left, lambda_right, f_left, f_right);
            }
            lambda[i] = lambda_left;
        }
    }

    //  Message("Lambdas are calculated\n");

    //  fout = fopen("lambda.txt", "w");
    //  for (i = 0; i < N_Lambda; i++)
    //      fprintf(fout, "%d\t%20.19f\t%e\n", i + 1, lambda[i], lambda[i] * cos(lambda[i]) + h_0*sin(lambda[i]));
    //  fclose(fout);
    //  Message("Lambdas are printed in lambda.txt\n");
    return 0;
}
// END VAP functions

/* convection diffusion controlled vaporisation model as implemented into Fluent
   p    ... tracked particle struct
   Cp   ... particle heat capacity
   hgas ... enthalpy of formation for gas species
   hvap ... vaporization enthalpy
   cvap_surf ... molar concentration of vapor at surface
   Z    ... compressibility factor
   dydt ... temperature / component mass source terms for particle 
            dydt[0] = temperature equation, 
            dydt[1+ns] = mass equation of particle species ns
   dzdt ... temperature / species mass source terms for gas phase
            dzdt->energy = convective heat transfer in Watt
            dzdt->species[n_gas] = mass source for gas species n_gas in kg/s
 */

// 4*N_component + 7 (x_i, y_i, dm_i, Mw_i, y_tot, dm_tot, D, BM, BT, diam1, diam2) + N_INT+1 
// for temperature distribution inside a droplet USER_REAL variables 
// 116 for single component n-dodecane
DEFINE_DPM_HEAT_MASS(multivap_conv_diffusion_new, p, Cp, hgas, hvap, cvap_surf, Z, dydt, dzdt)
{
    //-------------------------------------------------------------------------
    /* molecular weight of gas species */
    /* molwt_bulk is first the reciproke of the real molecular weight to avoid additional divisions */
    real molwt_bulk = 0.;
    int ns;
    Thread *t0 = P_CELL_THREAD(p);
    Material *gas_mix = THREAD_MATERIAL(DPM_THREAD(t0, p));
	cphase_state_t *c = &(p->cphase);  /* continuous phase struct, caching variables of the cell */

    mixture_species_loop_i(gas_mix, ns)
    {
        molwt_bulk += c->yi[ns] / solver_par.molWeight[ns]; /* c->yi[ns] contains mass fraction of fluid cell */
        //P_USER_REAL(p, 3 * nc + ns) = solver_par.molWeight[ns];
    }
    molwt_bulk = MAX(molwt_bulk, DPM_SMALL);
    molwt_bulk = 1. / molwt_bulk;
    

    /* when not using Runge Kutta solver, increase limiting time for the next integration step */
    // ANSYS stuff
    if (!p->in_rk) {
        p->limiting_time = P_DT(p)*1.01;
    }

    //-------------------------------------------------------------------------
    // Calculate molar fractions of the components at the droplet surface
    real xs_tot = 0.e-15;
    real xsM_tot = 0.e-15;
    real P_sat = 0.0; // Saturation pressure
    real x_surf =0.0; // molar fraction of component at droplet surface

    int nc = TP_N_COMPONENTS(p);
	real Tp = P_USER_REAL(p, 4 * nc + 7 + N_INT); //Dropet temperature at the surface
	if (nc != NCOMPONENTS) {
        Message("ALARM!!! nc != NCOMPONENTS.");
    }
    for (int ns = 0; ns < nc; ns++) {
        int gas_index = TP_COMPONENT_INDEX_I(p, ns); /* gas species index of vaporization */
        if (gas_index >= 0) {
            // FIXME {{{
            // Saturation pressure for n-dodecane vapour
            P_sat = exp(8.1948 - 7.8099*(300.0 / Tp) - 9.0098*(300.0 / Tp)*(300.0 / Tp))*1.e5;//Abramzon, B. and S. Sazhin, Convective vaporization of a fuel droplet with thermal radiation absorption. Fuel, 2006. 85(1): p. 32-46.
            if (Tp > 0.99*659.0) {
                P_sat = P_sat*exp(15.0*(Tp / 0.99 / 659.0 - 1.0));
            }
            // }}}
            x_surf = P_sat / c->pressure; //Saturation pressure for n-Dodecane from Abramzon&Sazhin 2006
            //above for x_surf will be modified for multicoponent droplet case
            P_USER_REAL(p, ns) = x_surf;
            xs_tot += x_surf*solver_par.molWeight[ns];
            xsM_tot += x_surf;
        }
    }
    
    xs_tot += (1.0 - xsM_tot)*28.967; //air
    
    //-------------------------------------------------------------------------
    // Calculate mass fractions of the components at the droplet surface
    real Ys = 0.e-15;
    real Y_inf = 0.e-15;
    real L_eff = 0.e-15;
    real Ys_tot = 0.e-15;
    for (int ns = 0; ns < nc; ns++) {
        /* gas species index of vaporization */
        int gas_index = TP_COMPONENT_INDEX_I(p, ns);
        if (gas_index >= 0) {
            Ys = P_USER_REAL(p, ns)* solver_par.molWeight[gas_index] / xs_tot;//!!
            Y_inf += c->yi[gas_index];
            //L_eff += Ys * p->hvap[gas_index]; // TODO Try for water
            
            // FIXME {{{
            // Effective latent heat for n-dodecane vapour
            L_eff += Ys * 37.44*pow(659.0 - Tp, 0.38)*1000.0;//Abramzon, B. and S. Sazhin, Convective vaporization of a fuel droplet with thermal radiation absorption. Fuel, 2006. 85(1): p. 32-46.
            if (Tp>0.99*659.0) {
                L_eff = Ys * 37.44*pow(659.0 - 653.0, 0.38)*1000.0;
            }
            // }}}
            //Latent heat as above will be calculated separately for multicomponent droplet
            Ys_tot += Ys;
            P_USER_REAL(p, nc + ns) = Ys;
        }
    }
    L_eff = L_eff / Ys_tot;
    P_USER_REAL(p, 4 * nc) = Ys_tot;

    //-------------------------------------------------------------------------
    // Calculate Nusselt number and total evaporation rate
    real T_ref = (c->temp + 2.0*Tp) / 3.0; //Sazhin, Progress in Energy and Combustion Science 32 (2006) 162–214 
    real rho_gas_s = c->pressure / (287.01625988193461525183829875375*T_ref); // ideal gas law
    // FIXME {{{
    real c_p_die = (0.2979 + 1.4394*(T_ref / 300.0) - 0.1351*(T_ref / 300.0)*(T_ref / 300.0))*1000.0; //n-dodecane vapour Abramzon, B. and S. Sazhin, Convective vaporization of a fuel droplet with thermal radiation absorption. Fuel, 2006. 85(1): p. 32-46.
    //D = 0.527*pow(T_ref / 300.0, 1.583) / c->pressure; //2015.10.29
    // }}}

	Material *cond_mix = P_MATERIAL(p);
	Material * cond_c = MIXTURE_COMPONENT(cond_mix, 0);
	real D = DPM_BINARY_DIFFUSIVITY(p, cond_c, Tp);
    real Sc = c->mu / (rho_gas_s * D);  //Schmidt number

    real kgas = c->tCond;
    real Re = p->Re;
    real Pr = c->sHeat * c->mu / kgas;
    //  BM = (Ys_tot - Y_inf) / (1.0 - Ys_tot);
    real BM = (Ys_tot) / (1.0 - Ys_tot); //assuming zero mass fraction in the ambient gas
    real FBM = pow(1.0 + BM, 0.7)*log(1.0 + BM) / BM;
    real Sh_Star = 2.0 + (pow(1.0 + Re*Sc, 1.0 / 3.0)*MAX(1.0, pow(Re, 0.077)) - 1.0) / FBM;
    real Sh = log(1.0 + BM)*Sh_Star;
    //Sh = log(1.0 + BM)*(2.0 + 0.6*sqrt(Re)*pow(Sc, 1.0 / 3.0));
	real Dp = P_DIAM(p);
	real Ap = DPM_AREA(Dp);
	real tot_vap_rate = Ap * D * rho_gas_s * Sh / Dp; // total evaporation rate
    
    P_USER_REAL(p, 4 * nc + 1) = tot_vap_rate;

    real BT = BM;
    real BT_i = BT;
    real dif = 1.0;
    real coef = c_p_die * rho_gas_s * D / kgas * Sh_Star;
    real phi = 0.e-15;
	real FBT, Nu_star;
    // find BT iteratively
    while (dif > ACCURACY) {
        FBT = pow(1.0 + BT, 0.7)*log(1.0 + BT) / BT;
        Nu_star = 2.0 + (pow(1.0 + Re*Pr, 1.0 / 3.0)*MAX(1.0, pow(Re, 0.077)) - 1.0) / FBT;
        phi = coef / Nu_star;
        BT = pow(1.0 + BM, phi) - 1.0;
        dif = abs(BT - BT_i);
        BT_i = BT;
    }

    real Nu = log(1.0 + BT) * Nu_star / BT; // Nusselt number

    //-------------------------------------------------------------------------
    // Temperature distribution calculations
    real T_av = P_USER_REAL(p, 4 * nc + 6);
    // FIXME {{{
    // _l stands for "liquid"
    real Visc_l = 1.e-3*exp(2.0303*(300.0 / T_av)*(300.0 / T_av) + 1.1769*(300.0 / T_av) - 2.929);//Abramzon, B. and S. Sazhin, Convective vaporization of a fuel droplet with thermal radiation absorption. Fuel, 2006. 85(1): p. 32-46.
    real k_l = 0.1405 - 0.00022*(T_av - 300.0);//Abramzon, B. and S. Sazhin, Convective vaporization of a fuel droplet with thermal radiation absorption. Fuel, 2006. 85(1): p. 32-46.
    //k_l *= 1000.0;  // for test purposes, to model infinitely high thermal conductivity
    real C_pl = (2.18 + 0.0041*(T_av - 300.0))*1000.0;//Abramzon, B. and S. Sazhin, Convective vaporization of a fuel droplet with thermal radiation absorption. Fuel, 2006. 85(1): p. 32-46.
    //C_pl = p->Cp;
    // }}}

    real rel_vel = sqrt((c->V[0] - P_VEL(p)[0])*(c->V[0] - P_VEL(p)[0]) + (c->V[1] - P_VEL(p)[1])*(c->V[1] - P_VEL(p)[1]) + (c->V[2] - P_VEL(p)[2])*(c->V[2] - P_VEL(p)[2]));
    real Pe = 12.69 / 16.0*P_RHO(p)*0.5*Dp* C_pl / k_l*rel_vel*c->mu / Visc_l*pow(Re, 1.0 / 3.0) / (1.0 + BM);
    real k_eff = (1.86 + 0.86*tanh(2.225*log10(Pe / 30.0)))*k_l;  // effective thermal conductivity to take into account recirculation Abramzon B, Sirignano WA. Int J Heat Mass Transfer 1989;32:1605–18.
    if (abs(Pe) < 1.e-12) {
        k_eff = k_l;
    }

    real T_eff = c->temp - tot_vap_rate*L_eff / PI / Dp / Nu / kgas;
    real h0 = kgas*Nu*0.5 / k_eff - 1.0;
    real zeta = (h0 + 1.0)*T_eff;
    real kappa = k_eff / (C_pl*P_RHO(p)*0.25*Dp*Dp);

	real lambda[N_Lambda];
	for (int i = 0; i < N_Lambda; i++) { lambda[i] = -1.0; }
    int err = Lambda(h0, lambda);

	real series[N_Lambda];
    for (int i = 0; i < N_Lambda; i++) { series[i] = 0.e-15; }

	real I_n, b_n;
	for (int i = 0; i < N_Lambda; i++)  {
        b_n = 0.5*(1.0 + h0 / (h0*h0 + lambda[i] * lambda[i]));
        I_n = 0.e-15;
        I_n = P_USER_REAL(p, 4 * nc + 7 + N_INT)*sin(lambda[i]);
        for (int j = 1; j < N_INT; j += 2) {
            I_n += 4.0 * P_USER_REAL(p, 4 * nc + 7 + j)*(((double)j)*Delta_R)*sin(lambda[i] * ((double)j)*Delta_R);
        }
        for (int j = 2; j < N_INT; j += 2) {
            I_n += 2.0 * P_USER_REAL(p, 4 * nc + 7 + j)*(((double)j)*Delta_R)*sin(lambda[i] * ((double)j)*Delta_R);
        }
        I_n = I_n*Delta_R / 3.0;
        series[i] = (I_n - sin(lambda[i]) / lambda[i] / lambda[i] * zeta)*exp(0.0 - kappa*lambda[i] * lambda[i] * P_DT(p)) / b_n;
    }

    for (int j = 0; j < N_INT + 1; j++) { P_USER_REAL(p, 4 * nc + 7 + j) = T_eff; }
    for (int i = 0; i < N_Lambda; i++)  {
        P_USER_REAL(p, 4 * nc + 7) += series[i] * lambda[i];
        for (int j = 1; j < N_INT + 1; j++) P_USER_REAL(p, 4 * nc + 7 + j) += series[i] * sin(lambda[i] * ((double)j)*Delta_R) / (((double)j)*Delta_R);
    }
    // Now we know temperature at each layer

    // Re-calculate droplet avarage temperature T_av
    Tp = P_USER_REAL(p, 4 * nc + 7 + N_INT);
    T_av = Tp;
    for (int j = 1; j < N_INT; j += 2) {
        T_av += 4.0 * P_USER_REAL(p, 4 * nc + 7 + j)*(((double)j)*Delta_R)*(((double)j)*Delta_R);
    }
    for (int j = 2; j < N_INT; j += 2) {
        T_av += 2.0 * P_USER_REAL(p, 4 * nc + 7 + j)*(((double)j)*Delta_R)*(((double)j)*Delta_R);
    }
    T_av = T_av*Delta_R;


    //-------------------------------------------------------------------------
    // update Fluent variables using our values
    p->state.temp = T_av;
    //p->source.htc = Nu*kgas*P_DIAM(p);
    p->source.htc = 0.e-15;  // htc - heat transfer coefficient
    
    // evaporation rates - source terms, droplet mass
    for (int ns = 0; ns < nc; ns++) {
        /* gas species index of vaporization */
        int gas_index = TP_COMPONENT_INDEX_I(p, ns);
        if (gas_index >= 0) {
            real vap_rate = P_USER_REAL(p, nc + ns) * tot_vap_rate / Ys_tot; //!!

            // ANSYS stuff
            if ((!p->in_rk) && (ABS(vap_rate)>0.)) {
                p->limiting_time = MIN(p->limiting_time, dpm_par.fractional_change_factor_mass*P_MASS(p) / vap_rate*TP_COMPONENT_I(p, ns));
            }

            P_USER_REAL(p, 2 * nc + ns) = vap_rate;
            dydt[1 + ns] -= vap_rate;

            {
                int source_index = injection_par.yi2s[gas_index];
                if (source_index >= 0) {
                    dzdt->species[source_index] += vap_rate;
                    //p->source.mtc[source_index] = c->rho * Ap * Sh_Star * D / Dp;
                    p->source.mtc[source_index] = c->rho * PI * Dp * Sh_Star * D;
                }
            }
        }
    }
    
    // Keep particle temperature independent form Fluent's source term.
    // source terms for energy equations: zero, as Temperature is calculated explicitly
    dydt[0] = 0.e-15;


    real dh_dt = Nu * kgas * Ap / Dp * (c->temp - T_av);
    dzdt->energy -= dh_dt;

    //-------------------------------------------------------------------------
    // ANSYS stuff
    real h = Nu * kgas / Dp;
	real mp = P_MASS(p);
    real convective_heating_rate = h * Ap / (mp * p->Cp);

    /* limit for higher heating rate */

    if ((!p->in_rk) && (ABS(convective_heating_rate)>DPM_SMALL)) {
            real factor = dpm_par.fractional_change_factor_heat;
            if (ABS(c->temp - Tp)>Tp) {
                factor = dpm_par.fractional_change_factor_heat*Tp / (c->temp - Tp);
            }
            p->limiting_time = MIN(p->limiting_time, factor / ABS(convective_heating_rate));
    }
    
    //-------------------------------------------------------------------------
    // update user reals
    P_USER_REAL(p, 4 * nc + 2) = BM;
    P_USER_REAL(p, 4 * nc + 3) = BT;
    P_USER_REAL(p, 4 * nc + 4) = L_eff; // used in temperature calculations
    //P_USER_REAL(p, 4 * nc + 5) = P_DIAM(p);
    //P_USER_REAL(p, 4 * nc + 6) = Dp;
    P_USER_REAL(p, 4 * nc + 5) = Nu; //used in temperature calculations
    P_USER_REAL(p, 4 * nc + 6) = T_av;

    P_USER_REAL(p, 4 * nc + 7 + N_INT + 1) = coef;
    P_USER_REAL(p, 4 * nc + 7 + N_INT + 2) = Nu_star;
    P_USER_REAL(p, 4 * nc + 7 + N_INT + 3) = D;
    P_USER_REAL(p, 4 * nc + 7 + N_INT + 4) = kgas;
    P_USER_REAL(p, 3 * nc + 0) = h;
    
}

DEFINE_DPM_SCALAR_UPDATE(Diesel_droplet, cell, thread, initialize, p)
{
    int nc = TP_N_COMPONENTS(p);
    if (nc != NCOMPONENTS) {
        Message("ALARM!!! nc != NCOMPONENTS.");
    }
    real Tp = P_T(p);
    cphase_state_t *c = &(p->cphase);
    clock_t t;
    t = clock();
    // Message("I'm in scalar update\n");
    // FIXME Only works in steady state case
    if (initialize) {
        for (int i = 0; i < N_INT + 1; i++) { P_USER_REAL(p, 4 * nc + 7 + i) = Tp; }
        P_USER_REAL(p, 4 * nc + 6) = Tp;
        P_USER_REAL(p, 4 * nc + 5) = 2.0;
        P_USER_REAL(p, 4 * nc + 4) = p->hvap[0];
        P_USER_REAL(p, 4 * nc + 2) = 0.0; //!!!
        P_USER_REAL(p, 4 * nc + 3) = 0.0; //!!!
        P_USER_REAL(p, 4 * nc + 1) = 0.e-15;
        P_USER_REAL(p, 4 * nc + 7 + N_INT + 1) = P_DIAM(p);
        P_USER_REAL(p, 4 * nc + 7 + N_INT + 2) = DPM_DIAM_FROM_VOL(P_MASS(p) / P_RHO(p));
        // P_USER_REAL(p, 4 * nc + 7 + N_INT + 2) = 8.74856E-06;
        // P_USER_REAL(p, 4 * nc + 7 + N_INT + 4) = 0.0;
        //P_USER_REAL(p, 4 * nc + 7 + N_INT + 3) = ((real) t)/CLOCKS_PER_SEC;
        //P_USER_REAL(p, 4 * nc + 7 + N_INT + 4) = 0.0;
        // Message("Temperature initiated\n");
    } else {
        //P_USER_REAL(p, 4 * nc + 7 + N_INT + 4) = ((real)t) / CLOCKS_PER_SEC - P_USER_REAL(p, 4 * nc + 7 + N_INT + 3);
        //
        // IMPORTANT for heating and evaporation
        //
        p->state.temp = P_USER_REAL(p, 4 * nc + 6);
    }
}

DEFINE_DPM_TIMESTEP(Constant_dt, p, dt)
{
    return 1.e-6;
}

// BEGIN n-dodecane properties
DEFINE_DPM_PROPERTY(Diesel_liquid_density, c, t, p, T)
{
    return 744.11 - 0.771*(P_T(p) - 300.0);
    //or set constant density
}

DEFINE_DPM_PROPERTY(Diesel_liquid_specific_heat, c, t, p, T)
{
    return (2.18 + 0.0041*(P_T(p) - 300.0))*1000.0;
}

DEFINE_DPM_PROPERTY(Diesel_latent_heat, c, t, p, T)
{
    real L;

    //L = 37.44*pow(659.0 - P_T(p), 0.38)*1000.0;
    L = 37.44*pow(659.0 - T, 0.38)*1000.0;
    if (P_T(p) >0.99*659.0) L = 37.44*pow(659.0 - 653.0, 0.38)*1000.0;
    return L;
}

DEFINE_DPM_PROPERTY(Diesel_binary_diffusivity, c, t, p, T)
{
    real D, T_ref;
    cphase_state_t *carr = &(p->cphase);  /* continuous phase struct, caching variables of the cell */
    int nc = TP_N_COMPONENTS(p);

    if (P_USER_REAL(p, 4 * nc + 7 + N_INT) < P_T(p)) T_ref = (2.0*P_T(p) + carr->temp) / 3.0;
    else T_ref = (2.0*P_USER_REAL(p, 4 * nc + 7 + N_INT) + carr->temp) / 3.0;
    D = 0.527*pow(T_ref / 300.0, 1.583) / carr->pressure;
    return D;
}

DEFINE_DPM_PROPERTY(Diesel_saturation_vapour_pressure, c, t, p, T)
{
    real P_sat;

    P_sat = exp(8.1948 - 7.8099*(300.0 / P_T(p)) - 9.0098*(300.0 / P_T(p))*(300.0 / P_T(p)))*1.e5;
    if (P_T(p) >0.99*659.0) P_sat = P_sat*exp(15.0*(P_T(p) / 0.99 / 659.0 - 1.0));
    return P_sat;
}
// END n-dodecane properties

// same as before, only k_l = k_l*1000
DEFINE_DPM_HEAT_MASS(multivap_conv_diffusion_kl, p, Cp, hgas, hvap, cvap_surf, Z, dydt, dzdt)
{
    int ns;
    real Bt = 0.;
    real tot_vap_rate = 0.;
    real vap_rate;
    int nc = TP_N_COMPONENTS(p);
    Thread *t0 = P_CELL_THREAD(p);
    Material *gas_mix = THREAD_MATERIAL(DPM_THREAD(t0, p));
    Material *cond_mix = P_MATERIAL(p);
    cphase_state_t *c = &(p->cphase);  /* continuous phase struct, caching variables of the cell */
    real Tp;
    real mp = P_MASS(p);
    /* molecular weight in bulk gas */
    real molwt_bulk = 0.;

    //real Dp = DPM_DIAM_FROM_VOL(mp / P_RHO(p));
    real Dp = P_DIAM(p);
    real Ap = DPM_AREA(Dp);
    real Re, Pr, Nu, Sh, Sc, Nu_star, Sh_Star;
    real Ys, Y_inf, Ys_tot, rho_gas_s, xs_tot, xsM_tot;
    real BM, BT, BT_i, FBM, FBT;
    real x_surf, P_sat, Visc_l, k_l, C_pl;
    real kgas;
    real T_ref;
    int i, ret, j;
    real c_p_die;
    Material * cond_c = MIXTURE_COMPONENT(cond_mix, 0);
    real D;
    real L_eff;
    real convective_heating_rate = 0.;
    real T_av;
    real rel_vel, Pe, k_eff;
    real h0, zeta, kappa;
    real lambda[N_Lambda];
    real series[N_Lambda];
    real I_n, b_n;
    real T_eff, dif, coef, phi, dh_dt, h, factor;



    Tp = P_USER_REAL(p, 4 * nc + 7 + N_INT);
    T_ref = (c->temp + 2.0*Tp) / 3.0;

    //  c_p_die = c->sHeat;

    //  kgas = 0.02667*((c->temp + 2.0*Tp) / 900.0) - 0.02087;
    kgas = c->tCond;
    //  real molwt_vap, molwt_vap;

    Tp = P_USER_REAL(p, 4 * nc + 7 + N_INT);

    rho_gas_s = c->pressure / (287.01625988193461525183829875375*T_ref);


    Re = p->Re;
    Pr = c->sHeat * c->mu / kgas;
    /* Heat transfer coefficient */

    /* molecular weight of gas species */
    /* molwt_bulk is first the reciproke of the real molecular weight to avoid additional divisions */
    mixture_species_loop_i(gas_mix, ns)
    {
        molwt_bulk += c->yi[ns] / solver_par.molWeight[ns]; /* c->yi[ns] contains mass fraction of fluid cell */
        //P_USER_REAL(p, 3 * nc + ns) = solver_par.molWeight[ns];
    }
    molwt_bulk = MAX(molwt_bulk, DPM_SMALL);
    molwt_bulk = 1. / molwt_bulk;

    /* when not using Runge Kutta solver, increase limiting time for the next integration step */
    if (!p->in_rk)
        p->limiting_time = P_DT(p)*1.01;

    T_av = Tp;
    for (j = 1; j < N_INT; j += 2) T_av += 4.0 * P_USER_REAL(p, 4 * nc + 7 + j)*(((double)j)*Delta_R)*(((double)j)*Delta_R);
    for (j = 2; j < N_INT; j += 2) T_av += 2.0 * P_USER_REAL(p, 4 * nc + 7 + j)*(((double)j)*Delta_R)*(((double)j)*Delta_R);
    T_av = T_av*Delta_R;



    /* loop over all components of multi component particle */
    Ys = 0.e-15;
    Y_inf = 0.e-15;
    L_eff = 0.e-15;
    Ys_tot = 0.e-15;
    xs_tot = 0.e-15;
    xsM_tot = 0.e-15;
    for (ns = 0; ns < nc; ns++)
    {
        /* gas species index of vaporization */
        int gas_index = TP_COMPONENT_INDEX_I(p, ns);
        if (gas_index >= 0)
        {
            //P_sat = exp(8.1948 - 7.8099*(300.0 / Tp) - 9.0098*(300.0 / Tp)*(300.0 / Tp))*1.e5;
            //if (Tp>0.99*659.0) P_sat = P_sat*exp(15.0*(Tp / 0.99 / 659.0 - 1.0));
            P_sat = DPM_VAPOR_PRESSURE(p, cond_c, Tp);
            //real x_surf = cvap_surf[ns] * Z * UNIVERSAL_GAS_CONSTANT * Tp / c->pressure;
            x_surf = P_sat / c->pressure; //Saturation pressure for n-Dodecane from Abramzon&Sazhin 2006
            P_USER_REAL(p, ns) = x_surf;
            xs_tot += x_surf*solver_par.molWeight[ns];
            xsM_tot += x_surf;
        }
    }
    xs_tot += (1.0 - xsM_tot)*28.967; //air
    for (ns = 0; ns < nc; ns++)
    {
        /* gas species index of vaporization */
        int gas_index = TP_COMPONENT_INDEX_I(p, ns);
        if (gas_index >= 0)
        {
            //Ys = x_surf / ( x_surf + (1.-x_surf) * molwt_bulk / solver_par.molWeight[gas_index]);
            Ys = P_USER_REAL(p, ns)* solver_par.molWeight[gas_index] / xs_tot;//!!
            Y_inf += c->yi[gas_index];
            //L_eff += Ys * p->hvap[gas_index];
            L_eff += Ys * p->hvap[gas_index];
            //if (Tp>0.99*659.0) L_eff = Ys * 37.44*pow(659.0 - 653.0, 0.38)*1000.0;;
            Ys_tot += Ys;
            P_USER_REAL(p, nc + ns) = Ys;
        }
    }
    L_eff = L_eff / Ys_tot;
    Tp = P_USER_REAL(p, 4 * nc + 7 + N_INT);
    T_ref = (c->temp + 2.0*Tp) / 3.0;
    c_p_die = (0.2979 + 1.4394*(T_ref / 300.0) - 0.1351*(T_ref / 300.0)*(T_ref / 300.0))*1000.0; //vapour
    //c_p_die = c->Cp;
    //D = 0.527*pow(T_ref / 300.0, 1.583) / c->pressure;
    D = DPM_BINARY_DIFFUSIVITY(p, cond_c, Tp);
    rho_gas_s = c->pressure / (287.01625988193461525183829875375*T_ref);
    Sc = c->mu / (rho_gas_s * D);

    //  BM = (Ys_tot - Y_inf) / (1.0 - Ys_tot);
    BM = (Ys_tot) / (1.0 - Ys_tot);
    FBM = pow(1.0 + BM, 0.7)*log(1.0 + BM) / BM;
    Sh_Star = 2.0 + (pow(1.0 + Re*Sc, 1.0 / 3.0)*MAX(1.0, pow(Re, 0.077)) - 1.0) / FBM;
    Sh = log(1.0 + BM)*Sh_Star;
    //Sh = log(1.0 + BM)*(2.0 + 0.6*sqrt(Re)*pow(Sc, 1.0 / 3.0));
    tot_vap_rate = Ap * D * rho_gas_s * Sh / Dp;
    P_USER_REAL(p, 4 * nc) = Ys_tot;
    P_USER_REAL(p, 4 * nc + 1) = tot_vap_rate;


    BT = BM;
    BT_i = BT;
    dif = 1.0;
    coef = c_p_die * rho_gas_s * D / kgas * Sh_Star;
    phi = 0.e-15;
    while (dif > ACCURACY)
    {
        FBT = pow(1.0 + BT, 0.7)*log(1.0 + BT) / BT;
        Nu_star = 2.0 + (pow(1.0 + Re*Pr, 1.0 / 3.0)*MAX(1.0, pow(Re, 0.077)) - 1.0) / FBT;
        phi = coef / Nu_star;
        BT = pow(1.0 + BM, phi) - 1.0;
        dif = abs(BT - BT_i);
        BT_i = BT;
    }

    Nu = log(1.0 + BT) * Nu_star / BT;

    //dydt[0] -= L_eff * tot_vap_rate / (mp * p->Cp);
    //dzdt->energy += L_eff * tot_vap_rate;





    Visc_l = 1.e-3*exp(2.0303*(300.0 / T_av)*(300.0 / T_av) + 1.1769*(300.0 / T_av) - 2.929);
    k_l = 0.1405 - 0.00022*(T_av - 300.0);
    k_l *= 1000.0;
    C_pl = p->Cp;

    rel_vel = sqrt((c->V[0] - P_VEL(p)[0])*(c->V[0] - P_VEL(p)[0]) + (c->V[1] - P_VEL(p)[1])*(c->V[1] - P_VEL(p)[1]) + (c->V[2] - P_VEL(p)[2])*(c->V[2] - P_VEL(p)[2]));
    Pe = 12.69 / 16.0*P_RHO(p)*0.5*Dp* C_pl / k_l*rel_vel*c->mu / Visc_l*pow(Re, 1.0 / 3.0) / (1.0 + BM);
    k_eff = (1.86 + 0.86*tanh(2.225*log10(Pe / 30.0)))*k_l;
    if (abs(Pe) < 1.e-12) k_eff = k_l;

    T_eff = c->temp - tot_vap_rate*L_eff / PI / Dp / Nu / kgas;
    h0 = kgas*Nu*0.5 / k_eff - 1.0;
    zeta = (h0 + 1.0)*T_eff;
    kappa = k_eff / (C_pl*P_RHO(p)*0.25*Dp*Dp);


    for (i = 0; i < N_Lambda; i++) lambda[i] = -1.0;
    ret = Lambda(h0, lambda);
    for (i = 0; i < N_Lambda; i++) series[i] = 0.e-15;

    for (i = 0; i < N_Lambda; i++)
    {
        b_n = 0.5*(1.0 + h0 / (h0*h0 + lambda[i] * lambda[i]));
        I_n = 0.e-15;
        I_n = P_USER_REAL(p, 4 * nc + 7 + N_INT)*sin(lambda[i]);
        for (j = 1; j < N_INT; j += 2) I_n += 4.0 * P_USER_REAL(p, 4 * nc + 7 + j)*(((double)j)*Delta_R)*sin(lambda[i] * ((double)j)*Delta_R);
        for (j = 2; j < N_INT; j += 2) I_n += 2.0 * P_USER_REAL(p, 4 * nc + 7 + j)*(((double)j)*Delta_R)*sin(lambda[i] * ((double)j)*Delta_R);
        I_n = I_n*Delta_R / 3.0;
        series[i] = (I_n - sin(lambda[i]) / lambda[i] / lambda[i] * zeta)*exp(0.0 - kappa*lambda[i] * lambda[i] * P_DT(p)) / b_n;
    }

    for (j = 0; j < N_INT + 1; j++) P_USER_REAL(p, 4 * nc + 7 + j) = T_eff;
    for (i = 0; i < N_Lambda; i++)
    {
        P_USER_REAL(p, 4 * nc + 7) += series[i] * lambda[i];
        for (j = 1; j < N_INT + 1; j++) P_USER_REAL(p, 4 * nc + 7 + j) += series[i] * sin(lambda[i] * ((double)j)*Delta_R) / (((double)j)*Delta_R);
    }


    Tp = P_USER_REAL(p, 4 * nc + 7 + N_INT);
    T_ref = (c->temp + 2.0*Tp) / 3.0;
    c_p_die = (0.2979 + 1.4394*(T_ref / 300.0) - 0.1351*(T_ref / 300.0)*(T_ref / 300.0))*1000.0;
    D = 0.527*pow(T_ref / 300.0, 1.583) / c->pressure;
    rho_gas_s = c->pressure / (287.01625988193461525183829875375*T_ref);
    Sc = c->mu / (rho_gas_s * D);


    for (ns = 0; ns < nc; ns++)
    {
        /* gas species index of vaporization */
        int gas_index = TP_COMPONENT_INDEX_I(p, ns);
        if (gas_index >= 0)
        {
            vap_rate = P_USER_REAL(p, nc + ns) * tot_vap_rate / Ys_tot; //!!

            if (!p->in_rk)
                if (ABS(vap_rate)>0.)
                    p->limiting_time = MIN(p->limiting_time, dpm_par.fractional_change_factor_mass*P_MASS(p) / vap_rate*TP_COMPONENT_I(p, ns));

            P_USER_REAL(p, 2 * nc + ns) = vap_rate;
            dydt[1 + ns] -= vap_rate;

            {
                int source_index = injection_par.yi2s[gas_index];
                if (source_index >= 0)
                {
                    dzdt->species[source_index] += vap_rate;
                    //p->source.mtc[source_index] = c->rho * Ap * Sh_Star * D / Dp;
                    p->source.mtc[source_index] = c->rho * PI * Dp * Sh_Star * D;
                }
            }
        }
    }


    //p->source.htc = Nu*kgas*P_DIAM(p);
    p->source.htc = 0.e-15;
    dh_dt = Nu * kgas * Ap / Dp * (c->temp - T_av);

    h = Nu * kgas / Dp;
    convective_heating_rate = h * Ap / (mp * p->Cp);

    //dydt[0] += dh_dt / (mp * p->Cp);
    //dzdt->energy -= dh_dt;
    /* limit for higher heating rate */
    if (!p->in_rk)
    {
        if (ABS(convective_heating_rate)>DPM_SMALL)
        {
            factor = dpm_par.fractional_change_factor_heat;
            if (ABS(c->temp - Tp)>Tp)
                factor = dpm_par.fractional_change_factor_heat*Tp / (c->temp - Tp);
            p->limiting_time = MIN(p->limiting_time, factor / ABS(convective_heating_rate));
        }
    }
    P_USER_REAL(p, 4 * nc + 2) = BM;
    P_USER_REAL(p, 4 * nc + 3) = BT;
    P_USER_REAL(p, 4 * nc + 4) = L_eff; // used in temperature calculations
    //P_USER_REAL(p, 4 * nc + 5) = P_DIAM(p);
    //P_USER_REAL(p, 4 * nc + 6) = Dp;
    P_USER_REAL(p, 4 * nc + 5) = Nu; //used in temperature calculations
    P_USER_REAL(p, 4 * nc + 6) = T_av;
    p->state.temp = T_av;
    P_USER_REAL(p, 4 * nc + 7 + N_INT + 1) = P_DIAM(p);
    P_USER_REAL(p, 4 * nc + 7 + N_INT + 2) = DPM_DIAM_FROM_VOL(mp / P_RHO(p));
    //P_USER_REAL(p, 4 * nc + 7 + N_INT + 3) = k_l;
    //P_USER_REAL(p, 4 * nc + 7 + N_INT + 4) = C_pl;
    P_USER_REAL(p, 3 * nc + 0) = Pe;
    dydt[0] = 0.e-15;
}
