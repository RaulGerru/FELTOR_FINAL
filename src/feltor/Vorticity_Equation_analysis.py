# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 16:41:10 2021
ANALYSIS PROGRAM FOR VORTICITY EQUATION IN FELTOR: 
@author: raulgerru
"""

import netCDF4 as nc
# import numpy as np
import math
import matplotlib.pyplot as plt
import json
import plotly.graph_objects as go
#import plotly.io as pio
#pio.renderers.default = 'svg'  # 'svg'

t = 7


def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


def edge_plot(magnitude, title):
    plt.figure()
    plt.pcolor(rho[(rho > rho_min) & (rho < rho_max)], (eta - math.pi) / math.pi,
               magnitude[:, (rho > rho_min) & (rho < rho_max)], cmap='jet')
    plt.axvline(x=1, color='k', linestyle='--')
    plt.axhline(-0.5, color='w', linestyle='--')
    plt.axhline(0, color='w', linestyle='--')
    plt.axhline(0.5, color='w', linestyle='--')
    plt.autoscale(enable=True)
    plt.colorbar()
    plt.xlabel('$\\rho $')
    ylabels = ('DOWN', 'HFS', 'UP', 'LFS', 'DOWN')
    y_pos = [-1, -0.5, 0, 0.5, 1]
    plt.yticks(y_pos, ylabels)
    plt.title(title)
    plt.show()


def radial_plot(magnitude):
    plt.figure()
    plt.plot(rho, magnitude[:][0])
    plt.plot(rho, magnitude[:][480])
    plt.plot(rho, magnitude[:][480 * 2])
    plt.plot(rho, magnitude[:][480 * 3])
    plt.axvline(x=1, color='k', linestyle='--')
    plt.legend(['DOWN', 'HFS', 'UP', 'LFS'])
    plt.autoscale(enable=True, axis='y')
    plt.xlim([rho_min, rho_max])


def edge_plot_alt(magnitude):
    fig = go.Figure(data=go.Contour(magnitude))
    fig.show()


# fn="Final_Test_1X_simple_diagRaul_FINAL2.nc"
fn = "TA_test_corrected_2_diag.nc"
ds = nc.Dataset(fn)
inputfile = ds.inputfile
inputfile_json = json.loads(inputfile)

# Physical constants
e = 1.60218e-19  # C
m_H = 1.6726219e-27  # Kg
m_e = 9.10938356e-31  # Kg
eps_0 = 8.85418781e-12  #
mu_0 = 1.25663706e-6  #

# INPUT PARAMETERS (Include by hand)
R_0 = 0.56  # m
B_0 = 0.4  # T
n0 = 1.5e19  # m^-3
m_i = 2  # m_H
T0 = 10  # eV
# tau=tau=inputfile_json['physical']['tau']
tau = 7
Ti = T0 * tau

Omega_0 = e * B_0 / (m_i * m_H)
C = e * n0 * Omega_0

rho = ds['rho'][:]
eta = ds['eta'][:]  # Poloidal direction (from 0 to 2pi)
time = 1e3 * ds['time'][:] / Omega_0
density = ds['electrons_2dX'][:][t]

rho_min = 0.6
rho_max = 1.1

# IF THERE IS EQUILIBRIUM, LET'S CHECK THE EQUATION

vort_elec = ds['v_vort_E_2dX'][:][t] - ds['v_vort_E_2dX'][:][t - 1]
vort_dielec = ds['v_vort_D_2dX'][:][t] - ds['v_vort_D_2dX'][:][t - 1]
dt_Omega = vort_elec + vort_dielec

edge_plot(dt_Omega, r'$\partial_t\Omega$')

electric_adv = ds['v_adv_E_tt_2dX'][:][t]
electric_adv_main = ds['v_adv_E_main_tt_2dX'][:][t]
electric_adv_alt = ds['v_adv_E_alt_tt_2dX'][:][t]
dielectric_adv = ds['v_adv_D_tt_2dX'][:][t]
dielectric_adv_main = ds['v_adv_D_main_tt_2dX'][:][t]
dielectric_adv_alt = ds['v_adv_D_alt_tt_2dX'][:][t]
advection = electric_adv + dielectric_adv
advection_main = electric_adv_main + dielectric_adv_main
advection_alt = electric_adv_alt + dielectric_adv_alt

edge_plot(electric_adv, r'$\nabla \cdot \nabla \cdot (\omega_E u_E)$')
# edge_plot(electric_adv_main)
# plt.title(r'$\nabla \cdot (\nabla \cdot \omega_E u_E)$')
# edge_plot(electric_adv_alt)
# plt.title(r'$\nabla \cdot(\omega_E \cdot \nabla u_E)$')

edge_plot(dielectric_adv, r'$\nabla \cdot \nabla \cdot (\omega_D u_E)$')
# edge_plot(dielectric_adv_main)
# plt.title(r'$\nabla \cdot (\nabla \cdot \omega_D u_E)$')
# edge_plot(dielectric_adv_alt)
# plt.title(r'$\nabla \cdot(\omega_D \cdot \nabla u_E)$')


LHS = dt_Omega + advection

J_par = ds['v_J_par_tt_2dX'][:][t]

edge_plot(J_par, r'$\nabla \cdot (J_\parallel \hat{b})$')

fluct_1 = ds['v_J_perp_tt_2dX'][:][t]
fluct_2 = ds['v_J_mag_tt_2dX'][:][t]
fluct_3 = ds['v_M_em_tt_2dX'][:][t]
J_b_perp = fluct_1 + fluct_2 - fluct_3

edge_plot(J_b_perp, r'$\nabla \cdot J_{b_\perp}$')

curv_1 = ds['v_J_D_divNK_tt_2dX'][:][t]
curv_2 = ds['v_J_JAK_tt_2dX'][:][t]
curv_3 = ds['v_J_NUK_tt_2dX'][:][t]
J_curv = curv_1 + curv_2 + curv_3

edge_plot(J_curv, r'$\nabla \cdot J_{curv}$')

E_r = ds['RFB_E_r_tt_2dX'][:][t]
# edge_plot(E_r)
# plt.title(r'$E_r$')

dP_dr = ds['RFB_GradPi_tt_2dX'][:][t]
# edge_plot(dP_dr)
# plt.title(r'$\partial P_i/\partial r$')

# edge_plot(E_r+dP_dr)
# plt.title('Difference from easy RFB')
'''
RHS=J_par+J_b_perp+J_curv
RHS_main=J_par+J_b_perp_main+J_curv

edge_plot(LHS_main)
plt.title('LHS_alt')
edge_plot(RHS_main)
plt.title('RHS_alt')

edge_plot(LHS_main-RHS_main)


u_E_tor=ds['u_E_tor_tt_2dX'][:][t]
u_E=ds['u_E_tt_2dX'][:][t]
u_E_pol=np.sqrt(ds['u_E_tt_2dX'][:][5]**2-ds['u_E_tor_tt_2dX'][:][t]**2)

edge_plot(u_E)
plt.title('u_E')
edge_plot(u_E_tor)
plt.title('u_E_tor')
edge_plot(u_E_pol)
plt.title('u_E_pol')

   
    elec_source= C*ds['v_S_E_tt_2dX'][:]
    dielec_source=C*tau*ds['v_S_D_tt_2dX'][:]
    Sources=elec_source+dielec_source
    
    edge_plot(Sources)
    plt.title(r'$\Omega_S$')
 
    
    plt.pcolor(rho[(rho>rho_min) & (rho<rho_max)], (eta-math.pi)/math.pi, curv_1[0][:, (rho>rho_min) & (rho<rho_max)], cmap='jet', shading='auto' )
    plt.axvline(x=1, color='k', linestyle='--')
    plt.axhline(-0.5, color='w', linestyle='--')
    plt.axhline(0, color='w', linestyle='--')
    plt.axhline(0.5, color='w', linestyle='--')
    plt.autoscale(enable=True)
    plt.colorbar()
    plt.xlabel('$\\rho $')
    #plt.ylabel('$\\theta/\\pi$')
    #plt.clim(-1e10, 1e10)
    ylabels=('X-point', 'HFS', 'UP', 'LFS', 'X-point')
    y_pos=[-1, -0.5, 0, 0.5, 1]
    plt.yticks(y_pos, ylabels)
    #plt.xlim(-0.5, 0.5)
'''

'''
#LIST OF CONCLUSIONS

1. VORTICITY RADIAL COMPONENT IS THE MAIN ONE (Except for little fluctuations)
2. Same for ADVECTION terms
3. The main components of the advections are small compared with the complete, although they have structure in the order of 0.03




CONCLUSIONS FROM ADVECTION
1. w_E*nabla u_E is small compared with W_E u_E, and therefore, it's divergence too, so the alt electric term makes sense to be so small
2. the radial component of W_E u_E is small, which makes sense, as u_E will be bigger in the perpendicular direction, not the radial.
3. In the diamagnetic case, the w_D*nabla u_E is smaller than the other component, but larger than w_E*nabla u_E (order 10^-3 compared with 10^-7)
4. W_D u_E is smaller in the radial direction than the total, but still relevant, with a strong poloidal distribution
5. the dielectric advective main term seems to be the radial direction of the radial componend of adv_WD_UE_r, so we can trust it
6. Mem nabla b is extremely small (order 10^-9)





*LIttle difference between grad B and curv term


'''
