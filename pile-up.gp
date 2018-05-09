# pile-up.gp
# version 0.1
# created 2018-04-20, last modification 2018-05-09
# Plot torque components and total torque acting on migrating hot Jupiters
# by Rene' Heller, heller@mps.mpg.de, Max Planck Institute for Solar System Research, Goettingen
# published at https://github.com/reneheller/pile-up/ under MIT License

# i=1; call "pile-up.gp" 1000 pdf

reset

NUM = ARG1        # [integer]
TERM = ARG2       # pdf / X / eps

if (i==1) count_loss=0; count_live = 0; Hist_max = 1; bin_010_less = 0; array Hist[30]; bin_020_021 = 0; bin_021_022 = 0; bin_022_023 = 0; bin_023_024 = 0; bin_024_025 = 0; bin_025_026 = 0; bin_026_027 = 0; bin_027_028 = 0; bin_028_029 = 0; bin_029_030 = 0; bin_030_031 = 0; bin_031_032 = 0; bin_032_033 = 0; bin_033_034 = 0; bin_034_035 = 0; bin_035_036 = 0; bin_036_037 = 0; bin_037_038 = 0; bin_038_039 = 0; bin_039_040 = 0; bin_040_041 = 0; bin_041_042 = 0; bin_042_043 = 0; bin_043_044 = 0; bin_044_045 = 0; bin_045_046 = 0; bin_046_047 = 0; bin_047_048 = 0; bin_048_049 = 0; bin_049_050 = 0

#### Natural constants and astrophysical quantities ####
M_sun = 1.98892e+30     # [kg]
M_jup = 1.8986e+27      # [kg]
R_sun = 696000000.0     # [m]
AU    = 149598000000.0  # [m]
G     = 6.673e-11       # [m^3 / kg / s^2]
k_B      = 1.3806503 * 10.**(-23.)                         # [m^2 * kg / s^2 / K], Boltzmann constant
h        = 6.626068 * 10.**(-34.)                          # [m^2 * kg / s], Planck constant
c        = 299792458.                                      # [m / s], speed of light
sigma_SB = 2. * pi**5. * k_B**4. / (15. * c**2. * h**3.)   # [kg / s^3 / K^4], Stefan-Boltzmann constant
m_pro    = 1.672621777 * 10**(-27.)                        # [kg] mass of a proton

#### System paramterization ####
Ms  = M_sun
Mp  = M_jup
Rs  = 2*R_sun
kappa = 10**(-7.)                    # [m^2 / kg]

#### Randomized system paramterization ####
X = 1.* invnorm( rand(0) ) # Gaussian distribution of X around 0
Sigma_p0_RAND = 10**(4 + X)     # [kg / m^2]  |  1 kg/m^2 = 0.1 g / cm^2
Y = 0.5 * invnorm( rand(0) )
D_omega_RAND = 10**(-3.25 + Y)
Z = 1. * invnorm( rand(0) )
alpha_RAND = 10**(-2. + Z)
W = 0.55 * invnorm( rand(0) )
mu_RAND = (1.85+W)*m_pro

#### Astrophysical model equations ####
Omega(a) = sqrt( G*(Mp+Ms)/a**3 )                       # [1/s], orbital frequency
Sigma_p(a) = Sigma_p0_RAND * (a / (1*AU))**(-3./2)            # [kg / m^2] , disk gas surface density
tau(a) = kappa*Sigma_p(a)/2.
Tc(a) = ( 3./4 * (tau(a)/2. + 1./sqrt(3) + 1./(3*tau(a))) * 9./(8*sigma_SB)*alpha_RAND*k_B/mu_RAND*Sigma_p(a)*Omega(a) )**(1./3)    # [K], disk midplane temperature
c_s(a) = sqrt( k_B * Tc(a) / mu_RAND )                                         # [m/s], sound speed in the disk
nu(a) = alpha_RAND * (c_s(a))**2 / Omega(a)                                    #      , disk viscosity
L(a) = Mp * sqrt( a*G*(Mp+Ms) )                          #    , planetary orbital angular momentum
adot(a) = -3*nu(a) / (2*a) * 2*Sigma_p(a) * a**2 / Mp    # [m/s], = da / dt
P(a) = sqrt(4*pi**2*a**3/(G*(Mp+Ms)))

#v_K(a) = a**(-1./2)*sqrt(G*(Ms+Mp))                     # [m/s], Keplerian orbital velocity
#h(a) = c_s(a)/v_K(a)                                    # [1], disk scale height (= H/a)

adot(a) = -3*nu(a) / (2*a) * 2*Sigma_p(a) * a**2 / Mp       # [m/s], = da / dt
Gam_t(a) = 3./2*G*Mp**2*Rs**5/a**6*D_omega_RAND
Gam_II(a) = adot(a) / (2*a) * L(a)

Torque_min = 1e35
do for [j=1:1000] { Torque_j = Gam_II(j*0.0001*AU) + Gam_t(j*0.0001*AU); if (abs(Torque_j)<Torque_min) {Torque_min=abs(Torque_j); a_tmb=j*0.0001*AU}  }
if (a_tmb < 0.02*AU) a_tmb = 0.0005*AU

#### Plot layout ####
if (TERM eq "X" && i==1) set term x11 enhanced font "Times,18" size 800,600; set output
if (TERM eq "eps" && i==1) set term postscript lw 2 color enhanced portrait "Times,12" size 6,4.5; PLOT = "torqueRandHist_typeII.eps"; set output PLOT
if (TERM eq "pdf" && i==1) set term pdfcairo color enhanced font "Times,30" size 6,4.5; PLOT = "torqueRandHist_typeII.pdf"; set output PLOT
set border 31 linewidth 0.2
set samples 1000
set tics in scale 1.5 offset 0,0.22 front mirror; unset ztics; unset y2tics; unset x2tics
set xtics ("0.01" 0.01 0, 0.02 1, 0.03 1, 0.04 1, "0.05" 0.05 0, 0.06 1, 0.07 1, 0.08 1, 0.09 1,\
           "0.1" 0.1 0)border in mirror textcolor rgb "#000000"
set ytics 2 border in mirror textcolor rgb "#000000"
set mytics 4
set format x "%.2f" #"10^{%L}"
set format y "%.f" #  y "%.t{/Symbol \327}10^{%02T}"
set xlabel "semimajor axis [AU]" offset 0,0
set ylabel "torque [10^{30} kg m^2 s^{-2}]" offset 0.5,0
set log x
set origin 0,0
set size 1,0.89
set key spacing 1.4 samplen 0.4 L reverse at 0.192, 6.4

set arrow 1 from Rs/AU,-5.23 to Rs/AU,-0.03 nohead dt 2
set arrow 2 from 0.02 ,0.9 to 0.02 ,6.1 nohead dt 2
set label 1 "star surface " at 0.00845, -5.25 rotate by 90
set label 2 "co-rotation" at 0.018,1 rotate by 90
set label 3 "0.02" at 0.018, graph -.052
set label 4 "0.2" at 0.18, graph -.052
set label 5 "{/*1. {/Symbol \260}}" at a_tmb/AU-0.00049-log10(a_tmb/AU / 0.02)*.002 ,-0.013 textcolor rgb "#99000000"
set label 6 "(b)" at 0.0085, 5.9

if (i==1) set multiplot

plot [.0076:.2][-10:7]\
(x >  0.02 ?  Gam_t(x*AU)/1e30 : 1/0) title "{/Symbol G}_t" w l dt 1 lw 0.03 lc rgb "#DD0000FF",\
((x <= 0.02 && x >= Rs/AU) ? -Gam_t(x*AU)/1e30 : 1/0) notitle w l dt 1 lw 0.03 lc rgb "#DD0000FF",\
(x >= Rs/AU ?  Gam_II(x*AU)/1e30 : 1/0) title "{/Symbol G}_{II}" w l dt 1 lw 0.04 lc rgb "#99FF0000",\
(x >  0.02 ?  ( Gam_t(x*AU)+Gam_II(x*AU))/1e30 : 1/0) title "{/Symbol G}_t + {/Symbol G}_{II}" w l dt 1 lw 0.03 lc rgb "#99000000",\
((x <= 0.02 && x >= Rs/AU) ?  (-Gam_t(x*AU)+Gam_II(x*AU))/1e30 : 1/0) notitle w l dt 1 lw 0.04 lc rgb "#99000000"

if (a_tmb/AU <= 0.02) count_loss = count_loss + 1
if (a_tmb/AU >  0.02) count_live = count_live + 1

if (a_tmb/AU >= 0.020 && a_tmb/AU < 0.021 ) bin_020_021 = bin_020_021 + 1
if (a_tmb/AU >= 0.021 && a_tmb/AU < 0.022 ) bin_021_022 = bin_021_022 + 1
if (a_tmb/AU >= 0.022 && a_tmb/AU < 0.023 ) bin_022_023 = bin_022_023 + 1
if (a_tmb/AU >= 0.023 && a_tmb/AU < 0.024 ) bin_023_024 = bin_023_024 + 1
if (a_tmb/AU >= 0.024 && a_tmb/AU < 0.025 ) bin_024_025 = bin_024_025 + 1
if (a_tmb/AU >= 0.025 && a_tmb/AU < 0.026 ) bin_025_026 = bin_025_026 + 1
if (a_tmb/AU >= 0.026 && a_tmb/AU < 0.027 ) bin_026_027 = bin_026_027 + 1
if (a_tmb/AU >= 0.027 && a_tmb/AU < 0.028 ) bin_027_028 = bin_027_028 + 1
if (a_tmb/AU >= 0.028 && a_tmb/AU < 0.029 ) bin_028_029 = bin_028_029 + 1
if (a_tmb/AU >= 0.029 && a_tmb/AU < 0.030 ) bin_029_030 = bin_029_030 + 1
if (a_tmb/AU >= 0.030 && a_tmb/AU < 0.031 ) bin_030_031 = bin_030_031 + 1
if (a_tmb/AU >= 0.031 && a_tmb/AU < 0.032 ) bin_031_032 = bin_031_032 + 1
if (a_tmb/AU >= 0.032 && a_tmb/AU < 0.033 ) bin_032_033 = bin_032_033 + 1
if (a_tmb/AU >= 0.033 && a_tmb/AU < 0.034 ) bin_033_034 = bin_033_034 + 1
if (a_tmb/AU >= 0.034 && a_tmb/AU < 0.035 ) bin_034_035 = bin_034_035 + 1
if (a_tmb/AU >= 0.035 && a_tmb/AU < 0.036 ) bin_035_036 = bin_035_036 + 1
if (a_tmb/AU >= 0.036 && a_tmb/AU < 0.037 ) bin_036_037 = bin_036_037 + 1
if (a_tmb/AU >= 0.037 && a_tmb/AU < 0.038 ) bin_037_038 = bin_037_038 + 1
if (a_tmb/AU >= 0.038 && a_tmb/AU < 0.039 ) bin_038_039 = bin_038_039 + 1
if (a_tmb/AU >= 0.039 && a_tmb/AU < 0.040 ) bin_039_040 = bin_039_040 + 1
if (a_tmb/AU >= 0.040 && a_tmb/AU < 0.041 ) bin_040_041 = bin_040_041 + 1
if (a_tmb/AU >= 0.041 && a_tmb/AU < 0.042 ) bin_041_042 = bin_041_042 + 1
if (a_tmb/AU >= 0.042 && a_tmb/AU < 0.043 ) bin_042_043 = bin_042_043 + 1
if (a_tmb/AU >= 0.043 && a_tmb/AU < 0.044 ) bin_043_044 = bin_043_044 + 1
if (a_tmb/AU >= 0.044 && a_tmb/AU < 0.045 ) bin_044_045 = bin_044_045 + 1
if (a_tmb/AU >= 0.045 && a_tmb/AU < 0.046 ) bin_045_046 = bin_045_046 + 1
if (a_tmb/AU >= 0.046 && a_tmb/AU < 0.047 ) bin_046_047 = bin_046_047 + 1
if (a_tmb/AU >= 0.047 && a_tmb/AU < 0.048 ) bin_047_048 = bin_047_048 + 1
if (a_tmb/AU >= 0.048 && a_tmb/AU < 0.049 ) bin_048_049 = bin_048_049 + 1
if (a_tmb/AU >= 0.049 && a_tmb/AU < 0.050 ) bin_049_050 = bin_049_050 + 1

if (i==NUM) set origin 0.1218,0.798; set size 0.8785,0.2; unset ylabel; unset xlabel; unset arrow; unset label; unset xtics; set border 0; unset ytics;\
Hist[1] = bin_020_021;\
Hist[2] = bin_021_022;\
Hist[3] = bin_022_023;\
Hist[4] = bin_023_024;\
Hist[5] = bin_024_025;\
Hist[6] = bin_025_026;\
Hist[7] = bin_026_027;\
Hist[8] = bin_027_028;\
Hist[9] = bin_028_029;\
Hist[10] = bin_029_030;\
Hist[11] = bin_030_031;\
Hist[12] = bin_031_032;\
Hist[13] = bin_032_033;\
Hist[14] = bin_033_034;\
Hist[15] = bin_034_035;\
Hist[16] = bin_035_036;\
Hist[17] = bin_036_037;\
Hist[18] = bin_037_038;\
Hist[19] = bin_038_039;\
Hist[20] = bin_039_040;\
Hist[21] = bin_040_041;\
Hist[22] = bin_041_042;\
Hist[23] = bin_042_043;\
Hist[24] = bin_043_044;\
Hist[25] = bin_044_045;\
Hist[26] = bin_045_046;\
Hist[27] = bin_046_047;\
Hist[28] = bin_047_048;\
Hist[29] = bin_048_049;\
Hist[30] = bin_049_050;\
do for [j=1:30] {  if (Hist[j]>Hist_max) {Hist_max=1.*Hist[j]}  };\
set object 12 rect from 0.020,0 to 0.021,bin_020_021/Hist_max fc rgb "#000000";\
set object 13 rect from 0.021,0 to 0.022,bin_021_022/Hist_max fc rgb "#000000";\
set object 14 rect from 0.022,0 to 0.023,bin_022_023/Hist_max fc rgb "#000000";\
set object 15 rect from 0.023,0 to 0.024,bin_023_024/Hist_max fc rgb "#000000";\
set object 16 rect from 0.024,0 to 0.025,bin_024_025/Hist_max fc rgb "#000000";\
set object 17 rect from 0.025,0 to 0.026,bin_025_026/Hist_max fc rgb "#000000";\
set object 18 rect from 0.026,0 to 0.027,bin_026_027/Hist_max fc rgb "#000000";\
set object 19 rect from 0.027,0 to 0.028,bin_027_028/Hist_max fc rgb "#000000";\
set object 20 rect from 0.028,0 to 0.029,bin_028_029/Hist_max fc rgb "#000000";\
set object 21 rect from 0.029,0 to 0.030,bin_029_030/Hist_max fc rgb "#000000";\
set object 22 rect from 0.030,0 to 0.031,bin_030_031/Hist_max fc rgb "#000000";\
set object 23 rect from 0.031,0 to 0.032,bin_031_032/Hist_max fc rgb "#000000";\
set object 24 rect from 0.032,0 to 0.033,bin_032_033/Hist_max fc rgb "#000000";\
set object 25 rect from 0.033,0 to 0.034,bin_033_034/Hist_max fc rgb "#000000";\
set object 26 rect from 0.034,0 to 0.035,bin_034_035/Hist_max fc rgb "#000000";\
set object 27 rect from 0.035,0 to 0.036,bin_035_036/Hist_max fc rgb "#000000";\
set object 28 rect from 0.036,0 to 0.037,bin_036_037/Hist_max fc rgb "#000000";\
set object 29 rect from 0.037,0 to 0.038,bin_037_038/Hist_max fc rgb "#000000";\
set object 30 rect from 0.038,0 to 0.039,bin_038_039/Hist_max fc rgb "#000000";\
set object 31 rect from 0.039,0 to 0.040,bin_039_040/Hist_max fc rgb "#000000";\
set object 32 rect from 0.040,0 to 0.041,bin_040_041/Hist_max fc rgb "#000000";\
set object 33 rect from 0.041,0 to 0.042,bin_041_042/Hist_max fc rgb "#000000";\
set object 34 rect from 0.042,0 to 0.043,bin_042_043/Hist_max fc rgb "#000000";\
set object 35 rect from 0.043,0 to 0.044,bin_043_044/Hist_max fc rgb "#000000";\
set object 36 rect from 0.044,0 to 0.045,bin_044_045/Hist_max fc rgb "#000000";\
set object 37 rect from 0.045,0 to 0.046,bin_045_046/Hist_max fc rgb "#000000";\
set object 38 rect from 0.046,0 to 0.047,bin_046_047/Hist_max fc rgb "#000000";\
set object 39 rect from 0.047,0 to 0.048,bin_047_048/Hist_max fc rgb "#000000";\
set object 40 rect from 0.048,0 to 0.049,bin_048_049/Hist_max fc rgb "#000000";\
set object 41 rect from 0.049,0 to 0.050,bin_049_050/Hist_max fc rgb "#000000";\
set label 1 sprintf("%.1f %%", 100.*(1. - 1.*count_loss/NUM)) at 0.055,0.65;\
plot [.0076:.2][0:1] -x notitle

if (i < NUM) print i; i = i+1; reread;\
else unset multiplot; set output; print "\nloss rate ="; print 1.*count_loss/NUM; print sprintf("%.1f", 100.*count_loss/NUM)
