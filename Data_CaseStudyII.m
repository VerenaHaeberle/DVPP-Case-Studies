%% Data for CaseStudyII_9_Bus.slx
clear all;
clc;
close all;

s = tf('s');

%load initial operating states to start simulations at steady-state
load('InitialState_CaseStudyII')

%load participation factors and feedback gains from H-infinity model matching
load('ParticipationMatrices_SISO_DVPP1.mat')         %participation matrices DVPP 1 (hydro supplementation)
load('ParticipationMatrices_MIMO_DVPP3.mat')         %participation matrices DVPP 3 (SG replacement)
load('ControlGains_DVPP1.mat')                       %control gains DVPP 1
load('ControlGains_DVPP3_LPV.mat')                   %control gains DVPP 3 for LPV control

%% Simulation Settings
T_s_power = 1e-4;                                    %powergui sample time
Tend=1100;                                           %simulation horizon in seconds
T_load_b3=1010;                                      %load step at bus 3 at time T_load_b3
T_cloud = 2010;                                      %cloud at bus 3 at time T_cloud
T_ms=0.001;                                          %data sampling rate

%% Network Base Values 
Sb = 100*10^6;                         %base apparent power [VA]
Pb = 100*10^6;                         %base active power [W]
Qb = 100*10^6;                         %base reactive power [VAr]
VbHV = 230*10^3;                       %base high voltage, l-l rms [V]
Zb = VbHV^2/Sb;                        %base impedance [Ohm]
fb = 50;                               %base frequency [Hz]
w_b = fb*2*pi;                         %base angular frequency [rad/s]
Lb = Zb/w_b;                           %base inductance [H]
Cb = 1/(Zb*w_b);                       %base capacitance [F]

%base voltages for MV regime, i.e the SG and converter injection levels
VbMV1 = 16.5*10^3;                     %base medium voltage, l-l rms [V]
VbMV2 = 18.0*10^3;                     %base medium voltage, l-l rms [V]
VbMV3 = 13.8*10^3;                     %base medium voltage, l-l rms [V]

%% Bus Generation Ratios 
r_b1 = 2.5;                            %rating of hydro 1: 2.5*100e6 MVA
r_b2 = 1.92*0.5;                       %rating of SG 2: 1.92*0.5*100e6 MVA
r_b3 = 1.28*0.5;                       %rating of SG 3: 1.28*0.5*100e6 MVA
r_sum = r_b1+r_b2+r_b3;                %total installed system rating: (2.5+1.92*0.5+1.28*0.5)*100e6 MVA

%rated power synchronous generators
SbSG1 = r_b1*100*10^6;                 %rated power SG 1 [VA]
SbSG2 = r_b2*100*10^6;                 %rated power SG 2 [VA]
SbSG3 = r_b3*100*10^6;                 %rated power SG 3 [VA]

%% Pi-Model Transmission Lines 
%line resistance [pu] 
r_s14 = 0;                              %line 1->4 [pu]
r_s46 = 0.0170;                         %line 4->5 [pu]
r_s69 = 0.039;                          %line 5->6 [pu]
r_s39 = 0;                              %line 3->6 [pu]
r_s98 = 0.0119;                         %line 6->7 [pu]
r_s87 = 0.0085;                         %line 7->8 [pu]
r_s72 = 0;                              %line 8->2 [pu]
r_s75 = 0.032;                          %line 8->9 [pu]
r_s54 = 0.01;                           %line 9->4 [pu]

%line resistance [Ohm]
R_s14 = r_s14*Zb;                       %line 1->4 [Ohm]
R_s46 = r_s46*Zb;                       %line 4->5 [Ohm]
R_s69 = r_s69*Zb;                       %line 5->6 [Ohm]
R_s39 = r_s39*Zb;                       %line 3->6 [Ohm]
R_s98 = r_s98*Zb;                       %line 6->7 [Ohm]
R_s87 = r_s87*Zb;                       %line 7->8 [Ohm]
R_s72 = r_s72*Zb;                       %line 8->2 [Ohm]
R_s75 = r_s75*Zb;                       %line 8->9 [Ohm]
R_s54 = r_s54*Zb;                       %line 9->4 [Ohm]

%line reactance [pu] 
l_s14 = 0.0576;                         %line 1->4 [pu]
l_s46 = 0.0920;                         %line 4->5 [pu]
l_s69 = 0.1700;                         %line 5->6 [pu]
l_s39 = 0.0586;                         %line 3->6 [pu]
l_s98 = 0.1008;                         %line 6->7 [pu]
l_s87 = 0.0720;                         %line 7->8 [pu]
l_s72 = 0.0625;                         %line 8->2 [pu]
l_s75 = 0.1610;                         %line 8->9 [pu]
l_s54 = 0.0850;                         %line 9->4 [pu]

%line inductance [H]
L_s14 = l_s14*Lb;                       %line 1->4 [H]
L_s46 = l_s46*Lb;                       %line 4->5 [H]
L_s69 = l_s69*Lb;                       %line 5->6 [H]
L_s39 = l_s39*Lb;                       %line 3->6 [H]
L_s98 = l_s98*Lb;                       %line 6->7 [H]
L_s87 = l_s87*Lb;                       %line 7->8 [H]
L_s72 = l_s72*Lb;                       %line 8->2 [H]
L_s75 = l_s75*Lb;                       %line 8->9 [H]
L_s54 = l_s54*Lb;                       %line 9->4 [H]

%shunt capacitance [pu]
b_s14 = 0;                              %line 1->4 [pu]
b_s46 = 0.158;                          %line 4->5 [pu]
b_s69 = 0.358;                          %line 5->6 [pu]
b_s39 = 0;                              %line 3->6 [pu]
b_s98 = 0.209;                          %line 6->7 [pu]
b_s87 = 0.149;                          %line 7->8 [pu]
b_s72 = 0;                              %line 8->2 [pu]
b_s75 = 0.306;                          %line 8->9 [pu]
b_s54 = 0.176;                          %line 9->4 [pu]

%shunt capacitance [F]
B_s14 = b_s14*Cb;                       %line 1->4 [F]
B_s46 = b_s46*Cb;                       %line 4->5 [F]
B_s69 = b_s69*Cb;                       %line 5->6 [F]
B_s39 = b_s39*Cb;                       %line 3->6 [F]
B_s98 = b_s98*Cb;                       %line 6->7 [F]
B_s87 = b_s87*Cb;                       %line 7->8 [F]
B_s72 = b_s72*Cb;                       %line 8->2 [F]
B_s75 = b_s75*Cb;                       %line 8->9 [F]
B_s54 = b_s54*Cb;                       %line 9->4 [F]

%% MV/HV Transformers
%MV/HV transformer 1
St_MVHV1 = Sb;                          %rated power [W]
r1_MVHV1 = 1e-6;                        %winding resistance 1 [pu]
r2_MVHV1 = 1e-6;                        %winding resistance 2 [pu]
l1_MVHV1 = 0;                           %winding inductance 1 [pu]
l2_MVHV1 = l_s14;                       %winding inductance 2 [pu]
rm_MVHV1 = 500;                         %magnetization resistance [pu]
lm_MVHV1 = 500;                         %magnetization inductance [pu]

%MV/HV transformer 2
St_MVHV2 = Sb;                           %rated power [W]
r1_MVHV2 = 1e-6;                         %winding resistance 1 [pu]
r2_MVHV2 = 1e-6;                         %winding resistance 2 [pu]
l1_MVHV2 = 0;                            %winding inductance 1 [pu]
l2_MVHV2 = l_s72;                        %winding inductance 2 [pu]
rm_MVHV2 = 500;                          %magnetization resistance [pu]
lm_MVHV2 = 500;                          %magnetization inductance [pu]

%MV/HV transformer 3
St_MVHV3 = Sb;                           %rated power [W]
r1_MVHV3 = 1e-6;                         %winding resistance 1 [pu]
r2_MVHV3 = 1e-6;                         %winding resistance 2 [pu]
l1_MVHV3 = 0;                            %winding inductance 1 [pu]
l2_MVHV3 = l_s39;                        %winding inductance 2 [pu]
rm_MVHV3 = 500;                          %magnetization resistance [pu]
lm_MVHV3 = 500;                          %magnetization inductance [pu]

%% Impedance loads
%total active power load
pL = 3;                               %total active power system load [pu]  
P_L = pL*Sb;                          %total active system load [W]

%reactive power load (positive Var)
qL = 0.825;                           %total reactive power system load, positive Var [pu]
Q_L = qL*Sb;                          %total reactive power system load, positive Var [Var]

%% Setpoints at generation buses
P_balance = 2.4e6;                                            %additional power required for steady-state frequency of 50 Hz
P_set_bus1 = r_b1/r_sum*P_L;                                  %active power generation bus 1 [W]
P_set_bus2 = r_b2/r_sum*P_L+P_balance;                        %active power generation bus 2 [W]
P_set_bus3 = r_b3/r_sum*P_L;                                  %active power generation bus 3 [W]

p_set_bus1 = P_set_bus1/SbSG1;                                %active power generation bus 1 [pu]
p_set_bus2 = P_set_bus2/SbSG2;                                %active power generation bus 2 [pu]
p_set_bus3 = P_set_bus3/SbSG3;                                %active power generation bus 3 [pu]

%% Turbine parameters
%hydro govenor and turbine parameters 
Rg = 0.03;                                  %droop control gain [pu]
Tg = 0.2;                                   %governor time constant [sec]
Rt = 0.38;                                  %transient droop gain [%]
Tr = 5;                                     %reset time constant [sec]
Tw = 1;                                     %water ways time constant [sec]
HGov_delay = tf([1],[Tg 1]);                %govenor delay transer function
H_tr_droop = tf([Tr 1],[Rt/Rg*Tr 1]);       %transient droop transfer function
HGov_dyn = HGov_delay*H_tr_droop;           %govenor and droop dynamics
H_water_way  = tf([-Tw 1],[0.5*Tw 1]);      %transfer function hydro turbine dynamics
HTurb_tot = H_water_way;                    %LTI system hydro turbine dynamics

%steam govenor and reheat turbine parameters 
Rg = 0.03;                                  %droop control gain [pu]
Tg = 0.2;                                   %governor time constant [sec]
Tch = 0.5;                                  %steam chest time constant [sec]
F = 0.3;                                    %turbine power fraction factor [%]
Trh = 7;                                    %reheat time-constant [sec]
SGov_dyn = tf([1],[Tg 1]);                  %govenor dynamics
STurb_ch = tf([1],[Tch 1]);                 %steam chest turbine dynamics
STurb_rh = tf([F*Trh 1],[Trh 1]);           %reheat turbine dynamics
STurb_reheat = STurb_ch*STurb_rh;           %steam chest and reheat turbine dynamics

%% Converter Parameters
%PLL
kp_pll = 150;                               %proportional gain PLL 
ki_pll = 1700;                              %integral gain PLL
xpll0_pu = w_b/(ki_pll);                    %initial state PLL

%inner current control loop 
kp_i_pu = 0.73;                             %proportional gain of PI
ki_i_pu = 1.19;                             %integral gain of PI 
i_star_max = 1.3;                            %reference current transient magnitude limitation [pu]

%reactive power PI-droop controller gains for DVPP 1 converters
kp_q_pu = 0.05;                             %proportional gain of PI
ki_q_pu = 0.005;                            %integral gain of PI

%% DVPP 1: hydro supplementation (hydro, battery, supercapacitor)

%DVPP ratios
cDVPP_1_1 = 1;                                         %device 1 (hydro), percentage of total DVPP rating 
cDVPP_1_2 = 0.2;                                       %device 2 (battery), percentage of total DVPP rating
cDVPP_1_3 = 0.1;                                       %device 3 (supercapacitor), percentage of total DVPP rating
cDVPP_1 = 1;                                           %c_DVPP_1=1, total rating of DVPP

%%Base values for grid-side converter control
Sb_c = 500*1e3;                                        %rating single converter module [VA]

%Device 2 in DVPP 1: battery converter
Sb_DVPP_1_2 = r_b1*100*10^6*cDVPP_1_2;                 %desired rating aggregated VSC converter [VA]
n_DVPP_1_2=Sb_DVPP_1_2/Sb_c;                           %total number of converter modules required for desired power rating
VbLV_DVPP_1_2 = 10^3;                                  %Base low voltage, l-l rms [V]
V_m_DVPP_1_2 = VbLV_DVPP_1_2*sqrt(2/3);                %Base low voltage, p-n peak [V]
I_m_DVPP_1_2=2/3*Sb_DVPP_1_2/V_m_DVPP_1_2;             %L-n peak Converter base current [A]
Zb_DVPP_1_2 = 3/2*V_m_DVPP_1_2^2/Sb_DVPP_1_2;          %Converter base impedance [Ohm]
R_f_DVPP_1_2=0.02/n_DVPP_1_2;                          %output filter resistance [Ohm] 
L_f_DVPP_1_2=(1/(n_DVPP_1_2))*600*10^-6;               %output filter inductance [H]
r_f_DVPP_1_2=R_f_DVPP_1_2/Zb_DVPP_1_2;                 %output filter resistance [pu]
l_f_DVPP_1_2=L_f_DVPP_1_2*w_b/Zb_DVPP_1_2;             %output filter inductance [pu]

%Device 3 in DVPP 1: supercapacitor converter
Sb_DVPP_1_3 = r_b1*100*10^6*cDVPP_1_3;                 %desired rating aggregated VSC converter [VA]
n_DVPP_1_3=Sb_DVPP_1_3/Sb_c;                           %total number of converter modules required for desired power rating
VbLV_DVPP_1_3 = 10^3;                                  %base low voltage, l-l rms [V]
V_m_DVPP_1_3 = VbLV_DVPP_1_3*sqrt(2/3);                %base low voltage, p-n peak [V]
I_m_DVPP_1_3=2/3*Sb_DVPP_1_3/V_m_DVPP_1_3;             %L-n peak Converter base current [A]
Zb_DVPP_1_3 = 3/2*V_m_DVPP_1_3^2/Sb_DVPP_1_3;          %converter base impedance [Ohm]
R_f_DVPP_1_3=0.02/n_DVPP_1_3;                          %output filter resistance [Ohm] 
L_f_DVPP_1_3=(1/(n_DVPP_1_3))*600*10^-6;               %output filter inductance [H]
r_f_DVPP_1_3=R_f_DVPP_1_3/Zb_DVPP_1_3;                 %output filter resistance [pu]
l_f_DVPP_1_3=L_f_DVPP_1_3*w_b/Zb_DVPP_1_3;             %output filter inductance [pu]

%%LV/MV transformer parameters
St = 1e6;                                              %rating of single transformer module [VA]

%Trafo for device 2 in DVPP 1
m_DVPP_1_2=Sb_DVPP_1_2/St;                             %number of aggregated modules for 1 VSC-transformer
St_MVLV_DVPP_1_2 = St*m_DVPP_1_2;                      %transformer rating [VA]
R1_DVPP_1_2_pu=0.00734/m_DVPP_1_2;                     %low voltage-side resistance [pu]
L1_DVPP_1_2_pu=0.0186/m_DVPP_1_2;                      %low voltage-side inductance [pu]
R2_DVPP_1_2_pu=R1_DVPP_1_2_pu;                         %high voltage-side resistance [pu]
L2_DVPP_1_2_pu=L1_DVPP_1_2_pu;                         %high voltage-side inductance [pu]
Rm_DVPP_1_2_pu=347;                                    %magnetizing resistance [pu]
Lm_DVPP_1_2_pu=156;                                    %magnetizing inductance [pu]


%Trafo for device 3 in DVPP 1
m_DVPP_1_3=Sb_DVPP_1_3/St;                             %number of aggregated modules for 1 VSC-transformer
St_MVLV_DVPP_1_3 = St*m_DVPP_1_3;                      %transformer rating [VA]
R1_DVPP_1_3_pu=0.00734/m_DVPP_1_3;                     %low voltage-side resistance [pu]
L1_DVPP_1_3_pu=0.0186/m_DVPP_1_3;                      %low voltage-side inductance [pu]
R2_DVPP_1_3_pu=R1_DVPP_1_3_pu;                         %high voltage-side resistance [pu]
L2_DVPP_1_3_pu=L1_DVPP_1_3_pu;                         %high voltage-side inductance [pu]
Rm_DVPP_1_3_pu=347;                                    %magnetizing resistance [pu]
Lm_DVPP_1_3_pu=156;                                    %magnetizing inductance [pu]

%%desired aggregate transfer function of DVPP 1 for P-f control
dp_DVPP_1 = 0.03;                                               %droop gain [pu]
Dp_DVPP_1_pu=1/dp_DVPP_1;                                       %inverse droop gain [pu]
tau_des_p_DVPP_1 = 0.2;                                         %time constant [sec]
T_des_DVPP_1_p = -Dp_DVPP_1_pu/((tau_des_p_DVPP_1*s+1));        %aggregated desired transfer function for P-f control
T_des_DVPP_1 = T_des_DVPP_1_p;                                  %T_des of DVPP 1

%%local reference models for each device 
M1T_des_DVPP_1 = (M1_DVPP_1.*T_des_DVPP_1);
M2T_des_DVPP_1 = (M2_DVPP_1.*T_des_DVPP_1);
M3T_des_DVPP_1 = (M3_DVPP_1.*T_des_DVPP_1);

ReferenceSS2_DVPP_1 = ss(M2T_des_DVPP_1,'minimal');
ReferenceSS3_DVPP_1 = ss(M3T_des_DVPP_1,'minimal');

%%DC side parameters
%Battery (device 2 in DVPP 1)
Vdc_n_DVPP_1_2=3*V_m_DVPP_1_2;                                              %DC-voltage reference [V]
C_dc_DVPP_1_2=0.008*(n_DVPP_1_2);                                           %DC-side capacitor [F]
R_dc_DVPP_1_2=(Vdc_n_DVPP_1_2/(0.05*(Sb_DVPP_1_2)/Vdc_n_DVPP_1_2));         %DC-side resistor [Ohm]
tau_dc_DVPP_1_2=0.2;                                                        %battery time constant [sec]
eta_DVPP_1_2= w_b/Vdc_n_DVPP_1_2;
m_p_DVPP_1_2=(2*pi*0.5)/(Sb_DVPP_1_2);
k_dc_DVPP_1_2=eta_DVPP_1_2/(Vdc_n_DVPP_1_2*m_p_DVPP_1_2);                   %proportional gain DC-voltage control

%Supercapacitor (device 3 in DVPP 1)
Vdc_n_DVPP_1_3=3*V_m_DVPP_1_3;                                              %DC-voltage reference [V]
C_dc_DVPP_1_3=0.008*(n_DVPP_1_3);                                           %DC-side capacitor [F]
R_dc_DVPP_1_3=(Vdc_n_DVPP_1_3/(0.05*(Sb_DVPP_1_3)/Vdc_n_DVPP_1_3));         %DC-side resistor [Ohm]
tau_dc_DVPP_1_3=0.01;                                                       %battery time constant [sec]
eta_DVPP_1_3= w_b/Vdc_n_DVPP_1_3;
m_p_DVPP_1_3=(2*pi*0.5)/(Sb_DVPP_1_3);
k_dc_DVPP_1_3=eta_DVPP_1_3/(Vdc_n_DVPP_1_3*m_p_DVPP_1_3);                   %proportional gain DC-voltage control


%%internal active power setpoints DVPP 1
pset_DVPP_1_1 = P_set_bus1/SbSG1;
pset_DVPP_1_2 = 0;
Pset_DVPP_1_2 = pset_DVPP_1_2*Sb_DVPP_1_2;
pset_DVPP_1_3 = 0;                       
Pset_DVPP_1_3 = pset_DVPP_1_3*Sb_DVPP_1_3;

%internal reactive power setpoints DVPP 1 
qset_DVPP_1_1 = 0;                          
qset_DVPP_1_2 = 0;                          
Qset_DVPP_1_2 = qset_DVPP_1_2*Sb_DVPP_1_2;
qset_DVPP_1_3 = 0;                          
Qset_DVPP_1_3 = qset_DVPP_1_3*Sb_DVPP_1_3;

%%Initialization of DVPP 1 devices
T0 = [0 1 1; sin(-2*pi/3) cos(-2*pi/3) 1;sin(2*pi/3) cos(2*pi/3) 1 ];       %park transformation matrix

%Battery (device 2 of DVPP 1)
%output filter voltage v 
vd0_DVPP_1_2 = V_m_DVPP_1_2;                           
vq0_DVPP_1_2 = 0;                                      
vd0_DVPP_1_2_pu = vd0_DVPP_1_2/V_m_DVPP_1_2;          
vq0_DVPP_1_2_pu = vq0_DVPP_1_2/V_m_DVPP_1_2;           
vdq0_DVPP_1_2_0 = [vd0_DVPP_1_2;vq0_DVPP_1_2;0];       
vabc_DVPP_1_2_0 = T0*vdq0_DVPP_1_2_0;             
vabc_DVPP_1_2_0_pu = T0*vdq0_DVPP_1_2_0/V_m_DVPP_1_2;

%output filter current i
id0_DVPP_1_2 = 2/3*Pset_DVPP_1_2/vd0_DVPP_1_2;
iq0_DVPP_1_2 = -2/3*Qset_DVPP_1_2/vd0_DVPP_1_2;
id0_DVPP_1_2_pu = id0_DVPP_1_2/I_m_DVPP_1_2;
iq0_DVPP_1_2_pu = iq0_DVPP_1_2/I_m_DVPP_1_2;
idq0_DVPP_1_2 = [id0_DVPP_1_2;iq0_DVPP_1_2;0];
iabc_DVPP_1_2_0 = T0*idq0_DVPP_1_2;
iabc_DVPP_1_2_0_pu = T0*idq0_DVPP_1_2/I_m_DVPP_1_2;

%reference converter output current i*
id_star0_DVPP_1_2 = id0_DVPP_1_2;
id_star0_DVPP_1_2_pu = id0_DVPP_1_2_pu;
iq_star0_DVPP_1_2 = iq0_DVPP_1_2;
iq_star0_DVPP_1_2_pu = iq0_DVPP_1_2_pu;

%integrator current controller
xid0_DVPP_1_2 = 0;
xiq0_DVPP_1_2 = 0;
xid0_DVPP_1_2_pu = 0;
xiq0_DVPP_1_2_pu = 0;

%integrator power controller
xq0_DVPP_1_2 = iq_star0_DVPP_1_2/ki_q_pu;
xq0_DVPP_1_2_pu = xq0_DVPP_1_2/I_m_DVPP_1_2;

%Supercapacitor (device 3 of DVPP 1)
%output filter voltage v 
vd0_DVPP_1_3 = V_m_DVPP_1_3;                           
vq0_DVPP_1_3 = 0;                                       
vd0_DVPP_1_3_pu = vd0_DVPP_1_3/V_m_DVPP_1_3;            
vq0_DVPP_1_3_pu = vq0_DVPP_1_3/V_m_DVPP_1_3;            
vdq0_DVPP_1_3_0 = [vd0_DVPP_1_3;vq0_DVPP_1_3;0];       
vabc_DVPP_1_3_0 = T0*vdq0_DVPP_1_3_0;                  
vabc_DVPP_1_3_0_pu = T0*vdq0_DVPP_1_3_0/V_m_DVPP_1_3;  

%output filter current i
id0_DVPP_1_3 = 2/3*Pset_DVPP_1_3/vd0_DVPP_1_3;
iq0_DVPP_1_3 = -2/3*Qset_DVPP_1_3/vd0_DVPP_1_3;
id0_DVPP_1_3_pu = id0_DVPP_1_3/I_m_DVPP_1_3;
iq0_DVPP_1_3_pu = iq0_DVPP_1_3/I_m_DVPP_1_3;
idq0_DVPP_1_3 = [id0_DVPP_1_3;iq0_DVPP_1_3;0];
iabc_DVPP_1_3_0 = T0*idq0_DVPP_1_3;
iabc_DVPP_1_3_0_pu = T0*idq0_DVPP_1_3/I_m_DVPP_1_3;

%reference converter output current i*
id_star0_DVPP_1_3 = id0_DVPP_1_3;
id_star0_DVPP_1_3_pu = id0_DVPP_1_3_pu;
iq_star0_DVPP_1_3 = iq0_DVPP_1_3;
iq_star0_DVPP_1_3_pu = iq0_DVPP_1_3_pu;

%integrator current controller
xid0_DVPP_1_3 = 0;
xiq0_DVPP_1_3 = 0;
xid0_DVPP_1_3_pu = 0;
xiq0_DVPP_1_3_pu = 0;

%integrator power controller
xq0_DVPP_1_3 = iq_star0_DVPP_1_3/ki_q_pu;
xq0_DVPP_1_3_pu = xq0_DVPP_1_3/I_m_DVPP_1_3;

%% DVPP 3: SG replacement (wind, PV, STATCOM-battery)
%ratings of DVPP converters
S_DVPP_3 = 1.28*0.5*100e6;                          %total DVPP base value
Sb_DVPP_3_1 = 70.5e6;                               %rating aggregated VSC converter 1 (wind)[VA]
Sb_DVPP_3_2 = 53e6;                                 %rating aggregated VSC converter 2 (PV)[VA]
Sb_DVPP_3_3 = 80e6;                                 %rating aggregated VSC converter 3 (STATCOM & battery)[VA]

%DVPP ratios
cDVPP_3_1 = Sb_DVPP_3_1/S_DVPP_3;                   %device 1, percentage of total DVPP base value 
cDVPP_3_2 = Sb_DVPP_3_2/S_DVPP_3;                   %device 2, percentage of total DVPP base value
cDVPP_3_3 = Sb_DVPP_3_3/S_DVPP_3;                   %device 3, percentage of total DVPP base value
cDVPP_3 = 1;  

scaling1_DVPP_3 = cDVPP_3/(cDVPP_3_1);              %device 1, scaling factor for dynamic participation
scaling2_DVPP_3 = cDVPP_3/(cDVPP_3_2);              %device 1, scaling factor for dynamic participation
scaling3_DVPP_3 = cDVPP_3/(cDVPP_3_3);              %device 1, scaling factor for dynamic participation

%%nominal active power ratings of DVPP devices (primary sources)
Pb_DVPP_3_1 = 37e6;                                 %active power rating primary source 1 (wind) [W]
Pb_DVPP_3_2 = 28e6;                                 %active power rating primary source 2 (PV)[W]
Pb_DVPP_3_3 = 0e6;                                  %active power rating primary source 3 (STATCOM & battery)[W]
Pb_DVPP_3 = Pb_DVPP_3_1+Pb_DVPP_3_2 + Pb_DVPP_3_3;

%nominal active power DC-gains
m_DVPP_3_1p_0_n = Pb_DVPP_3_1/(Pb_DVPP_3);          %device1 (wind)
m_DVPP_3_2p_0_n = Pb_DVPP_3_2/(Pb_DVPP_3);          %device2 (PV)

%%nominal reactive power ratings of DVPP converters
Qb_DVPP_3_1 = 60e6;                                 %reactive power rating primary source 1 [Var]
Qb_DVPP_3_2 = 45e6;                                 %reactive power rating primary source 2 [Var]
Qb_DVPP_3_3 = 80e6;                                 %reactive power rating primary source 3 [Var]
Qb_DVPP_3 = Qb_DVPP_3_1+Qb_DVPP_3_2+Qb_DVPP_3_3;

%nominal reactive power DC-gains
m_DVPP_3_1q_0_n = Qb_DVPP_3_1/(Qb_DVPP_3);          %device1 (wind)
m_DVPP_3_2q_0_n = Qb_DVPP_3_2/(Qb_DVPP_3);          %device2 (PV)
m_DVPP_3_3q_0_n = Qb_DVPP_3_3/(Qb_DVPP_3);          %device3 (STATCOM & Battery)

%%Base values for grid-side converter control
Sb_c = 500*1e3;                                     %rating single converter module [VA]

%Device 1 in DVPP 3: wind turbine converter
n_DVPP_3_1 = Sb_DVPP_3_1/Sb_c;                      %total number of converter modules required for desired power rating
VbLV_DVPP_3_1 = 10^3;                               %base low voltage, l-l rms [V]
V_m_DVPP_3_1 = VbLV_DVPP_3_1*sqrt(2/3);             %base low voltage, p-n peak [V]
I_m_DVPP_3_1 = 2/3*Sb_DVPP_3_1/V_m_DVPP_3_1;        %L-n peak Converter base current [A]
Zb_DVPP_3_1 = 3/2*V_m_DVPP_3_1^2/Sb_DVPP_3_1;       %converter base impedance [Ohm]
R_f_DVPP_3_1=0.02/n_DVPP_3_1;                       %output filter resistance [Ohm] 
L_f_DVPP_3_1=(1/(n_DVPP_3_1))*600*10^-6;            %output filter inductance [H]
r_f_DVPP_3_1=R_f_DVPP_3_1/Zb_DVPP_3_1;              %output filter resistance [pu]
l_f_DVPP_3_1=L_f_DVPP_3_1*w_b/Zb_DVPP_3_1;          %output filter inductance [pu]

%Device 2 in DVPP 3: PV power plant converter
n_DVPP_3_2=Sb_DVPP_3_2/Sb_c;                        %total number of converter modules required for desired power rating
VbLV_DVPP_3_2 = 10^3;                               %base low voltage, l-l rms [V]
V_m_DVPP_3_2 = VbLV_DVPP_3_2*sqrt(2/3);             %base low voltage, p-n peak [V]
I_m_DVPP_3_2=2/3*Sb_DVPP_3_2/V_m_DVPP_3_2;          %L-n peak Converter base current [A]
Zb_DVPP_3_2 = 3/2*V_m_DVPP_3_2^2/Sb_DVPP_3_2;       %converter base impedance [Ohm]
R_f_DVPP_3_2=0.02/n_DVPP_3_2;                       %output filter resistance [Ohm]
L_f_DVPP_3_2=(1/(n_DVPP_3_2))*600*10^-6;            %output filter inductance [H]
r_f_DVPP_3_2=R_f_DVPP_3_2/Zb_DVPP_3_2;              %output filter resistance [pu]
l_f_DVPP_3_2=L_f_DVPP_3_2*w_b/Zb_DVPP_3_2;          %output filter inductance [pu]

%Device 3 in DVPP 3: STATCOM-battery converter
n_DVPP_3_3=Sb_DVPP_3_3/Sb_c;                        %total number of converter modules required for desired power rating
VbLV_DVPP_3_3 = 10^3;                               %base low voltage, l-l rms [V]
V_m_DVPP_3_3 = VbLV_DVPP_3_3*sqrt(2/3);             %base low voltage, p-n peak [V]
I_m_DVPP_3_3=2/3*Sb_DVPP_3_3/V_m_DVPP_3_3;          %L-n peak Converter base current [A]
Zb_DVPP_3_3 = 3/2*V_m_DVPP_3_3^2/Sb_DVPP_3_3;       %converter base impedance [Ohm]
R_f_DVPP_3_3=0.02/n_DVPP_3_3;                       %output filter resistance [Ohm] for R/L = 0.1326
L_f_DVPP_3_3=(1/(n_DVPP_3_3))*600*10^-6;            %output filter inductance [H]
r_f_DVPP_3_3=R_f_DVPP_3_3/Zb_DVPP_3_3;              %output filter resistance [pu]
l_f_DVPP_3_3=L_f_DVPP_3_3*w_b/Zb_DVPP_3_3;          %output filter inductance [pu]


%%LV/MV transformer parameters
St = 1e6;                                          %rating of single transformer module [VA]

%Trafo for device 1 of DVPP 3 
m_DVPP_3_1=Sb_DVPP_3_1/St;                         %number of aggregated modules for 1 VSC-transformer
St_MVLV_DVPP_3_1 = St*m_DVPP_3_1;                  %transformer rating [VA]
R1_DVPP_3_1_pu=0.00734/m_DVPP_3_1;                 %low voltage-side resistance [pu]
L1_DVPP_3_1_pu=0.0186/m_DVPP_3_1;                  %low voltage-side inductance [pu]
R2_DVPP_3_1_pu=R1_DVPP_3_1_pu;                     %high voltage-side resistance [pu]
L2_DVPP_3_1_pu=L1_DVPP_3_1_pu;                     %high voltage-side inductance [pu]
Rm_DVPP_3_1_pu=347;                                %magnetizing resistance [pu]
Lm_DVPP_3_1_pu=156;                                %magnetizing inductance [pu]


%Trafo for device 3 of DVPP 3 
m_DVPP_3_2=Sb_DVPP_3_2/St;                         %number of aggregated modules for 1 VSC-transformer
St_MVLV_DVPP_3_2 = St*m_DVPP_3_2;                  %transformer rating [VA]
R1_DVPP_3_2_pu=0.00734/m_DVPP_3_2;                 %low voltage-side resistance [pu]
L1_DVPP_3_2_pu=0.0186/m_DVPP_3_2;                  %low voltage-side inductance [pu]
R2_DVPP_3_2_pu=R1_DVPP_3_2_pu;                     %high voltage-side resistance [pu]
L2_DVPP_3_2_pu=L1_DVPP_3_2_pu;                     %high voltage-side inductance [pu]
Rm_DVPP_3_2_pu=347;                                %magnetizing resistance [pu]
Lm_DVPP_3_2_pu=156;                                %magnetizing inductance [pu]

%Trafo for device 3 of DVPP 3 
m_DVPP_3_3=Sb_DVPP_3_3/St;                         %number of aggregated modules for 1 VSC-transformer
St_MVLV_DVPP_3_3 = St*m_DVPP_3_3;                  %transformer rating [VA]
R1_DVPP_3_3_pu=0.00734/m_DVPP_3_3;                 %low voltage-side resistance [pu]
L1_DVPP_3_3_pu=0.0186/m_DVPP_3_3;                  %low voltage-side inductance [pu]
R2_DVPP_3_3_pu=R1_DVPP_3_3_pu;                     %high voltage-side resistance [pu]
L2_DVPP_3_3_pu=L1_DVPP_3_3_pu;                     %high voltage-side inductance [pu]
Rm_DVPP_3_3_pu=347;                                %magnetizing resistance [pu]
Lm_DVPP_3_3_pu=156;                                %magnetizing inductance [pu]

%%desired aggregate transfer function of DVPP 3
%specified aggregated transfer function active power
dp_DVPP_3 = 0.03;                                                           %droop gain [pu]
Dp_DVPP_3_pu=1/dp_DVPP_3;                                                   %inverse droop gain [pu]
Mp_DVPP_3_pu = 13;                                                          %virtual inertia gain [sec]
tau_des_p_DVPP_3 = 0.2;                                                     %time constant [sec]
T_des_p_DVPP_3 = (-Dp_DVPP_3_pu-Mp_DVPP_3_pu*s)/((s*tau_des_p_DVPP_3+1));   %aggregated desired transfer function P-f control

%specified aggregated transfer function reactive power 
Dq1_DVPP_3_pu=100;                                                          %inverse droop gain [pu]
tau_des_q_DVPP_3 = 0.2;                                                     %time constant [sec]
T_des_q_DVPP_3 = (-Dq1_DVPP_3_pu)/((s*tau_des_q_DVPP_3+1));                 %aggregated desired transfer function Q-v control

%MIMO transfer function of DVPP 3
T_des_DVPP_3 = [T_des_p_DVPP_3,0;0,T_des_q_DVPP_3];                         %total desired transfer function DVPP 3

%%DC side parameters
%Wind (device 1 in DVPP 3)
Vdc_n_DVPP_3_1=3*V_m_DVPP_3_1;
C_dc_DVPP_3_1=0.008*(n_DVPP_3_1);
R_dc_DVPP_3_1=(Vdc_n_DVPP_3_1/(0.05*(Sb_DVPP_3_1)/Vdc_n_DVPP_3_1));
tau_dc_DVPP_3_1=1.5;
eta_DVPP_3_1= w_b/Vdc_n_DVPP_3_1;
m_p_DVPP_3_1=(2*pi*0.5)/(Sb_DVPP_3_1);
k_dc_DVPP_3_1=eta_DVPP_3_1/(Vdc_n_DVPP_3_1*m_p_DVPP_3_1);

%PV (device 2 in DVPP 3)
Vdc_n_DVPP_3_2=3*V_m_DVPP_3_2;
C_dc_DVPP_3_2=0.008*(n_DVPP_3_2);
R_dc_DVPP_3_2=(Vdc_n_DVPP_3_2/(0.05*(Sb_DVPP_3_2)/Vdc_n_DVPP_3_2));
tau_dc_DVPP_3_2=0.6;
eta_DVPP_3_2= w_b/Vdc_n_DVPP_3_2;
m_p_DVPP_3_2=(2*pi*0.5)/(Sb_DVPP_3_2);
k_dc_DVPP_3_2=eta_DVPP_3_2/(Vdc_n_DVPP_3_2*m_p_DVPP_3_2);

%STATCOM & Battery (Device 3 in DVPP 3)
Vdc_n_DVPP_3_3=3*V_m_DVPP_3_3;
C_dc_DVPP_3_3=0.008*(n_DVPP_3_3);
R_dc_DVPP_3_3=(Vdc_n_DVPP_3_3/(0.05*(Sb_DVPP_3_3)/Vdc_n_DVPP_3_3));
tau_dc_DVPP_3_3=0.2;
eta_DVPP_3_3= w_b/Vdc_n_DVPP_3_3;
m_p_DVPP_3_3=(2*pi*0.5)/(Sb_DVPP_3_3);
k_dc_DVPP_3_3=eta_DVPP_3_3/(Vdc_n_DVPP_3_3*m_p_DVPP_3_3);

%%internal active power setpoints DVPP 3
pset_DVPP_3_1 = (0.6*P_set_bus3)/Sb_DVPP_3_1;
Pset_DVPP_3_1 = pset_DVPP_3_1*Sb_DVPP_3_1;
pset_DVPP_3_2 = (0.4*P_set_bus3)/Sb_DVPP_3_2;
Pset_DVPP_3_2 = pset_DVPP_3_2*Sb_DVPP_3_2;
pset_DVPP_3_3 = 0;                       
Pset_DVPP_3_3 = pset_DVPP_3_3*Sb_DVPP_3_3;

%internal reactive power setpoints DVPP 3 
qset_DVPP_3_1 = 0;           
Qset_DVPP_3_1 = qset_DVPP_3_1*Sb_DVPP_3_1;
qset_DVPP_3_2 = 0;                          
Qset_DVPP_3_2 = qset_DVPP_3_2*Sb_DVPP_3_2;
qset_DVPP_3_3 = 0;                          
Qset_DVPP_3_3 = qset_DVPP_3_3*Sb_DVPP_3_3;


%%Initialization of DVPP 3 devices
T0 = [0 1 1; sin(-2*pi/3) cos(-2*pi/3) 1;sin(2*pi/3) cos(2*pi/3) 1 ];      %Park transformation matrix

%Wind (device 1 of DVPP 3)
%output filter voltage v 
vd0_DVPP_3_1 = V_m_DVPP_3_1;                            
vq0_DVPP_3_1 = 0;                                       
vd0_DVPP_3_1_pu = vd0_DVPP_3_1/V_m_DVPP_3_1;            
vq0_DVPP_3_1_pu = vq0_DVPP_3_1/V_m_DVPP_3_1;           
vdq0_DVPP_3_1_0 = [vd0_DVPP_3_1;vq0_DVPP_3_1;0];        
vabc_DVPP_3_1_0 = T0*vdq0_DVPP_3_1_0;                  
vabc_DVPP_3_1_0_pu = T0*vdq0_DVPP_3_1_0/V_m_DVPP_3_1;   

%output filter current i
id0_DVPP_3_1 = 2/3*Pset_DVPP_3_1/vd0_DVPP_3_1;
iq0_DVPP_3_1 = -2/3*Qset_DVPP_3_1/vd0_DVPP_3_1;
id0_DVPP_3_1_pu = id0_DVPP_3_1/I_m_DVPP_3_1;
iq0_DVPP_3_1_pu = iq0_DVPP_3_1/I_m_DVPP_3_1;
idq0_DVPP_3_1 = [id0_DVPP_3_1;iq0_DVPP_3_1;0];
iabc_DVPP_3_1_0 = T0*idq0_DVPP_3_1;
iabc_DVPP_3_1_0_pu = T0*idq0_DVPP_3_1/I_m_DVPP_3_1;

%reference converter output current i*
id_star0_DVPP_3_1 = id0_DVPP_3_1;
id_star0_DVPP_3_1_pu = id0_DVPP_3_1_pu;
iq_star0_DVPP_3_1 = iq0_DVPP_3_1;
iq_star0_DVPP_3_1_pu = iq0_DVPP_3_1_pu;

%integrator current controller
xid0_DVPP_3_1 = 0;
xiq0_DVPP_3_1 = 0;
xid0_DVPP_3_1_pu = 0;
xiq0_DVPP_3_1_pu = 0;

%PV (device 2 of DVPP 3)
%output filter voltage v 
vd0_DVPP_3_2 = V_m_DVPP_3_2;                            
vq0_DVPP_3_2 = 0;                                       
vd0_DVPP_3_2_pu = vd0_DVPP_3_2/V_m_DVPP_3_2;            
vq0_DVPP_3_2_pu = vq0_DVPP_3_2/V_m_DVPP_3_2;            
vdq0_DVPP_3_2_0 = [vd0_DVPP_3_2;vq0_DVPP_3_2;0];        
vabc_DVPP_3_2_0 = T0*vdq0_DVPP_3_2_0;                   
vabc_DVPP_3_2_0_pu = T0*vdq0_DVPP_3_2_0/V_m_DVPP_3_2;   

%output filter current i
id0_DVPP_3_2 = 2/3*Pset_DVPP_3_2/vd0_DVPP_3_2;
iq0_DVPP_3_2 = -2/3*Qset_DVPP_3_2/vd0_DVPP_3_2;
id0_DVPP_3_2_pu = id0_DVPP_3_2/I_m_DVPP_3_2;
iq0_DVPP_3_2_pu = iq0_DVPP_3_2/I_m_DVPP_3_2;
idq0_DVPP_3_2 = [id0_DVPP_3_2;iq0_DVPP_3_2;0];
iabc_DVPP_3_2_0 = T0*idq0_DVPP_3_2;
iabc_DVPP_3_2_0_pu = T0*idq0_DVPP_3_2/I_m_DVPP_3_2;

%reference converter output current i*
id_star0_DVPP_3_2 = id0_DVPP_3_2;
id_star0_DVPP_3_2_pu = id0_DVPP_3_2_pu;
iq_star0_DVPP_3_2 = iq0_DVPP_3_2;
iq_star0_DVPP_3_2_pu = iq0_DVPP_3_2_pu;

%integrator current controller
xid0_DVPP_3_2 = 0;
xiq0_DVPP_3_2 = 0;
xid0_DVPP_3_2_pu = 0;
xiq0_DVPP_3_2_pu = 0;


%STATCOM & Battery (device 3 of DVPP 3)
%output filter voltage v 
vd0_DVPP_3_3 = V_m_DVPP_3_3;                           
vq0_DVPP_3_3 = 0;                                       
vd0_DVPP_3_3_pu = vd0_DVPP_3_3/V_m_DVPP_3_3;            
vq0_DVPP_3_3_pu = vq0_DVPP_3_3/V_m_DVPP_3_3;            
vdq0_DVPP_3_3_0 = [vd0_DVPP_3_3;vq0_DVPP_3_3;0];        
vabc_DVPP_3_3_0 = T0*vdq0_DVPP_3_3_0;                   
vabc_DVPP_3_3_0_pu = T0*vdq0_DVPP_3_3_0/V_m_DVPP_3_3;   

%output filter current i
id0_DVPP_3_3 = 2/3*Pset_DVPP_3_3/vd0_DVPP_3_3;
iq0_DVPP_3_3 = -2/3*Qset_DVPP_3_3/vd0_DVPP_3_3;
id0_DVPP_3_3_pu = id0_DVPP_3_3/I_m_DVPP_3_3;
iq0_DVPP_3_3_pu = iq0_DVPP_3_3/I_m_DVPP_3_3;
idq0_DVPP_3_3 = [id0_DVPP_3_3;iq0_DVPP_3_3;0];
iabc_DVPP_3_3_0 = T0*idq0_DVPP_3_3;
iabc_DVPP_3_3_0_pu = T0*idq0_DVPP_3_3/I_m_DVPP_3_3;

%reference converter output current i*
id_star0_DVPP_3_3 = id0_DVPP_3_3;
id_star0_DVPP_3_3_pu = id0_DVPP_3_3_pu;
iq_star0_DVPP_3_3 = iq0_DVPP_3_3;
iq_star0_DVPP_3_3_pu = iq0_DVPP_3_3_pu;

%integrator current controller
xid0_DVPP_3_3 = 0;
xiq0_DVPP_3_3 = 0;
xid0_DVPP_3_3_pu = 0;
xiq0_DVPP_3_3_pu = 0;

%%LPV DVPP controller for DVPP 3
tau1_p_DVPP_3 = 1.5;                                        %LPF-time constant device 1 (wind), m_pf
tau2_p_DVPP_3 = 0.6;                                        %LPF-time constant device 2 (PV), m_pf
tau3_p_DVPP_3 = 0.081;                                      %PLL-time constant for HPF cutoff (STATCOM), m_pf

tau1_q_DVPP_3 = 0.081;                                      %PLL-time constant for m_qv (wind)
tau2_q_DVPP_3 = 0.081;                                      %PLL-time constant for m_qv (pv) 
tau3_q_DVPP_3 = 0.081;                                      %PLL-time constant for m_qv (STATCOM)

%extreme points of parameter polytope
%device1
theta1_p_1 = 0;
theta1_q_1 = 0.3109;
theta1_p_2 = 1;
theta1_q_2 = 0.3603;

%device2
theta2_p_1 = 0;
theta2_q_1 = 0.2302;
theta2_p_2 = 1;
theta2_q_2 = 0.2746;

%device3
theta3_q_1 = 0.3931;
theta3_q_2 = 0.4324;

%% Simulation events 
%Load step at bus 3
Load_step_b3 = Pset_DVPP_3_2;
