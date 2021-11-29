%% Data for InitialSystem_9_Bus.slx 
clear all;
clc;
close all;

s = tf('s');

%load initial operating states to start simulations at steady-state
load('InitialState_InitialSystem.mat')

%% Simulation Settings
T_s_power = 1e-4;                                    %powergui sample time
Tend=1100;                                           %simulation horizon in seconds
T_load_b6=1010;                                      %load step at bus 6 at time T_load_b6
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

%% Simulation events 
%Load step at bus 6
Load_step_b6 = 30e6;
