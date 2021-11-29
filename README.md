# DVPP-9-bus
This repository contains the Simulink/MATLAB models of the IEEE 9 bus system used in the paper:

V. Häberle, M. W. Fisher, E. Prieto-Araujo and F. Dörfler, "Control Design of Dynamic Virtual Power Plants - A Divide-and-Conquer Approach" 

There are three different simulation models:
1. InitialSystem_9_Bus.slx contains the Simulink model of the initial 9 bus system 
2. CaseStudyI_9_Bus.slx contains the Simulink model of case study I with DVPP 1 at bus 1
3. CaseStudyII_9_Bus.slx contains the Simulink model of case study II with DVPP1 at bus 1 and DVPP 3 at bus 3

The underlying data for the above three Simulink models is generated in the MATLAB scripts:
1. Data_InitialSystem.m
2. Data_CaseStudyI.m
3. Data_CaseStudyII.m

Each simulation model can be initialized by a predefined initial system state to start the simulation at steady-state. This initial system state is loaded in the beginning of each simulation's data file. 

This source code is distributed in the hope that it will be useful, but without any warranty.

We do request that publications in which this testbed is adopted, explicitly acknowledge that fact by citing the above mentioned paper:

@ARTICLE{...}

@misc{simulink_model, title={{D}ynamic {V}irtual {P}ower {P}lant Case Studies}, howpublished={GitLab repository}, author={V. Häberle}, note = {https://git.ee.ethz.ch/verenhae/dvppcasestudies}, year={2021}}

For further information do not hesitate to contact me: verenhae@ethz.ch
