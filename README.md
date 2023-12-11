# Data-Driven-Safety-Preserving-Control-Architecture-for-Constrained-CPS

## Description 
This repository contains the codes for "A Data-Driven Safety Preserving Control Architecture for Constrained Cyber-Physical Systems" by [Mehran Attar](https://scholar.google.com/citations?user=nnLTy-oAAAAJ&hl=en) and [Walter Lucia](https://users.encs.concordia.ca/~wlucia/index.html)

## Abstracat 
In this work, we propose a data-driven networked control architecture for unknown constrained linear cyber-physical systems subject to bounded disturbances capable of detecting and ensuring plant safety in the presence of cyber-attacks on the communication channels. To this end, first, by using data-driven forward reachability analysis, a passive anomaly detector that is local to the controller is designed to detect network attacks. Then, to ensure that intelligent and undetectable attacks will be detected before causing safety risks, a safety verification module, which is local to the plant, is designed based on an outer approximation of the one-step evolution set at each time step. This unit is in charge of activating the emergency controller whenever the control input is deemed unsafe or under the effect of cyber-attacks. Finally, an emergency controller is designed by exploiting a family of data-driven Robust One-Step Controllable (ROSC) Sets and dual-mode set-theoretic model predictive controller to safely confine the system into the closest robust control invariant region in a finite number of steps. Numerical simulations involving a two-tank water system is performed to clarify the proposed solution's capabilities.

## Running
1- Download [CORA release 2022](https://tumcps.github.io/CORA/) and [MPT3 release 2022](https://www.mpt3.org/) 

2- Add CORA and MPT folder and subfolders to the Matlab (in this work, [MATLAB R2021-a](https://www.mathworks.com/products/new_products/release2021a.html) has been used) path.

3- Add the repo folder and subfolders to the Matlab path.

## Files Description
Two scenarios have been considered in this work. 
### Scenario 1: Attack on the actuation channel 
1- To simulate this scenario, please run "data_driven_architecture_attack_on_actuation.m"  

### Scenario 2: Attack on the measurement channel 
2- To simulate this scenario, please run "data_driven_architecture_attack_on_measurement.m" 

#### Function Descriptions
- "compute_AB.m": computes all possible system matrices that are consistent with the data
- "computeRPI.m": computes a model-based RCI set based on the proposed method in **"Invariant approximations of the minimal robust positively invariant set", by Rakovic et al.
- "computing_ROSC_sets.m": computes the family of ROSC sets by considering a target set $\hat{\mathcal{T}}^0$
- "compute_intersec.m": computes the intersection of polyhedrons
- "compute_presets_approx.m": computes the data-driven ROSC sets in the extended space of $(x,u)$.
- "indx_finder.m": computes the set memebership index of an state for the model-based ROSC sets. 
- "one_step_ctrl.m" computes the data-driven ST-MPC control commands. 
- "poly_approx.m": computes a zonotopic inner approximation of a polyhedron 
- "set_index.m": computes the set memebership index of an state for the data-driven ROSC sets.
- 


## Animated test 
### Scenario 1: Attack on the actuation channel 
1- please run "data_driven_architecture_attack_on_actuation.m"  

### Scenario 2: Attack on the measurement channel 
2- please run "data_driven_architecture_attack_on_measurement.m" 
