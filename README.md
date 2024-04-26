# AiNU 
AiNU is a material parameter extraction platform by physics-based machine learning (latest version: v1.0) <br>
© Hualin Zhan, [The ANU perovskite PV group](https://www.perovskitegroup.com.au/), Australian National University <br> 
Acknowledgements: [ACAP](https://www.acap.org.au/), [ARENA](https://arena.gov.au/)

This program enables material parameter extraction from the analysis of experiments/characterizations. In this program, different experiments can be analyzed and the users can select their favourite theoretical model for experiment analysis. 

This program is free to use for academic purposes only. However, a valid licenses for MATLAB is required. For the full-physics model, a COMSOL (≥ v5.5) license is also required.

* Version v1.0 features:
  1. Added the Time-Resolved Photoluminescence (TRPL) analysis module of a semiconductor (1D material model). <br>
  2. Three theoretical models are available to choose for TRPL analysis:
     * the kinetic (or ABC) model;
     * the full-physics model (requires a COMSOL license), which describes the drift-diffusion of carriers and the dynamic occupation of defects. <br>
  4. Added a TRPL data pre-process tool, which includes the removal of the noise and interpolation for analysis.
