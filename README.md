# AiNU 
AiNU is a material parameter extraction platform by physics-based machine learning (latest version: v1.0.0) <br>
© Hualin Zhan, [The ANU perovskite PV group](https://www.perovskitegroup.com.au/), Australian National University <br> 
Acknowledgements: [ACAP](https://www.acap.org.au/), [ARENA](https://arena.gov.au/)

This platform, which is packaged as software here for convenient implementation, enables material parameter extraction from the analysis of experiments/characterizations. In this software, different experiments can be analyzed and the users can select their favorite theoretical model for experiment analysis. 

This software is free to use. However, a valid COMSOL license (≥ 5.5) is required for the full-physics model.

* Version v1.0.0 features:
  1. Added the Time-Resolved Photoluminescence (TRPL) analysis module of a semiconductor (1D material model). <br>
  2. Three theoretical models are available to choose for TRPL analysis:
     * the bi-exponential equation;
     * the kinetic (or ABC) model;
     * the full-physics model (require a COMSOL license), which describes the drift-diffusion of carriers and the dynamic occupation of defects. <br>
  4. Added a TRPL data pre-process tool, which includes the removal of the noise and interpolation for analysis.
