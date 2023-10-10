# ainu
AiNU: a material parameter extraction platform by physics-based machine learning (latest version: 0.1)

© Hualin Zhan, Australian National University

This platform, which is packaged as software here for convenient implementation, enables material parameter extraction from the analysis of experiments/characterizations. In this software, different experiments can be analyzed and the users can select their favorite theoretical model for experiment analysis. 

This software is free to use. However, a valid COMSOL license (≥ 5.5) may be required.

Version 0.1 features:
1, Added the TRPL analysis module of a semiconductor (1D material model).
2, Three theoretical models are available to choose for TRPL analysis: bi-exponential equation, the kinetic (or ABC) model, or the full-physics model (require a COMSOL license). Here the full-physics model describes the drift-diffusion of carriers and the dynamic occupation of defects.
3, Added a TRPL data pre-process tool, which includes the removal of the noise and interpolation for analysis.
