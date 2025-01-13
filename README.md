# AiNU 
AiNU is a material parameter extraction platform by physics-based machine learning (latest version: v1.0) <br>
Developer: Dr Hualin Zhan, [The nexSAS group](https://www.nexsas.org/), Australian National University <br>
Acknowledgements: [ACAP](https://www.acap.org.au/), [ARENA](https://arena.gov.au/) <br>
Reference: [DOI: 10.1039/D4EE00911H](https://doi.org/10.1039/D4EE00911H) <br>
Manual: [DOI: 10.5281/zenodo.11097885](https://doi.org/10.5281/zenodo.11097885)

This program enables material parameter extraction from the analysis of experiments/characterizations. In this program, different experiments can be analyzed and the users can select their favourite theoretical model for experiment analysis. 

This program is **free for academic purposes only**. However, a valid license for MATLAB is required. For the full-physics model, a COMSOL (≥ v5.5) license is also required.

The source code of the program will be released soon. Early access to the source code for collaboration is possible by reaching out to us!

## Change log

* Version v1.0 features:
  1. Added the Time-Resolved Photoluminescence (TRPL) analysis module of a semiconductor (1D material model). <br>
  2. Three theoretical models are available to choose for TRPL analysis:
     * the kinetic (or ABC) model;
     * the full-physics model (requires a COMSOL license), which describes the drift-diffusion of carriers and the dynamic occupation of defects. <br>
  4. Added tools for TRPL experimental data preparation, which includes the removal of the noise and interpolation for analysis, and AiNU result visualization.

## Disclaimer

The program is provided "as is", without warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement. In no event shall the authors or copyright holders be liable for any claim, damages or other liability, whether in an action of contract, tort or otherwise, arising from, out of or in connection with the program or the use or other dealings in the program.
	
## Copyright
	
© Hualin Zhan. All rights reserved.
	
Redistribution and use in source and binary forms, with or without modification, for academic purposes only are permitted provided that the following conditions are met:
  1. This program is used for academic purposes only, which includes higher education and academic research. Use for commercial purposes is not permitted.  <br>
  2. Redistributions in binary form must reproduce the above disclaimer, the above copyright notice, and this list of conditions in documentation such as the AiNU manual and/or other materials provided with the distribution. <br>
  3. Redistributions of source code must retain the above disclaimer, the above copyright notice, and this list of conditions.
