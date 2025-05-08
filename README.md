**Infant Sulcal Morphology**

This repository contains the code and data for the paper: **"How infant brains fold: Sulcal deepening is linked to development of sulcal span, thickness, curvature, and microstructure."** Tung, S.S., Yan, X., Fascendini, B., Tyagi, C., Reyes, C.M., Ducre, K., Perez, K., Allen, A., Horenziak, J., Wu, H., Keil, B., Natu, V.S., Grill-Spector, K. (2025).

https://www.biorxiv.org/content/10.1101/2025.04.18.649604v1#:~:text=Cortical%20folding%20begins%20in%20utero,postnatal%20sulcal%20development%20remain%20unknown.

The repository is organized as follows:

``code/``: Contains analysis scripts for each figure in the paper
- MATLAB scripts for all figures except Figure 4B
- R script for Figure 4B

``data/``: Contains the necessary data files (in CSV format) to run each analysis
- Data files are organized by figure number

**Requirements**
- MATLAB
- RStudio

**Dependencies**

MATLAB Version: Code was written and tested using MATLAB R2021b and R2023a, but it should be compatible with both older and more recent MATLAB versions.

Optional Tools: 
- If you wish to process your own raw data using the scripts (e.g., from NIfTI files), you will need to install Vistasoft, a MATLAB toolbox for neuroimaging data analysis: https://github.com/vistalab/vistasoft
- For infant brain surface reconstruction and analysis, we used a specialized version of FreeSurfer: freesurfer-linux-centos7_x86_64-infant-dev-4a14499-20210109. Note that this is a developmental version specifically optimized for infant brains. More information about infant Freesurfer can be found at : https://surfer.nmr.mgh.harvard.edu/fswiki/infantFS
