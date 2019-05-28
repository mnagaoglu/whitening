This repo contains analysis scripts, figures, and some intermediate data for the paper **"Exploration of Functional Consequences of Fixational Eye Movements in the Absence of a Fovea"**.

# Contents
------------------------------------------------------------------------------
  [ConvertPPM2JPGand16to8Bit.m](utils/ConvertPPM2JPGand16to8Bit.m) : used once for format conversion <br>
  [ProcessRawEyeMovements.m](ProcessRawEyeMovements.m) : used once to process raw eye movements. <br>
  [WhiteningSpectralAnalysisGPU.m](WhiteningSpectralAnalysisGPU.m) : The MAIN analysis function. Uses GPU. <br>
  [ComputeAveragePSDs.m](ComputeAveragePSDs.m) : as the name suggests, it was used to compute average PSD. It must be run again if WhiteningSpectralAnalysisGPU is run again. <br>
  [FancyFinalPlots.m](FancyFinalPlots.m) : The final plotting script for publication quality images. <br>
  [DONOTUSE_WhiteningSpectralAnalysisPar.m] : OBSOLETE. do not use. it's here for only bookkeeping reasons. <br>
     
There are also several short scripts under [utils/](utils/) to generate some of the intermediate data or figures.  

We used 100 natural images from http://natural-scenes.cps.utexas.edu/db.shtml (~40deg field of view horizontally).

For questions/comments: mnagaoglu@gmail.com

# Change log
| Author        | Date           | Change  |
|:-------------:|:-------------:|:----- |
| MNA      | 1/2018 | initial commit |
| MNA      | 5/2019      |  re-organized folder structure. converted absolute paths to relative. merged diffusion simulation and data with this repo. |

