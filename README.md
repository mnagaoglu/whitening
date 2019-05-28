# whitening
Order of things:
------------------------------------------------------------------------------
  ConvertPPM2JPGand16to8Bit.m : used once for format conversion <br>
  ProcessRawEyeMovements.m : used once to process raw eye movements. <br>
  WhiteningSpectralAnalysisGPU.m : The MAIN analysis function. Uses GPU. <br>
  WhiteningSpectralAnalysisPar.m : OBSOLETE. do not use. it's here for only bookkeeping reasons. <br>
  ComputeAveragePSDs.m : as the name suggests, it was used to compute average PSD. It must be run again if WhiteningSpectralAnalysisGPU is run again. <br>
  FancyFinalPlots.m : The final plotting script for publication quality images. <br>


To get realistic and naturalistic content, 100 images from http://natural-scenes.cps.utexas.edu/db.shtml were used (~40deg field of view horizontally).


MNA 1/2018  mnagaoglu@gmail.com
MNA 5/2019  re-organized folder structure. converted absolute paths to relative. merged diffusion simulation and data with this repo.



