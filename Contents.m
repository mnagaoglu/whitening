function Contents
% Contents.m
%
% Order of things:
% ------------------------------------------------------------------------------
% ConvertPPM2JPGand16to8Bit.m : used once for format conversion
% ProcessRawEyeMovements.m : used once to process raw eye movements.
% WhiteningSpectralAnalysisGPU.m : The MAIN analysis function. Uses GPU.
% WhiteningSpectralAnalysisPar.m : OBSOLETE. do not use. it's here for only
%      bookkeeping reasons.
% ComputeAveragePSDs.m : as the name suggests, it was used to compute
%      average PSD. It must be run again if WhiteningSpectralAnalysisGPU is run
%      again.
% FancyFinalPlots.m : The final plotting script for publication quality
%      images.
%
%
% MNA 1/1/2018 