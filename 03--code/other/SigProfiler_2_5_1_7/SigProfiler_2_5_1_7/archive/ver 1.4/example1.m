%% Example for deciphering mutational signatures when the number of 
%% signatures is KNOWN
%
% Ludmil B. Alexandrov
% Cancer Genome Project
% Wellcome Trust Sanger Institute
% la2@sanger.ac.uk
%
% This software and its documentation are copyright 2012 by the
% Wellcome Trust Sanger Institute/Genome Research Limited. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever. 
% Neither the Wellcome Trust Sanger Institute nor Genome Research Limited 
% is responsible for its use, misuse, or functionality.

clear all;
addpath('source/');
addpath('plotting/');
clc;

%% Open matlabpool
if ( matlabpool('size') == 0 )
    matlabpool open; % opens the default matlabpool, if it is not already opened
end

%% Define parameters 
totalSignatures = 4;
iterationsPerCore = 10;
inputFile = 'input/21_WTSI_BRCA_whole_genome_substitutions.mat';
outputFile = 'output/res_example_1_21_WTSI_BRCA_whole_genome_substitutions.mat';

%% Decipher the signatures of mutational processes from catalogues of mutations
[input allProcesses allExposures idx processes exposures processStab processStabAvg] = ...
    decipherMutationalProcesses(iterationsPerCore, totalSignatures, inputFile, outputFile);

%% Plotting the signatures of mutational processes
plotSignatures(processes, input, allProcesses, idx, processStabAvg);

%% Plotting the signature exposures
plotSignaturesExposureInSamples(exposures, input);
