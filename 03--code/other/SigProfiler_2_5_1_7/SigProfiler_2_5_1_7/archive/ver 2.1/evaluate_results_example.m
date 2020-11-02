%% This example requires running the "run_example.m" script
%% Clearing workspaces
close all;
clearvars;
clc;
addpath('source/');
addpath('plotting/');

%% Defining input and output params for display
inputFolder = 'output/21_WGS_BRCA'; % for example strand-bias results use: output/100_WES_BRCA
inputFile   = '21_breast_WGS_substitutions'; % for example strand-bias results use: 100_breast_WES_substitutions_strand_bias; please do not use an extension
sigNum = 4; % plot a solution for a particular number of signatures

%% Reading input files
outputFileNameFull = [inputFolder filesep 'full' filesep 'res_' inputFile '_full_signatures_' num2str(sigNum) '.mat'];
outputFileNameSkinny = [inputFolder filesep 'skinny' filesep 'res_' inputFile '_skinny_signatures_' num2str(sigNum) '.mat'];
outputFileNameSummary = [inputFolder filesep 'summary' filesep 'res_' inputFile '_summary.mat'];

%% Plotting individual set of signatures
signaturesFile = load(outputFileNameFull);

if ( (length(signaturesFile.input.types) == 96) || (length(signaturesFile.input.types) == 99) || (length(signaturesFile.input.types) == 1536) )
    plotSignatures(signaturesFile.processes, signaturesFile.input, signaturesFile.allProcesses, signaturesFile.idx, signaturesFile.processStabAvg);
end
if ( (length(signaturesFile.input.types) == 192) )
    plotSignaturesWithStrandBias(signaturesFile.processes, signaturesFile.allProcesses, signaturesFile.idx, signaturesFile.input);
end
plotSignaturesExposureInSamples(signaturesFile.exposures, signaturesFile.input);

%% Plotting summary file
summary = load(outputFileNameSummary);
plotSignatureStabilityAndReconstruction(summary.minSignatures:summary.maxSignatures, summary.avgStability, summary.avgReconstructionError, summary.input);

