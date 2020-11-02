%% This example requires running the "run_example.m" script
%% Clearing workspaces
close all;
clearvars;
clc;
addpath('source/');
addpath('plotting/');

%% Defining input and output params for display
% Note that the display options require first running run_example.m or run_example_text_IO.m
% for displaying the results from the 21 breast WGS example analysis use: output/21_WGS_BRCA;
% for displaying the results from the 100 breat WES strand-bias analysis use: output/100_WES_BRCA;
% for displaying the results from the 35 billary WGS text analysis use: output/Biliary-AdenoCA-96;
inputFolder = 'output/21_WGS_BRCA'; 
    
% Note that the display options require first running run_example.m or run_example_text_IO.m
% for displaying the results from the 21 breast WGS example analysis use: 21_breast_WGS_substitutions;
% for displaying the results from the 100 breat WES strand-bias analysis use: 100_breast_WES_substitutions_strand_bias;
% for displaying the results from the 35 billary WGS text analysis use: Biliary-AdenoCA-96;
inputFile   = '21_breast_WGS_substitutions'; 

% plot a solution for a particular number of signatures
sigNum = 5; 

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

