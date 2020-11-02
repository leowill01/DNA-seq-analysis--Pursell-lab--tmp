%% Clearing workspaces
close all;
clearvars;
clc;
addpath('source/');
addpath('plotting/');

%% Defining input and output params for display
inputFolder = 'output/100_WES_BRCA'; % for example strand-bias results use: output/100_WES_BRCA
inputFile   = '100_breast_WES_substitutions_strand_bias'; % for example strand-bias results use: 100_breast_WES_substitutions_strand_bias
sigNum = 2;

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


%% Plotting summary file
summary = load(outputFileNameSummary);
plotSignatureStabilityAndReconstruction(summary.minSignatures:summary.maxSignatures, summary.avgStability, summary.avgReconstructionError, summary.input);

