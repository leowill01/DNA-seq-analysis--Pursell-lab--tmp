%% This example requires running the "run_example.m" script
%% Clearing workspaces
close all;
clearvars;
clc;
addpath('source/');

%% Starting the default cluster (PLEASE USE AT LEAST 100 CORES)
if ( ~isempty(gcp('nocreate')) )
    delete(gcp);
end
c = parcluster; 
job = parpool(c);

%% Parameters of hierarchical analysis
inputFile     = 'output/21_WGS_BRCA/full/res_21_breast_WGS_substitutions_full_signatures_4.mat'; % Identified file with optimal number of mutational signatures
input         = load(inputFile);
cosineCutoff  = 0.99; % minimum cosine similarity to consider a sample explained
minMutations  = 1000; % minimum number of mutations in a sample for subsequent examinations

%% Evaluating solution for unexplained samples
recon        = input.processes * input.exposures;
totalSamples = size(recon,2);
accuracy     = zeros(totalSamples, 1);

for i = 1 : totalSamples
  accuracy(i) = 1 - pdist([recon(:,i), input.input.originalGenomes(:,i)]', 'cosine');
end

unexplainedIDs = find( (accuracy<cosineCutoff)' & sum(input.input.originalGenomes)>minMutations );

%% Generate new input file with the unexplained samples
if ( length(unexplainedIDs) > 5 ) % do not examine datasets with less than 10 samples
    inputNew = input.input;
    inputNew.originalGenomes = inputNew.originalGenomes(:, unexplainedIDs);
    inputNew.sampleNames     = inputNew.sampleNames(unexplainedIDs);
    inputFolder  = 'temp';
    inputFile    =  char(java.util.UUID.randomUUID);

    if ( exist(inputFolder,'dir') == 0 )
        mkdir(inputFolder);
    end

    fileName = [inputFolder filesep inputFile];
    save(fileName, '-struct', 'inputNew');

    %% Evaluate mutational signatures in unexplained samples
    tic
        outputFolder = 'output/21_WGS_BRCA/L2';
        analysisSignatures(inputFolder, outputFolder, inputFile, ... % input folder; output folder; input file
                           1, 5, 10, ... % start number of signatures; end number of signatures; total iterations per core (please make at least 1,000 iterations)
                           job); % parallel job variable
    toc
else
    disp(['Only ' num2str(length(unexplainedIDs)) ' samples are unexplained! Not enough samples for re-analysis!']);
end