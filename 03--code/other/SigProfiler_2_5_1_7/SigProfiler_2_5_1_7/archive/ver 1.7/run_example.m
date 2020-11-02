%% Clearing workspaces
close all;
clearvars;
clc;
addpath('source/');

%% Starting the default cluster
if ( ~isempty(gcp('nocreate')) )
    delete(gcp);
end
c = parcluster;
job = parpool(c);

%% Perform mutational signatures analysis
tic
    inputFolder  = 'input';
    inputFile    = '21_breast_WGS_substitutions'; % please do not use an extension 
    outputFolder = 'output/21_WGS_BRCA';
    analysisSignatures(inputFolder, outputFolder, inputFile, 1, 5, 10, job);
toc

%% Perform mutational signatures analysis with strand bias
tic
    inputFolder  = 'input';
    inputFile    = '100_breast_WES_substitutions_strand_bias'; % please do not use an extension 
    outputFolder = 'output/100_WES_BRCA';
    analysisSignatures(inputFolder, outputFolder, inputFile, 1, 5, 10, job);
toc
                           