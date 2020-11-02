%% Clearing workspaces
close all;
clearvars;
clc;
addpath('source/');
addpath('source/text_import/');

%% Starting the default cluster (PLEASE USE AT LEAST 100 CORES)
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
    analysisSignatures(inputFolder, outputFolder, inputFile, ... % input folder; output folder; input file
                       1, 10, 10, ... % start number of signatures; end number of signatures; total iterations per core (please make at least 1,000 iterations)
                       job); % parallel job variable
toc

%% Perform mutational signatures analysis with strand bias
tic
    inputFolder  = 'input';
    inputFile    = '100_breast_WES_substitutions_strand_bias'; % please do not use an extension 
    outputFolder = 'output/100_WES_BRCA';
    analysisSignatures(inputFolder, outputFolder, inputFile, ... % input folder; output folder; input filee
                       1, 10, 10, ... % start number of signatures; end number of signatures; total iterations per core (please make at least 1,000 iterations)
                       job); % parallel job variable
toc