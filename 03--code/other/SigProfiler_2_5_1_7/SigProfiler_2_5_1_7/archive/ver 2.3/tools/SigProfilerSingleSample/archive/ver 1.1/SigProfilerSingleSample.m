%% Example for running SigProfilerSingleSample 
%% The example examines 10 bladder exomes from TCGA and
%% assigns mutational signatures using the COSMIC signatures, which were released in 2015
%% Clearing all data
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

%% Analysis of all signatures in individual samples
tic
    analysis_individual_samples('input/COSMIC_2015_release_with_30_signatures.mat', ... % set of signatures
                                'input/signatures_in_cancer_types.mat', ... % set of signatures in samples
                                'input/BLCA_example_cancer_samples.mat', ... % set of individual samples for examination
                                'output/', ... % output folder
                                'signatures_in_BLCA_example_cancer_samples.mat' ); % output file
toc