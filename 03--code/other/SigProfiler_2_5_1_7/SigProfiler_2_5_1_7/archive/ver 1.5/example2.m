%% Example of identifying the number of mutational signatures operative in
%% a set of mutational catalogues. The program will perform extraction with
%% different number of signatures and save a file for each of them. The files
%% could be used for plotting the data (see example3.m)
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

%% Starting the default parallel cluster
if ( ~isempty(gcp('nocreate')) )
    delete(gcp);
end
c = parcluster;
job = parpool(c);

%% Define parameters
iterationsPerCore = 10;
minNumberOfSignature = 1;
maxNumberOfSignature = ;
stability = zeros(maxNumberOfSignature, 1);
reconstructionError = zeros(maxNumberOfSignature, 1);
outputFolder = 'output';
inputFile = '21_WTSI_BRCA_whole_genome_substitutions.mat';
allOutputFile = 'res_example_2_all_21_WTSI_BRCA_whole_genome_substitutions.mat';

%% Sequentially deciphering signatures between minNumberOfSignature and maxNumberOfSignature
for totalSignatures = minNumberOfSignature : maxNumberOfSignature
    
    % Decipher the signatures of mutational processes from catalogues of mutations
    [input, allProcesses, allExposures, idx, processes, exposures, processStab, processStabAvg] = ...
        decipherMutationalProcesses(iterationsPerCore, totalSignatures, inputFile, ...
            ['res_example_2_21_WTSI_BRCA_whole_genome_substitutions_' num2str(totalSignatures) '_signatures.mat'], job);
    
    % Record the stability and average Frobenius reconstruction error
    stability(totalSignatures-minNumberOfSignature+1) = mean(processStabAvg);
    reconstructionError(totalSignatures-minNumberOfSignature+1) = norm(input.originalGenomes - processes*exposures, 'fro');
    
end

%% Plotting the stability and average Frobenius reconstruction error
try %% Some old versions of MATLAB plotyy has a bug under linux with -nodisplay -nosplash -nodesktop options
  plotSignatureStabilityAndReconstruction(minNumberOfSignature:maxNumberOfSignature, stability, reconstructionError, input);
catch ME
  %% Do not do anything - just ignore the plot in order to save the final output daya
end

%% Saving the data
if ( exist(outputFolder,'dir') == 0 )
       mkdir(outputFolder);
end
   
save([outputFolder filesep allOutputFile], 'minNumberOfSignature', ...
     'maxNumberOfSignature', 'stability', 'reconstructionError', 'input');
