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
    inputCSVFolder     = 'input/text'; % specify the folder containing the CSV file
    inputCSVFile       = 'Biliary-AdenoCA.96.csv'; % please do include the complete fileName including file extensions
    inputCSVFileType   = '96'; % possible options are '96', '192', '1536', 'dinucs', or 'indels'
    outputFolder       = 'output/Biliary-AdenoCA-96';
    
    analysisSignaturesTextIO(inputCSVFolder, outputFolder, inputCSVFile, inputCSVFileType, ... % input CSV folder; output folder; input SCV file; input file type
                             1, 10, 10, ... % start number of signatures; end number of signatures; total iterations per core (please make at least 1,000 iterations)
                             job); % parallel job variable
toc
