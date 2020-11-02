function analysis_individual_samples_IO(signaturesSet, signaturesInSamples, ...
                                        samplesForAnalysis96, samplesForAnalysis192, ...
                                        seqType, outputFolder, outputFile, useRules, ...
                                        allowSigsEverywhere, connected)

%% Read text input files and generate MATLAB input
templatesFolder   = ['source' filesep 'text_import' filesep 'templates' filesep];
inputTempFileName = ['input' filesep outputFile '_temp' '.mat'];
import_96_mutation_types_single_sample(samplesForAnalysis96, templatesFolder, ...
                                       inputTempFileName, seqType, samplesForAnalysis192);
disp('Completed the import of input data from text/csv files.');

%% Perform single sample analysis  
inputFile = inputTempFileName;
analysis_individual_samples(signaturesSet, ... % set of signatures
                            signaturesInSamples, ... % set of signatures in samples
                            inputFile, ... % set of individual samples for examination
                            outputFolder, ... % output folder
                            outputFile, ... % output file
                            useRules, ... % boolean variable indicating whether to use rules or not (1==use rules; 0==do not use rules)
                            allowSigsEverywhere, ... % IDs of signatures to be included in all samples regardless of rules or sparsity  
                            connected); % connected signatures (e.g., signatures SBS-2 and SBS-13)

%% Removing temporary input file
% delete(inputTempFileName);

end






