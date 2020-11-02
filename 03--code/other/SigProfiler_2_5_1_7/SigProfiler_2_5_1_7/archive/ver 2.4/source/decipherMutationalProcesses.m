function [input, allProcesses, allExposures, idx, processes, ...
          exposures, processStab, processStabAvg] = decipherMutationalProcesses( totalIterationsPerCore, ...
                                                                                 numberProcessesToExtract, ...
                                                                                 input, outputFileNameFull,...
                                                                                 outputFileNameSkinny, outputFileTextSignatures, ...
                                                                                 outputFileTextActivities, job)
   
   %% Define function specific constants
   processesDistance = 'cosine';
   removeWeakMutationTypes = 0.00; % removes weak mutation types, i.e. reduces the dimmensions
   totalCores = job.NumWorkers;
   
   %% Minization options
   saOptions = optimset( 'Display', 'off', 'TolFun', 1e-100, ...
                         'MaxFunEvals', Inf, 'MaxIter', 100000, ...
                         'Algorith', 'interior-point', 'FinDiffType', 'central', ...
                         'TolCon', 1e-100, 'TolX', 1e-100 );

   %% Simple validation of the input params
   if ( totalCores == 0)
       error( 'decipherMutationalProcesses: Please initialize a matlabpool!' );
   end
   
   if ( exist('totalIterationsPerCore', 'var') == 0 )
      error( 'decipherMutationalProcesses: Please specify the number of iterations that need to be performed for each lab!' );
   end
   
   if ( exist('numberProcessesToExtract', 'var') == 0 )
      error( 'decipherMutationalProcesses: Please specify the number of processes that need to be extracted!' );
   end
   
   if ( exist('input', 'var') == 0 )
      error( 'decipherMutationalProcesses: Please specify an input data structure!' );
   end
     
   if ( exist('outputFileNameFull', 'var') == 0 )
       error( 'decipherMutationalProcesses: Please specify an output file name for the full file!' );
   end
   
   if ( exist('outputFileNameSkinny', 'var') == 0 )
       error( 'decipherMutationalProcesses: Please specify an output file name for the skinny file!' );
   end
   
   if ( isfield(input, 'cancerType') == 0 || isfield(input, 'originalGenomes') == 0 ...
           || isfield(input, 'sampleNames') == 0 || isfield(input, 'subtypes') == 0 ...
           || isfield(input, 'types') == 0)
     error( 'decipherMutationalProcesses: Please specify an input file containing the variables: cancerType, originalGenomes, sampleNames, subtypes, types!' );
   end
   
   %% Extracting mutational signatures
   [allProcesses, allExposures, genomeErrors, genomesReconstructed, idx, ...
    idxS, processes, processesStd, exposures, exposureStd, processStab, ...
    processStabAvg, clusterCompactness] = ... ...
        deconvolute(input.originalGenomes, totalIterationsPerCore, numberProcessesToExtract, processesDistance, removeWeakMutationTypes, totalCores);
   
   parfor j = 1 : size(exposures, 2)
     exposures(:,j) = remove_all_single_signatures(exposures(:,j), processes, input.originalGenomes(:, j), saOptions);
   end
 
   %% Saving to MATLAB output files
   save(outputFileNameFull, ...
        'allProcesses', 'allExposures', 'genomeErrors', 'genomesReconstructed', 'idx', ...
        'idxS', 'processes', 'processesStd', 'exposures', 'exposureStd', 'processStab', ...
        'processStabAvg', 'clusterCompactness', 'input');

   save(outputFileNameSkinny, ...
        'processes', 'processesStd', 'exposures', ...
        'exposureStd', 'processStabAvg', 'input'); 
    
   %% Saving to text output files
   exportSignaturesPatterns(processes, input.types, input.subtypes, ...
                            outputFileTextSignatures, ',');
                        
   exportExposures(exposures, input.sampleNames, ...
                   outputFileTextActivities, ',')                          
                              
end