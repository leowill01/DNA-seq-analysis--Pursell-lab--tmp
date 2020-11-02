function [input, allProcesses, allExposures, idx, processes, ...
          exposures, processStab, processStabAvg] = decipherMutationalProcesses( totalIterationsPerCore, ...
                                                                                 numberProcessesToExtract, ...
                                                                                 input, outputFileNameFull,...
                                                                                 outputFileNameSkinny, job)
   
   %% Define function specific constants
   processesDistance = 'cosine';
   removeWeakMutationTypes = 0.001; % removes weak mutation types, i.e. reduces the dimmensions
   totalCores = job.NumWorkers;

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
   
   %% Saving to the output file
   save(outputFileNameFull, ...
        'allProcesses', 'allExposures', 'genomeErrors', 'genomesReconstructed', 'idx', ...
        'idxS', 'processes', 'processesStd', 'exposures', 'exposureStd', 'processStab', ...
        'processStabAvg', 'clusterCompactness', 'input');

   save(outputFileNameSkinny, ...
        'processes', 'processesStd', 'exposures', ...
        'exposureStd', 'processStabAvg', 'input'); 
end