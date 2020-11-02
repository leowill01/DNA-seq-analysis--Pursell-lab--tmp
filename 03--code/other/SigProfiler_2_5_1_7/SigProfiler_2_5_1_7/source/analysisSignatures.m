function analysisSignatures(inputFolder, outputFolder, inputFile, minSignatures, maxSignatures, iterationsPerCore, job) 

     %% Reading input file
       input = load([ inputFolder filesep inputFile '.mat']);
       if ( isfield(input, 'cancerType') == 0 || isfield(input, 'originalGenomes') == 0 ...
               || isfield(input, 'sampleNames') == 0 || isfield(input, 'subtypes') == 0 ...
               || isfield(input, 'types') == 0)
         error( 'analysisSignatures: Please specify an input file containing the variables: cancerType, originalGenomes, sampleNames, subtypes, types!' );
       end
     totalSamples    = size(input.originalGenomes, 2);
     
     %% Creating output folder
     if ( exist(outputFolder,'dir') == 0 )
       mkdir(outputFolder);
     end
     if ( exist([outputFolder filesep 'full'],'dir') == 0 )
       mkdir([outputFolder filesep 'full']);
     end
     if ( exist([outputFolder filesep 'skinny'],'dir') == 0 )  
       mkdir([outputFolder filesep 'skinny']);
     end
     if ( exist([outputFolder filesep 'summary'],'dir') == 0 )  
       mkdir([outputFolder filesep 'summary']);
     end
     if ( exist([outputFolder filesep 'text'],'dir') == 0 )  
       mkdir([outputFolder filesep 'text']);
     end
   
     %% Allocating memory for summary variables
     totalSignatures = maxSignatures - minSignatures + 1;
     avgReconstructionError = zeros(totalSignatures, 1);
     avgStability = zeros(totalSignatures, 1);
     avgReconstructionErrorPercentage  = zeros(totalSignatures, 1);
     reconError = zeros(totalSignatures, totalSamples);
     reconErrorPercentage = zeros(totalSignatures, totalSamples);
     simError = zeros(totalSignatures, totalSamples);
        
     %% Examine individual signatures
     p = 1;
     for i = minSignatures : maxSignatures
         
        totalSignatures = i;
        outputFileTextSignatures = [outputFolder filesep 'text' filesep 'res_' inputFile '_signature_patterns_for_' num2str(totalSignatures) '_sigs.csv'];
        outputFileTextActivities = [outputFolder filesep 'text' filesep 'res_' inputFile '_signature_activities_for_' num2str(totalSignatures) '_sigs.csv'];
        outputFileNameFull = [outputFolder filesep 'full' filesep 'res_' inputFile '_full_signatures_' num2str(totalSignatures) '.mat'];
        outputFileNameSkinny = [outputFolder filesep 'skinny' filesep 'res_' inputFile '_skinny_signatures_' num2str(totalSignatures) '.mat'];
        
        %% Decipher the signatures of mutational processes from catalogues of mutations
        [input, allProcesses, allExposures, idx, processes, exposures, processStab, processStabAvg] = ...
            decipherMutationalProcesses(iterationsPerCore, totalSignatures, input, ...
                                        outputFileNameFull, outputFileNameSkinny, ...
                                        outputFileTextSignatures, outputFileTextActivities, job);
        
        recon = processes * exposures;
        avgReconstructionError(p) = norm(input.originalGenomes - recon, 'fro');
        avgReconstructionErrorPercentage(p) = norm(input.originalGenomes - recon, 'fro') ./ norm(input.originalGenomes);
        avgStability(p) = mean(processStabAvg);
        
        for j = 1 : totalSamples
            reconError(p, j) = norm(input.originalGenomes(:,j) - recon(:,j), 'fro');
            reconErrorPercentage(p, j) = norm(input.originalGenomes(:,j) - recon(:,j), 'fro') ./ norm(input.originalGenomes(:,j), 'fro');
            simError(p, j) = 1 - pdist([input.originalGenomes(:,j) recon(:,j)]', 'cosine');
        end
        p = p + 1;
        
     end

     %% Generate summary file
     outputFileNameSummary = [outputFolder filesep 'summary' filesep 'res_' inputFile '_summary.mat'];
     save(outputFileNameSummary, ...
            'avgReconstructionError', 'avgStability', 'input', 'maxSignatures', 'minSignatures', ...
            'reconError', 'reconErrorPercentage', 'simError', 'avgReconstructionErrorPercentage');
end
   