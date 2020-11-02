function analysis_individual_samples(signaturesSet, signaturesInCancerTypes, samplesForAnalysis, outputFolder, outputFile, useRules, allowSigsEverywhere)

 %% Analysis of data
    % Loading samples
    input = load(samplesForAnalysis);
        totalSamples = size(input.originalGenomes, 2);
        accr_org = zeros(totalSamples, 1); kl_div_org = zeros(totalSamples, 1);  
        frob_rel_div_org = zeros(totalSamples, 1);  norm_one_dif_org = zeros(totalSamples, 1); 

    % Loading signatures
    newSignatures = load(signaturesSet);
        if ( exist('allowSigsEverywhere', 'var') == 1 )
            disp('Allow these signatures in all samples:');
            for i = 1 : length(allowSigsEverywhere)
                 disp(newSignatures.signatureNames{allowSigsEverywhere(i)});
            end
            idsToAdd = allowSigsEverywhere;
        else
            idsToAdd = 0;
        end
        allGenomeSignatures = newSignatures.genomeSignatures;   
        allExomeSignatures  = newSignatures.exomeSignatures;   
        sigNames = newSignatures.signatureNames;
        totalSignatures = size(allGenomeSignatures, 2);
        exposuresNew = zeros(totalSignatures, totalSamples);

    % Loading signatures in cancer types
    sigsInCanType = load(signaturesInCancerTypes);

    % Minization options
    saOptions = optimset( 'Display', 'off', 'TolFun', 1e-100, ...
                          'MaxFunEvals', Inf, 'MaxIter', 100000, ...
                          'Algorith', 'interior-point', 'FinDiffType', 'central', ...
                          'TolCon', 1e-100, 'TolX', 1e-100 );
                      
    parfor iSample = 1 : 10 %totalSamples

        %% Select cancer sample and set of signatures
        genome = input.originalGenomes(:,iSample);
        if ( strcmp(input.seqType(iSample), 'WGS') )
           allSignatures =  allGenomeSignatures;
        else
           allSignatures =  allExomeSignatures;
        end

        %% Identify signatures in cancer type and apply signature rules
        exposures = sigsInCanType.signaturesInCancerTypes(strcmpi(input.cancerType(iSample), sigsInCanType.cancerTypes), :)';    
        
        if ( useRules == 1 )
            exposures = check_signature_rules(exposures, sigNames, iSample, input.seqType, ...
                                              input.totalMutations, input.strandBias);
        end

        %% Remove all signatures in the sample: one by one
        exposures = remove_all_single_signatures(exposures, allSignatures, genome, saOptions, sigNames);

        %% Add all remaining signatures to the sample: one by one
        exposures = add_all_single_signatures(exposures, allSignatures, genome, saOptions, sigNames);

        %% Allow signatures 1 and 5 everywhere
        if ( idsToAdd ~= 0 )
            exposures(idsToAdd) = 1;
        end
        
        [accr_org(iSample), kl_div_org(iSample), frob_rel_div_org(iSample), norm_one_dif_org(iSample), exposures] = ...
                                                            eval_single_sample(exposures, allSignatures, genome, saOptions);
        exposuresNew(:, iSample) = exposures;

        %% Dispay summary output data
        disp(['Sample #' num2str(iSample) ': ' input.cancerType{iSample} ' ' input.sampleNames{iSample} ' with ' num2str(sum(exposuresNew(:, iSample)>0)) ' signatures and an accuracy of ' num2str(accr_org(iSample),'%.2f')])
    end

    if ( exist(outputFolder,'dir') == 0 )
       mkdir(outputFolder);
    end
    save([outputFolder filesep outputFile]);
end
