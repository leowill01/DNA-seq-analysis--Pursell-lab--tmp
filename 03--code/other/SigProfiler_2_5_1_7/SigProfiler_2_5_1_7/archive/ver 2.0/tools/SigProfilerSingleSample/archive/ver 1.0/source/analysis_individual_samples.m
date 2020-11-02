function analysis_individual_samples(signaturesSet, signaturesInCancerTypes, samplesForAnalysis, outputFolder, outputFile)

    %% Analysis of data
    % Loading samples
    input = load(samplesForAnalysis);
        totalSamples = size(input.originalGenomes, 2);
        accr_org = zeros(totalSamples, 1); kl_div_org = zeros(totalSamples, 1);  
        frob_rel_div_org = zeros(totalSamples, 1);  norm_one_dif_org = zeros(totalSamples, 1); 

    % Loading signatures
    newSignatures = load(signaturesSet);
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
                      
    parfor iSample = 1 : totalSamples

        %% Select cancer sample and set of signatures
        genome = input.originalGenomes(:,iSample);
        if ( strcmp(input.seqType(iSample), 'WGS') )
           allSignatures =  allGenomeSignatures;
        else
           allSignatures =  allExomeSignatures;
        end

        %% Identify signatures in cancer type and apply signature rules
        exposures = sigsInCanType.signaturesInCancerTypes(:, strcmpi(input.cancerType(iSample), sigsInCanType.cancerTypes))';                         
        exposures = check_signature_rules(exposures, sigNames, iSample, input.seqType, ...
                                          input.totalMutations, input.strandBias);

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
