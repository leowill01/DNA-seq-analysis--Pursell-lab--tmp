function analysis_individual_samples(signaturesSet, signaturesInSamples, samplesForAnalysis, outputFolder, outputFile, useRules, allowSigsEverywhere, connected)

 %% Analysis of data
    % Loading samples
    input = load(samplesForAnalysis);
        totalSamples = size(input.originalGenomes, 2);
        accr_org = zeros(totalSamples, 1); kl_div_org = zeros(totalSamples, 1);  
        frob_rel_div_org = zeros(totalSamples, 1);  norm_one_dif_org = zeros(totalSamples, 1); 

    % Loading signatures
    newSignatures = load(signaturesSet);
        if ( (exist('allowSigsEverywhere', 'var') == 1) && ~isempty(allowSigsEverywhere) )
            disp('Allow these signatures in all samples:');
            for i = 1 : length(allowSigsEverywhere)
                 disp(newSignatures.signatureNames{allowSigsEverywhere(i)});
            end
            idsToAdd = allowSigsEverywhere;
        else
            idsToAdd = 0;
        end
        
        if ( (exist('connected', 'var') == 1) && ~isempty(connected) )
            disp('These mutational signatures are connected:');
            for i = 1 : length(connected)
                dispStr = ['Set ' num2str(i) ':' newSignatures.signatureNames{connected{i}(1)}];
                for j = 2 : length(connected{i})
                    dispStr = horzcat(dispStr, '; ', newSignatures.signatureNames{connected{i}(j)});
                end
                disp(dispStr);
            end
        end
        
        allGenomeSignatures = newSignatures.genomeSignatures;   
        allExomeSignatures  = newSignatures.exomeSignatures;   
        sigNames = newSignatures.signatureNames;
        totalSignatures = size(allGenomeSignatures, 2);
        exposuresNew = zeros(totalSignatures, totalSamples);

    % Loading signatures in samples
    sigsInCanType = load(signaturesInSamples);
    for i = 1 : size(sigsInCanType.sampleNames, 1)
       sigsInCanType.longSampleNames{i} = [sigsInCanType.sampleCancerTypes{i} '::' sigsInCanType.sampleNames{i}];
    end
    
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
        longSampleName = [input.cancerType{iSample} '::' input.sampleNames{iSample}];
        exposures = sigsInCanType.signaturesInSamples(:, strcmpi(longSampleName, sigsInCanType.longSampleNames)); 
        if ( isempty(exposures) )
            exposures = zeros(totalSignatures, 1);
        end
        
        if ( useRules == 1 )
            exposures = check_signature_rules(exposures, sigNames, iSample, input.seqType, ...
                                              input.totalMutations, input.strandBias);
        end

        %% Remove all signatures in the sample: one by one
        % Add signatures that are allowed everywhere
        if ( ~isempty(idsToAdd) )
            exposures(idsToAdd) = 1;
        end
        
        % Add connected signatures
        if (~isempty(connected))
            for iConnect = 1 : length(connected)
                if ( sum(exposures(connected{iConnect})) > 0 )
                    exposures(connected{iConnect}) = 1;
                end
            end
        end
        
        exposures = remove_all_single_signatures(exposures, allSignatures, genome, saOptions, sigNames);

        %% Add all remaining signatures to the sample: one by one
        % Add signatures that are allowed everywhere
        if ( ~isempty(idsToAdd) )
            exposures(idsToAdd) = 1;
        end
        % Add connected signatures
        if (~isempty(connected))
            for iConnect = 1 : length(connected)
                if ( sum(exposures(connected{iConnect})) > 0 )
                    exposures(connected{iConnect}) = 1;
                end
            end
        end
        exposures = add_all_single_signatures(exposures, allSignatures, genome, saOptions, sigNames);

        %% Add signatures that are allowed everywhere
        % Add signatures that are allowed everywhere
        if ( ~isempty(idsToAdd) )
            exposures(idsToAdd) = 1;
        end
        % Add connected signatures
        if (~isempty(connected))
            for iConnect = 1 : length(connected)
                if ( sum(exposures(connected{iConnect})) > 0 )
                    exposures(connected{iConnect}) = 1;
                end
            end
        end
        
        [exposures, accr_org(iSample), kl_div_org(iSample), frob_rel_div_org(iSample), norm_one_dif_org(iSample)] = ...
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






