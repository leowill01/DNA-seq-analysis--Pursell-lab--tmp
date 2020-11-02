function [exposuresOutput, accr, kl_div, frob_rel_div, norm_one_dif] = add_all_single_signatures(exposures_init, allSignatures, genome, saOptions, sigNames)

    %% Check initial sample
    [exposuresOutput, accr, kl_div, frob_rel_div, norm_one_dif] = eval_single_sample(exposures_init, allSignatures, genome, saOptions);
    
    %% Removing singature one by one
    numSig = length(exposuresOutput) - sum(exposuresOutput>0);
    for j = 1 : numSig
        [accr_temp, kl_div_temp, frob_rel_div_temp, norm_one_dif_temp, exposuresSample_temp, fID] = ...
                                                add_one_sig(exposuresOutput, allSignatures, genome, saOptions);
               
        if ( (accr_temp - accr > 0.05) )
            exposuresOutput = exposuresSample_temp;     
            accr = accr_temp; kl_div = kl_div_temp; 
            frob_rel_div = frob_rel_div_temp; norm_one_dif =  norm_one_dif_temp;
        else 
            break;
        end
    end
end