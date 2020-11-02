function [exposuresOutput, accr, kl_div, frob_rel_div, norm_one_dif] = remove_all_single_signatures(exposures_init, allSignatures, genome, saOptions, sigNames)

    %% Check initial sample
    if ( sum(exposures_init) > 0)
        [exposures, accr_first, kl_div_first, frob_rel_div_first, norm_one_dif_first] = eval_single_sample(exposures_init, allSignatures, genome, saOptions);
        exposuresOutput = exposures;     
        accr = accr_first; kl_div = kl_div_first; 
        frob_rel_div = frob_rel_div_first; norm_one_dif =  norm_one_dif_first;

        %% Removing singature one by one
        numSig = sum(exposures>0);
        for j = 1 : ( numSig - 1)
            [accr_temp, kl_div_temp, frob_rel_div_temp, norm_one_dif_temp, exposuresSample_temp, fID] = ...
                                                    remove_one_sig(exposuresOutput, allSignatures, genome, saOptions);
            if ( (accr_first - accr_temp < 0.01) )
                exposuresOutput = exposuresSample_temp;     
                accr = accr_temp; kl_div = kl_div_temp; 
                frob_rel_div = frob_rel_div_temp; norm_one_dif =  norm_one_dif_temp;
            else 
                break;
            end
        end
    else
      exposuresOutput = exposures_init;
      accr = 0;
      kl_div = 0;
      frob_rel_div = 0;
      norm_one_dif = 0;
    end
    
end