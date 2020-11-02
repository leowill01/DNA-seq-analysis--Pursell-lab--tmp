function [accr, kl_div, frob_rel_div, norm_one_dif, exposuresSample, fID] = remove_one_sig(exposures, allSignatures, genome, saOptions)
    numberOfSignatures = sum(exposures>0);
    sig_IDs = find(exposures>0);
    totalSignatures = size(allSignatures, 2);
    
    if ( numberOfSignatures > 1 )
        exposuresSample_wo = zeros(totalSignatures, numberOfSignatures);
        accr_wo = zeros(numberOfSignatures, 1); kl_div_wo = zeros(numberOfSignatures, 1);
        frob_rel_div_wo = zeros(numberOfSignatures, 1); norm_one_dif_wo = zeros(numberOfSignatures, 1);
        for j = 1 : numberOfSignatures
            exposuresRemoveOne = exposures;
            exposuresRemoveOne(sig_IDs(j)) = 0; 
            [exposuresSample_wo(:,j), accr_wo(j), kl_div_wo(j), frob_rel_div_wo(j), norm_one_dif_wo(j)] = eval_single_sample(exposuresRemoveOne, allSignatures, genome, saOptions);
        end
        [fVal, fID] = min(harmmean([1-accr_wo frob_rel_div_wo],2));
        accr = accr_wo(fID);
        kl_div = kl_div_wo(fID);
        frob_rel_div = frob_rel_div_wo(fID);
        norm_one_dif = norm_one_dif_wo(fID);
        exposuresSample = exposuresSample_wo(:,fID);
        fID = sig_IDs(fID);
    else
        [exposuresSample, accr, kl_div, frob_rel_div, norm_one_dif] = eval_single_sample(exposures, allSignatures, genome, saOptions);
        fID = 0;
    end
    
end