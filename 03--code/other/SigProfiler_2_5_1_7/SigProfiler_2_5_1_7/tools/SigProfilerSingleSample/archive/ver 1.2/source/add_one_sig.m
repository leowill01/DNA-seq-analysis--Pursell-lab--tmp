function [accr, kl_div, frob_rel_div, norm_one_dif, exposuresSample, fID] = add_one_sig(exposures, allSignatures, genome, saOptions)
    numberOfSignatures = length(exposures) - sum(exposures>0);
    sig_IDs = find(exposures==0);
    totalSignatures = size(allSignatures, 2);
    
    if ( numberOfSignatures > 1 )
        exposuresSample_with = zeros(totalSignatures, numberOfSignatures);
        accr_with = zeros(numberOfSignatures, 1); kl_div_with = zeros(numberOfSignatures, 1);
        frob_rel_div_with = zeros(numberOfSignatures, 1); norm_one_dif_with = zeros(numberOfSignatures, 1);
        for j = 1 : numberOfSignatures
            exposuresAddOne = exposures;
            exposuresAddOne(sig_IDs(j)) = 1; 
            [accr_with(j), kl_div_with(j), frob_rel_div_with(j), norm_one_dif_with(j), exposuresSample_with(:,j)] = eval_single_sample(exposuresAddOne, allSignatures, genome, saOptions);
        end
        [fVal, fID] = min(harmmean([1-accr_with frob_rel_div_with],2));
        accr = accr_with(fID);
        kl_div = kl_div_with(fID);
        frob_rel_div = frob_rel_div_with(fID);
        norm_one_dif = norm_one_dif_with(fID);
        exposuresSample = exposuresSample_with(:,fID);
        fID = sig_IDs(fID);
    else
        [accr, kl_div, frob_rel_div, norm_one_dif, exposuresSample] = eval_single_sample(exposures, allSignatures, genome, saOptions);
        fID = 0;
    end
    
end