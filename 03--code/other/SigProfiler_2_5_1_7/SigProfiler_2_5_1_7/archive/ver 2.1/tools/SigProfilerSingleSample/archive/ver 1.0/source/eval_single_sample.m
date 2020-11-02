function [accr, kl_div, frob_rel_div, norm_one_dif, exposuresSample] = eval_single_sample(exposures, allSignatures, genome, saOptions)

    exposuresSample = zeros(size(allSignatures, 2), 1);
    numberOfSignatures = sum(exposures>0);
    maxMutations = sum(genome);
    
    x0 = maxMutations * rand(numberOfSignatures, 1); 
    x0 = x0 / sum(x0) * maxMutations;
    A = ones(1, numberOfSignatures );
    lb = zeros( numberOfSignatures , 1);
    ub = maxMutations * ones( numberOfSignatures, 1);
    
    subSignatures = allSignatures(:,exposures>0);

    ObjectiveFunction = @(x) parameterized_objective2_custom(x, subSignatures, ...
                                                                      genome);

    [x, minFunction] = fmincon(ObjectiveFunction, x0, [], [], A, maxMutations, lb, ub, [], saOptions); % simulannealbnd(ObjectiveFunction, X0, lb, ub, saOptions);
    x = round(x);
    if ( sum(x) ~= maxMutations)
         [A B] = max(x);
         x(B) = x(B) + maxMutations - sum(x);
    end

    
    exposuresSample(exposures>0) = x;
    recon = allSignatures * exposuresSample; 
    accr = 1 - pdist([recon genome]', 'cosine');
    kl_div = KLDiv(genome', recon');
    frob_rel_div = norm(recon - genome, 'fro') ./ norm(genome, 'fro');
    norm_one_dif = norm(recon - genome, 1);
end