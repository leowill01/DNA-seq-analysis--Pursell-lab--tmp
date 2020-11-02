%% The parameterization of objective function

function y = parameterized_objective2_custom(x, signatures, samples)

      rec = signatures * x; %% reconstructed genome using signatures and their activities
      
      y = norm(samples - rec, 2);        %% L^2-Norm 
%     y = KLDiv(samples', rec');         %% Kullback?Leibler divergence
%     y = norm(samples - rec, 'fro');    %% Frobenius Norm

end