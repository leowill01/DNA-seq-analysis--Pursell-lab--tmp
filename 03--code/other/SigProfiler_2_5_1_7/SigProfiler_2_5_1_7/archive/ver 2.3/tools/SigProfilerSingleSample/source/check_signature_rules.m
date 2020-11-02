function signaturesInSample = check_signature_rules(signaturesInSample, signatures, sampleID, ...
                                                    seqType, totalMutations, strandBias, indels)
    totalSignatures = length(signaturesInSample);
    strand_bias_cutoff = 10^-2;
    
    if ( strcmp(seqType(sampleID), 'WGS') )
        pole_subs_cutoff = 10^5;
        msi_subs_cutoff = 10^4;
    elseif ( strcmp(seqType(sampleID), 'WES') )
        pole_subs_cutoff = 2 * 10^3;
        msi_subs_cutoff = 2 * 10^2;
    else
        error('Invalid type of sequencing data!')
    end
    
    
    for i = 1 : totalSignatures

        %% Transcriptional Strand Bias
        % [C>A transcribed] Check C>A transcriptional strand bias for signatures 4, 8, 24, 35, 42, and 51
        if ( strcmp('Signature Subs-04', signatures(i)) || ...
             strcmp('Signature Subs-08', signatures(i)) || ...
             strcmp('Signature Subs-24', signatures(i)) || ...
             strcmp('Signature Subs-35', signatures(i)) || ...
             strcmp('Signature Subs-42', signatures(i)) || ...
             strcmp('Signature Subs-51', signatures(i)) )

           if ( (strandBias.C_to_A_p(sampleID) > strand_bias_cutoff) || ...
                (strandBias.C_to_A_d(sampleID) ~= -1) ) % damage on G
               signaturesInSample(i) = 0;
           end
        end    
        
        % [C>A untranscribed] Check C>A transcriptional strand bias for signatures 50
        if ( strcmp('Signature Subs-50', signatures(i)) )

           if ( (strandBias.C_to_A_p(sampleID) > strand_bias_cutoff) || ...
                (strandBias.C_to_A_d(sampleID) ~= 1) ) % damage on C
               signaturesInSample(i) = 0;
           end
        end            
                
        % [C>T transcribed] Check C>T transcriptional strand bias for signatures 19, 23, 31, 32, 42, and 51
        if ( strcmp('Signature Subs-19', signatures(i)) || ...
             strcmp('Signature Subs-23', signatures(i)) || ...
             strcmp('Signature Subs-31', signatures(i)) || ...
             strcmp('Signature Subs-32', signatures(i)) || ...
             strcmp('Signature Subs-42', signatures(i)) ) 
         
            if ( (strandBias.C_to_T_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.C_to_T_d(sampleID) ~= -1) ) % damage on G
               signaturesInSample(i) = 0;
            end
        end     
        
        % [C>T untranscribed] Check C>T transcriptional strand bias for signatures 7a and 7b
        if ( strcmp('Signature Subs-07a', signatures(i)) || ...
             strcmp('Signature Subs-07b', signatures(i)) || ...
             strcmp('Signature Subs-51',  signatures(i)) )
         
            if ( (strandBias.C_to_T_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.C_to_T_d(sampleID) ~= 1) ) % damage on C
               signaturesInSample(i) = 0;
            end
        end     
        
        % [T>A transcribed] Check T>A transcriptional strand bias for signature 22, 25, and 27
        if ( strcmp('Signature Subs-22', signatures(i)) || ...
             strcmp('Signature Subs-25', signatures(i)) || ...
             strcmp('Signature Subs-27', signatures(i)) || ...
             strcmp('Signature Subs-51', signatures(i)) )
         
            if ( (strandBias.T_to_A_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.T_to_A_d(sampleID) ~= -1) ) % damage on A
               signaturesInSample(i) = 0;
            end
        end  
        
        % [T>A untranscribed] Check T>A transcriptional strand bias for signature 7c
        if ( strcmp('Signature Subs-07c', signatures(i)) || ...
             strcmp('Signature Subs-47',  signatures(i)) )
            
            if ( (strandBias.T_to_A_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.T_to_A_d(sampleID) ~= 1) ) % damage on T
               signaturesInSample(i) = 0;
            end
        end  
        
        % [T>C transcribed] Check T>C transcriptional strand bias for signatures 5, 12, and 16
        if ( strcmp('Signature Subs-05', signatures(i)) || ...
             strcmp('Signature Subs-12', signatures(i)) || ...
             strcmp('Signature Subs-16', signatures(i)) || ...
             strcmp('Signature Subs-46', signatures(i)) ) 
         
            if ( (strandBias.T_to_C_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.T_to_C_d(sampleID) ~= -1) ) % damage on A
               signaturesInSample(i) = 0;
            end
        end  
        
        % [T>C untranscribed] Check T>C transcriptional strand bias for signatures 7c, 7d, 21, 26, 33, and 57
        if ( strcmp('Signature Subs-07c', signatures(i)) || ...
             strcmp('Signature Subs-07d', signatures(i)) || ...
             strcmp('Signature Subs-21',  signatures(i)) || ...
             strcmp('Signature Subs-26',  signatures(i)) || ...
             strcmp('Signature Subs-33',  signatures(i)) || ...
             strcmp('Signature Subs-57',  signatures(i))  ) 
         
            if ( (strandBias.T_to_C_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.T_to_C_d(sampleID) ~= 1) ) % damage on T
               signaturesInSample(i) = 0;
            end
        end  
        
        % [T>G transcribed] Check T>G transcriptional strand bias for signatures 5, 12, and 16
        if ( strcmp('Signature Subs-47', signatures(i))  ) 
         
            if ( (strandBias.T_to_G_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.T_to_G_d(sampleID) ~= -1) ) % damage on A
               signaturesInSample(i) = 0;
            end
        end  
        
        % [T>G untranscribed] Check T>G transcriptional strand bias for signatures 7d, 21, 26, and 33
        if ( strcmp('Signature Subs-44', signatures(i)) || ...
             strcmp('Signature Subs-57', signatures(i)) ) 
         
            if ( (strandBias.T_to_G_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.T_to_G_d(sampleID) ~= 1) ) % damage on T
               signaturesInSample(i) = 0;
            end
        end  

     
                
        %% Check large numbers of mutations for POLE signatures
        if ( strcmp('Signature Subs-10a', signatures(i)) || ...
             strcmp('Signature Subs-10b', signatures(i)) )
            
           % Checking numbers of single base mutations
           if ( totalMutations(sampleID) < pole_subs_cutoff )
               signaturesInSample(i) = 0;
           end
        end             
        
        %% Check short mutations for MSI signatures
        if ( strcmp('Signature Subs-06', signatures(i)) || ...
             strcmp('Signature Subs-14', signatures(i)) || ...
             strcmp('Signature Subs-15', signatures(i)) || ...
             strcmp('Signature Subs-20', signatures(i)) || ...
             strcmp('Signature Subs-21', signatures(i)) || ...
             strcmp('Signature Subs-26', signatures(i)) || ...
             strcmp('Signature Subs-60', signatures(i)) )

            % Checking numbers of subs
            if ( totalMutations(sampleID) < msi_subs_cutoff )
                signaturesInSample(i) = 0;
            end
        end
              
    end
end