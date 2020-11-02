function signaturesInSample = check_signature_rules(signaturesInSample, signatures, sampleID, ...
                                                    seqType, totalMutations, strandBias, indels)
    totalSignatures = length(signaturesInSample);
    strand_bias_cutoff = 10^-2;
    
    if ( strcmp(seqType(sampleID), 'WGS') )
        short_indels_cutoff = 5000;
        long_indels_cutoff = 100;
        pole_subs_cutoff = 10^5;
        msi_subs_cutoff = 10^4;
    elseif ( strcmp(seqType(sampleID), 'WES') )
        short_indels_cutoff = 100;
        long_indels_cutoff = 2;
        pole_subs_cutoff = 2 * 10^3;
        msi_subs_cutoff = 2 * 10^2;
    else
        error('Invalid type of sequencing data!')
    end
    
    for i = 1 : totalSignatures

        
        %% Check T>C at ATN transcriptional strand bias for signature 16
        if ( strcmp('Signature Subs-16', signatures(i)) ) 
            if ( (strandBias.T_to_C_ATN_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.T_to_C_ATN_d(sampleID) ~= -1) ) % damage on A
               signaturesInSample(i) = 0;
            end
        end  
        
        %% Check T>C transcriptional strand bias for signature 12
        if ( strcmp('Signature Subs-12', signatures(i)) ) 
            if ( (strandBias.T_to_C_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.T_to_C_d(sampleID) ~= -1) ) % damage on A
               signaturesInSample(i) = 0;
            end
        end  
        
        %% Check T>C transcriptional strand bias for signature 7d and 33
        if ( strcmp('Signature Subs-07d', signatures(i)) || ...
             strcmp('Signature Subs-33', signatures(i)) ) 
            if ( (strandBias.T_to_C_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.T_to_C_d(sampleID) ~= 1) ) % damage on T
               signaturesInSample(i) = 0;
            end
        end  
        
        %% Check C>T transcriptional strand bias for signature 11, 23, 31, 32, 35, and 42
        if ( strcmp('Signature Subs-11', signatures(i)) || ...
             strcmp('Signature Subs-23', signatures(i)) || ...
             strcmp('Signature Subs-31', signatures(i)) || ...
             strcmp('Signature Subs-32', signatures(i)) || ...
             strcmp('Signature Subs-35', signatures(i)) || ... 
             strcmp('Signature Subs-42', signatures(i)) ) 
            if ( (strandBias.C_to_T_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.C_to_T_d(sampleID) ~= -1) ) % damage on G
               signaturesInSample(i) = 0;
            end
        end     
        
        %% Check C>T transcriptional strand bias for signatures 7a and 7b
        if ( strcmp('Signature Subs-07a', signatures(i)) || ...
             strcmp('Signature Subs-07b', signatures(i)) ) 
            if ( (strandBias.C_to_T_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.C_to_T_d(sampleID) ~= 1) ) % damage on C
               signaturesInSample(i) = 0;
            end
        end     
        
        %% Check T>A transcriptional strand bias for signature 22, 25, 27
        if ( strcmp('Signature Subs-22', signatures(i)) || ...
             strcmp('Signature Subs-25', signatures(i)) || ...
             strcmp('Signature Subs-27', signatures(i)))
            if ( (strandBias.T_to_A_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.T_to_A_d(sampleID) ~= -1) ) % damage on A
               signaturesInSample(i) = 0;
            end
        end  
        
        %% Check T>A transcriptional strand bias for signature 7c
        if ( strcmp('Signature Subs-07c', signatures(i)) )
            if ( (strandBias.T_to_A_p(sampleID) > strand_bias_cutoff) || ...
                 (strandBias.T_to_A_d(sampleID) ~= 1) ) % damage on T
               signaturesInSample(i) = 0;
            end
        end     
        
        %% Check C>A transcriptional strand bias for signatures 4, 8, 24, and 29
        if ( strcmp('Signature Subs-04', signatures(i)) || ...
             strcmp('Signature Subs-08', signatures(i)) || ...
             strcmp('Signature Subs-24', signatures(i)) || ...
             strcmp('Signature Subs-29', signatures(i)) )
           
           % Checking for significant strand bias and correct direction
           if ( (strandBias.C_to_A_p(sampleID) > strand_bias_cutoff) || ...
                (strandBias.C_to_A_d(sampleID) ~= -1) ) % damage on G
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
        
        %% Check long indels for HRD signature
        if ( strcmp('Signature Subs-03', signatures(i)) )
            
           % Checking numbers of long indels
           if exist('indels','var')
               if ( sum(indels.originalGenomes([100 98], sampleID)) < long_indels_cutoff )
                   signaturesInSample(i) = 0;
               end
           end
        end      
        
        %% Check short indels for MSI signatures
        if ( strcmp('Signature Subs-06', signatures(i)) || ...
             strcmp('Signature Subs-14', signatures(i)) || ...
             strcmp('Signature Subs-15', signatures(i)) || ...
             strcmp('Signature Subs-20', signatures(i)) || ...
             strcmp('Signature Subs-21', signatures(i)) || ...
             strcmp('Signature Subs-26', signatures(i)) )
            
            % Checking numbers of short indels
            if exist('indels','var')
                if ( sum(indels.originalGenomes([99 97],sampleID)) < short_indels_cutoff )
                    signaturesInSample(i) = 0;
                end
            end
            
            % Checking numbers of subs
            if ( totalMutations(sampleID) < msi_subs_cutoff )
                signaturesInSample(i) = 0;
            end
        end
              
    end
end