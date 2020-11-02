function exportSignaturesPatterns(processes, types, subtypes, ...
                                  outputFile, dataSep)
    
    totalSignatures = size(processes, 2);
    totalMutationTypes =  length(subtypes);

    fileID = fopen(outputFile, 'w');
    fprintf(fileID, 'Mutation Type%sMutation Subtype', dataSep);
    
    % Print signature names
    for iSig = 1 : totalSignatures
        fprintf(fileID, '%s%s', dataSep, ['Signature ' num2str(iSig)]);
    end
    fprintf(fileID, '\n');
    

    for iType = 1 : totalMutationTypes 
        fprintf(fileID, '%s', types{iType});
        fprintf(fileID, '%s%s', dataSep, subtypes{iType});
        for iSig = 1 : totalSignatures
            fprintf(fileID, '%s%.10e', dataSep, processes(iType, iSig));
        end
                fprintf(fileID, '\n');
    end

    fclose(fileID);

end