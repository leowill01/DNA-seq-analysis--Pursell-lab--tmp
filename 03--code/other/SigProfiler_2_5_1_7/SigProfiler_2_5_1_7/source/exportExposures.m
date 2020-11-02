function exportExposures(exposures, sampleNames, outputFile, dataSep)

    totalSignatures = size(exposures, 1);
    totalSamples    = size(exposures, 2);

    fileID = fopen(outputFile, 'w');
    fprintf(fileID, 'Sample Names');
    
    % Print signature names
    for iSig = 1 : totalSignatures
        fprintf(fileID, '%s%s', dataSep, ['Signature ' num2str(iSig)]);
    end
    fprintf(fileID, '\n');
    
    for iSample = 1 : totalSamples
        fprintf(fileID, '%s', sampleNames{iSample});
        
        for iSig = 1 : totalSignatures
            fprintf(fileID, '%s%f', dataSep, exposures(iSig,iSample));
        end
        fprintf(fileID, '\n');
    end
    fclose(fileID);

end