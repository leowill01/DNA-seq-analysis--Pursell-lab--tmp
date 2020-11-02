function analysisSignaturesTextIO(inputCSVFolder, outputFolder, inputCSVFile, inputCSVFileType, minSignatures, maxSignatures, iterationsPerCore, job) 

     %% Generate MATLAB input files
     outputSubFolder = strsplit(outputFolder, filesep);
     if ( strcmp(inputCSVFileType, '96') == 1 )
         inputTempFileName = outputSubFolder{end};
         import_96_mutation_types([inputCSVFolder filesep inputCSVFile], 'source/text_import/templates/', ['input' filesep inputTempFileName '.mat']);
     elseif ( strcmp(inputCSVFileType, '192') == 1 )
         inputTempFileName = outputSubFolder{end};
         import_192_mutation_types([inputCSVFolder filesep inputCSVFile], 'source/text_import/templates/', ['input' filesep inputTempFileName '.mat']);
     elseif ( strcmp(inputCSVFileType, '1536') == 1 )   
         inputTempFileName = outputSubFolder{end};
         import_1536_mutation_types([inputCSVFolder filesep inputCSVFile], 'source/text_import/templates/', ['input' filesep inputTempFileName '.mat']);
     elseif ( strcmp(inputCSVFileType, 'dinucs') == 1 )
         inputTempFileName = outputSubFolder{end};
         import_dinuc_mutation_types([inputCSVFolder filesep inputCSVFile], 'source/text_import/templates/', ['input' filesep inputTempFileName '.mat']);
     elseif ( strcmp(inputCSVFileType, 'indels') == 1 )
         inputTempFileName = outputSubFolder{end};
         import_indel_mutation_types([inputCSVFolder filesep inputCSVFile], 'source/text_import/templates/', ['input' filesep inputTempFileName '.mat']);
     else
         error( 'analysisSignaturesTextIO: Please specify a valid inputCSVFileType! Valid formats include: ''96'', ''192'', ''1536'', ''dinucs'', or ''indels''.');
     end
     
    inputFolder  = 'input';
    inputFile    = inputTempFileName; % please do not use an extension 
    
    analysisSignatures(inputFolder, outputFolder, inputFile, ... % input folder; output folder; input file
                       minSignatures, maxSignatures, iterationsPerCore, ... % start number of signatures; end number of signatures; total iterations per core (please make at least 1,000 iterations)
                       job); % parallel job variable
    delete([inputFolder filesep inputFile '.mat']);
end
   