function strandBias = import_192_mutation_types_single_sample(inputFile, templateFolder)

    mutationTypeColumns = 3;
    input = readtable( inputFile,'ReadVariableNames',false);
            columnNames = input{1,:};
            input(1,:) = [];
    output = load([templateFolder filesep '192_template.mat']);
    output.originalGenomes = zeros(size(input,1), size(input,2) - mutationTypeColumns);

    sampleNames = cell(length(columnNames) - mutationTypeColumns, 1);
    cancerType = cell(length(columnNames)- mutationTypeColumns, 1);
    for i = (mutationTypeColumns + 1) : length(columnNames)
        splitCancerSample = strsplit(columnNames{i}, '::');
        cancerType(i - mutationTypeColumns) = splitCancerSample(1);
        sampleNames(i - mutationTypeColumns) = splitCancerSample(2);
    end

    output.cancerType = cancerType;
    output.sampleNames = sampleNames;
    indexId = zeros(size(input,1), 1);

    for i = 1 : size(input, 1)
        strandType = input{i, 1};
        mutType    = input{i, 2};
        mutSubType = input{i, 3};
        indexId(i) = find( strcmp(strandType, output.strandTypes) & ...
                           strcmp(mutType, output.types) & ...
                           strcmp(mutSubType, output.subtypes) );
    end

    for featureID = 1 : size(input,1)
        for sampleID = (mutationTypeColumns + 1) : length(columnNames)
            output.originalGenomes(indexId(featureID),sampleID - mutationTypeColumns) = str2double(cell2mat(input{featureID, sampleID}));
        end
    end
    
    strandBias = evaluate_all_TSB(output);
end