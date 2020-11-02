function import_96_mutation_types_single_sample(inputFile, templateFolder, outputFile, seqType, samplesForAnalysis192)

mutationTypeColumns = 2;
input = readtable(inputFile, 'ReadVariableNames', false);
        columnNames = input{1,:};
        input(1,:) = [];
output = load([templateFolder filesep '96_template.mat']);
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
    mutType    = input{i, 1};
    mutSubType = input{i, 2};
    indexId(i) = find(strcmp(mutType, output.types) & strcmp(mutSubType, output.subtypes));
end

for featureID = 1 : size(input,1)
    for sampleID = (mutationTypeColumns + 1) : length(columnNames)
        output.originalGenomes(indexId(featureID),sampleID - mutationTypeColumns) = str2double(cell2mat(input{featureID, sampleID}));
    end
end

output.totalMutations = sum(output.originalGenomes)';
output.seqType = repmat({seqType}, length(output.sampleNames), 1);

%% Generating strand bias information
output.strandBias = import_192_mutation_types_single_sample(samplesForAnalysis192, templateFolder);

%% Save output file
save(outputFile, '-struct', 'output');
end