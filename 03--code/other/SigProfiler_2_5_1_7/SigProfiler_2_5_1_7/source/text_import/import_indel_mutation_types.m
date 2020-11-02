function import_indel_mutation_types(inputFile, templateFolder, outputFile)

mutationTypeColumns = 4;
input = readtable( inputFile,'ReadVariableNames',false);
        columnNames = input{1,:};
input(1,:) = [];
output = load([templateFolder filesep 'indel_template.mat']);
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
    mutLabel   = [input{i, :}{1} '_' input{i, :}{2} '_' input{i, :}{3} '_' input{i, :}{4}];
    indexId(i) = find(strcmp(mutLabel, output.labels));
end

for featureID = 1 : size(input,1)
    for sampleID = (mutationTypeColumns + 1) : length(columnNames)
        output.originalGenomes(indexId(featureID),sampleID - mutationTypeColumns) = str2double(cell2mat(input{featureID, sampleID}));
    end
end

%% Save output file
output.sampleNames = output.sampleNames(sum(output.originalGenomes)>0);
output.cancerType  = output.cancerType(sum(output.originalGenomes)>0);
output.originalGenomes = output.originalGenomes(:,sum(output.originalGenomes)>0);

save(outputFile, '-struct', 'output');
end