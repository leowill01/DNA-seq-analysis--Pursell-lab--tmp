%% This function plots mutational samples based on types and subtypes
function plotSamplesWithStrandBias( samples, types, subtypes, strandTypes, sampleNames)
%
% Ludmil B. Alexandrov
% Cancer Genome Project
% Wellcome Trust Sanger Institute
% la2@sanger.ac.uk
%
% This software and its documentation are copyright 2012 by the
% Wellcome Trust Sanger Institute/Genome Research Limited. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever. 
% Neither the Wellcome Trust Sanger Institute nor Genome Research Limited 
% is responsible for its use, misuse, or functionality.
%
% Defining constants and initial sorting of types and subtypes
%

  originalSamples = samples;
  screen_size = get(0, 'ScreenSize');
  totalProcesses = size(samples, 2);
  totalMutationTypes = size(samples, 1);
  topElements = 5;
  [sortedType sortedIndex] = sort(types);
  uniqueTypes = unique(sortedType);
  sortedSubType = subtypes(sortedIndex);
  sortedSubSubType = strandTypes(sortedIndex);
  generateNames = 0;
  if ( exist('sampleNames', 'var') == 0 )
      generateNames = 1;
  end
  
  errors = eps*ones(totalMutationTypes, totalProcesses);
  
  for i = 1 : totalProcesses
      samples(:, i) = samples(sortedIndex, i);
      if ( generateNames == 1 )
          sampleNames{i} = ['Sample ' num2str(i) ' of ' num2str(totalProcesses)];
      end
  end
  
  edgeColor = [ 0.70 0.70 0.70 ];
  colors = zeros(10, 3);
  colors(1, :) = [ 0.68 0.92 1.00 ];
  colors(2, :) = [ 0.00 0.00 0.00 ];
  colors(3, :) = [ 0.93 0.84 0.84 ];
  colors(4, :) = [ 0.80 0.80 0.80 ];
  colors(5, :) = [ 0.76 0.87 0.78 ];
  colors(6, :) = [ 1.00 0.60 0.78 ];

  
  if ( length(uniqueTypes) == 6)
      subplotLeft = 2;
      subplotRight = 3;
  end
  
  % Plot mutational samples and their features
  for i = 1 : totalProcesses
      f1 = figure('XVisual','','InvertHardcopy','off','Color',[1 1 1]);
      set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
      [sortedProcess sortedProcessIndex] = sort(samples(:, i)); % sortedProcess is unused
      minTopElementValue = samples(sortedProcessIndex(totalMutationTypes-topElements), i);
      
      for j = 1 : length(uniqueTypes)
          subplot1 = subplot(subplotLeft, subplotRight, j,'LineWidth',1,'FontWeight','bold', 'FontSize',20);
          indxSubType = find ( strcmp(sortedType, uniqueTypes{j}) );
          [doubleSortedSubType doubleSortedSubTypeIndex] = sort(sortedSubType(indxSubType));
          sortedSubSubTypes = sortedSubSubType( indxSubType(doubleSortedSubTypeIndex));
          indxSubTypeSorted = indxSubType( doubleSortedSubTypeIndex );
          
          
          processesFixed = zeros( round(length(indxSubTypeSorted)/2), 2);
          processesFixed(:, 1) = samples( indxSubTypeSorted(strcmp(sortedSubSubType(indxSubTypeSorted), 'T')), i);
          processesFixed(:, 2) = samples( indxSubTypeSorted(strcmp(sortedSubSubType(indxSubTypeSorted), 'U')), i);
          
          processesErrors = zeros( round(length(indxSubTypeSorted)/2), 2);
          processesErrors(:, 1) = errors( indxSubTypeSorted(strcmp(sortedSubSubType(indxSubTypeSorted), 'T')), i);
          processesErrors(:, 2) = errors( indxSubTypeSorted(strcmp(sortedSubSubType(indxSubTypeSorted), 'U')), i);
          
          percentT = round(10000 * (sum(processesFixed(:,1))/sum(sum(processesFixed)) ) ) /100;
          precentU = 100 - percentT;
          
          % mutationByProcess = round(processesFixed *  sum(exposure(i,:)));
          
          handles.bars =  bar( processesFixed );%, 'FaceColor', colors(j, :), 'EdgeColor', edgeColor, 'BarWidth', 0.4);
          hold on;
          
          x =get(get(handles.bars(1),'children'), 'xdata');
          x = mean(x([1 3],:));
          handles.errors(i) = errorbar(x,  processesFixed(:,1),  processesErrors(:, 1), 'k', 'linestyle', 'none', 'linewidth', 0.5);
          
          x =get(get(handles.bars(2),'children'), 'xdata');
          x = mean(x([1 3],:));
          handles.errors(i) = errorbar(x,  processesFixed(:,2),  processesErrors(:, 2), 'k', 'linestyle', 'none', 'linewidth', 0.5);
          
          hold off;
          
          set(handles.bars(1),'DisplayName','Transcribed');
          set(handles.bars(2),'DisplayName','Untranscribed');

          axis([ 0.5 (round(length(indxSubTypeSorted)/2) +0.5) 0 (max(samples(:, i)+0.002)+5) ]);
          title( { uniqueTypes{j} [ThousandSep(sum(samples(indxSubTypeSorted, i))) ' mutations (T:' num2str(percentT) '% / U:' num2str(precentU) '%)' ]},'FontWeight','bold','FontSize',16 );
          
          
          set(subplot1,'XTick', 1:round(length(indxSubTypeSorted)/2) );
          set(subplot1,'XTickLabel',doubleSortedSubType((strcmp(sortedSubSubType(indxSubTypeSorted), 'T'))) );

          for iColorIndex = 1 : round(totalMutationTypes/12)
              strandBiasColor{iColorIndex} = 'k';     
          end
          
          rotateticklabelWithColor(subplot1, 90, strandBiasColor);
      end
      
      annotation(f1,'textbox',...
    [0.003 0.96 0.11 0.043],...
    'String',{sampleNames{i}, [ThousandSep(sum(samples(:,i))) ' mutations'], ['T:' ThousandSep(sum(sum(originalSamples(strcmp(strandTypes,'T'),i)))) ' / U:' ThousandSep(sum(sum(originalSamples(strcmp(strandTypes,'U'),i)))) ]},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor','none',...
    'Interpreter', 'none', ...
    'Color',[1 0 0]);
    
  end
  
  
end