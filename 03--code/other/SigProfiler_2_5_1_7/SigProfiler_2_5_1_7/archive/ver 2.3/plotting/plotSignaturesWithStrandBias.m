%% This function plots mutational processes based on types and subtypes
function plotSignaturesWithStrandBias( processes, allProcesses, idx, input )

% Defining constants and initial sorting of types and subtypes
%
  type = input.types;
  subtype = input.subtypes;
  subsubtype = input.strandTypes;
  
  screen_size = get(0, 'ScreenSize');
  totalProcesses = size(processes, 2);
  totalMutationTypes = size(processes, 1);
  topElements = 5;
  [sortedType sortedIndex] = sort(type);
  uniqueTypes = unique(sortedType);
  sortedSubType = subtype(sortedIndex);
  sortedSubSubType = subsubtype(sortedIndex);
  generateNames = 0;
  if ( exist('processNames', 'var') == 0 )
      generateNames = 1;
  end
  
  if ( isfield(input, 'cancerType') == 0 )
    cancerType = 'Unknown';
  else
    cancerType = input.cancerType;  
  end
  
  errors = zeros(totalMutationTypes, totalProcesses);
  for i = 1 : totalProcesses
      errors(:,i) = std(allProcesses(:, idx==i), [], 2);
  end
  
  for i = 1 : totalProcesses
      processes(:, i) = processes(sortedIndex, i);
      if ( generateNames == 1 )
          processNames{i} = ['Signature ' num2str(i) ' of ' num2str(totalProcesses)];
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
  
  % Plot mutational processes and their features
  for i = 1 : totalProcesses
      f1 = figure('InvertHardcopy','off','Color',[1 1 1]);
      set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
      [sortedProcess sortedProcessIndex] = sort(processes(:, i)); % sortedProcess is unused
      minTopElementValue = processes(sortedProcessIndex(totalMutationTypes-topElements), i);
      
      for j = 1 : length(uniqueTypes)
          subplot1 = subplot(subplotLeft, subplotRight, j,'LineWidth',1,'FontWeight','bold', 'FontSize',20);
          indxSubType = find ( strcmp(sortedType, uniqueTypes{j}) );
          [doubleSortedSubType doubleSortedSubTypeIndex] = sort(sortedSubType(indxSubType));
          sortedSubSubTypes = sortedSubSubType( indxSubType(doubleSortedSubTypeIndex));
          indxSubTypeSorted = indxSubType( doubleSortedSubTypeIndex );
          
          
          processesFixed = zeros( round(length(indxSubTypeSorted)/2), 2);
          processesFixed(:, 1) = processes( indxSubTypeSorted(strcmp(sortedSubSubType(indxSubTypeSorted), 'T')), i);
          processesFixed(:, 2) = processes( indxSubTypeSorted(strcmp(sortedSubSubType(indxSubTypeSorted), 'U')), i);
          
          processesErrors = zeros( round(length(indxSubTypeSorted)/2), 2);
          processesErrors(:, 1) = errors( indxSubTypeSorted(strcmp(sortedSubSubType(indxSubTypeSorted), 'T')), i);
          processesErrors(:, 2) = errors( indxSubTypeSorted(strcmp(sortedSubSubType(indxSubTypeSorted), 'U')), i);
          
          percentT = round(10000 * sum( processesFixed(:, 1))) /100;
          precentU = round(10000 * sum(processes(indxSubTypeSorted, i)))/100 - percentT;
          
          % mutationByProcess = round(processesFixed *  sum(exposure(i,:)));
          
          handles.bars =  bar( processesFixed );%, 'FaceColor', colors(j, :), 'EdgeColor', edgeColor, 'BarWidth', 0.4);
          hold on;
          
          if ~verLessThan('matlab', '8.4')
               x =  handles.bars(1).XData + handles.bars(1).XOffset;
          else
               x =get(get(handles.bars(1),'children'), 'xdata');
               x = mean(x([1 3],:));
          end
          
          handles.errors(i) = errorbar(x,  processesFixed(:,1),  processesErrors(:, 1), 'k', 'linestyle', 'none', 'linewidth', 0.5);
          
          if ~verLessThan('matlab', '8.4')
               x =  handles.bars(2).XData + handles.bars(2).XOffset;
          else
               x =get(get(handles.bars(1),'children'), 'xdata');
               x = mean(x([1 3],:));
          end
          
          handles.errors(i) = errorbar(x,  processesFixed(:,2),  processesErrors(:, 2), 'k', 'linestyle', 'none', 'linewidth', 0.5);
          
          hold off;
          
          set(handles.bars(1),'DisplayName','Transcribed');
          set(handles.bars(2),'DisplayName','Untranscribed');

          axis([ 0.5 (round(length(indxSubTypeSorted)/2) +0.5) max(min(processes(:, i)-0.002),0) min(max(processes(:, i)+0.002),1.0001) ]);
          title( { uniqueTypes{j} [num2str( round(10000 * sum(processes(indxSubTypeSorted, i)))/100) '% (T:' num2str(percentT) '% / U:' num2str(precentU) '%) of all mutations' ]},'FontWeight','bold','FontSize',16 );
          
          
          set(subplot1,'XTick', 1:round(length(indxSubTypeSorted)/2) );
          set(subplot1,'XTickLabel',doubleSortedSubType((strcmp(sortedSubSubType(indxSubTypeSorted), 'T'))) );

          for iColorIndex = 1 : round(totalMutationTypes/12)
              strandBiasColor{iColorIndex} = 'k';

                  
         end
          
          rotateticklabelWithColor(subplot1, 90, strandBiasColor);
      end
      annotation(f1,'textarrow',[0.0585937500000001 0.12890625],...
    [0.535290355247701 0.582683720176611],'TextEdgeColor','none',...
    'TextRotation',90,...
    'FontWeight','bold',...
    'FontSize',20,...
    'String', processNames{i},...
    'HeadStyle','none',...
    'LineStyle','none');
    legend1 = legend(subplot1,{'Transcribed', 'Untranscribed'});
    set(legend1,...
        'Position',[0.00179565853542754 0.00716395967609684 0.115391841464572 0.141894725532195],...
        'FontWeight','bold',...
        'FontSize',15);
    
    annotation(f1,'textbox',...
    [0.003 0.96 0.11 0.043],...
    'String',{[num2str(size(input.sampleNames, 1)) ' ' cancerType ' samples'], [ThousandSep(sum(sum(input.originalGenomes))) ' mutations']},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor','none',...
    'Color',[1 0 0]);
  end
  
  
end