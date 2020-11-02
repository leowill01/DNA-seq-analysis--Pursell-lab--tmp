function plotSignatures(processes, input, allProcesses, idx, processStabAvg, processNames)

  if ( isfield(input, 'cancerType') == 0 )
    cancerType = 'Unknown';
  else
    cancerType = input.cancerType;  
  end
  
  type = input.types;
  subtype = input.subtypes;
  
  screen_size = get(0, 'ScreenSize');
  totalProcesses = size(processes, 2);
  totalMutationTypes = size(processes, 1);
  topElements = 16;
  totalDisplayElements = 10;
  titleFont = 20;
  fontLabels = 4;
  [sortedType sortedIndex] = sort(type);
  uniqueTypes = unique(sortedType);
  sortedSubType = subtype(sortedIndex);
  generateNames = 0;
  if ( exist('processNames', 'var') == 0 )
      generateNames = 1;
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
  
  if ( size(type, 1) < 200 )
      fontLabels = 12;
      totalDisplayElements = 3;
      titleFont = 15;
  end
  
  
  if ( size(type, 1) == 96 )
      fontLabels = 14;
      totalDisplayElements = 3;
  end
  
  edgeColor = [ 0.70 0.70 0.70 ];
  colors = zeros(10, 3);
  colors(1, :) = [ 0.68 0.92 1.00 ];
  colors(2, :) = [ 0.85 0.70 1.00 ];
  colors(3, :) = [ 0.93 0.84 0.84 ];
  colors(4, :) = [ 0.80 0.80 0.80 ];
  colors(5, :) = [ 0.76 0.87 0.78 ];
  colors(6, :) = [ 1.00 0.60 0.78 ];
  colors(7, :) = [ 0.00 0.00 1.00 ];
  
  texboxesLocations(1, :) = [0.281770833333334 0.895427438057103 0.07109375 0.0238805970149254];
  texboxesLocations(2, :) = [0.5665625 0.897053769131643 0.07109375 0.0238805970149254];
  texboxesLocations(3, :) = [0.845520833333333 0.897692684910927 0.07109375 0.0238805970149254];
  texboxesLocations(4, :) = [0.285625 0.423346121503376 0.07109375 0.0238805970149254];
  texboxesLocations(5, :) = [0.564270833333334 0.42237806729234 0.07109375 0.0238805970149254];
  texboxesLocations(6, :) = [0.8465625 0.425282229925447 0.07109375 0.0238805970149254];
 

  
  if ( length(uniqueTypes) == 6)
      subplotLeft = 2;
      subplotRight = 3;
  end
  
 if ( length(uniqueTypes) == 7)
      subplotLeft = 2;
      subplotRight = 4;
      colors(5, :) = [ 0.80 0.80 0.80 ];
      colors(6, :) = [ 0.76 0.87 0.78 ];
      colors(7, :) = [ 1.00 0.60 0.78 ];
      colors(4, :) = [ 0.00 0.30 0.50 ];
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
          indxSubTypeSorted = indxSubType(doubleSortedSubTypeIndex);

          
          errorsPlotType = errors(indxSubTypeSorted, i);
          
          handles.bars = bar( 1:length(indxSubTypeSorted), processes(indxSubTypeSorted, i), 'FaceColor', colors(j, :), 'EdgeColor', edgeColor, 'BarWidth', 0.4);
   
          hold on;
          plot(0.5:0.5:(length(indxSubTypeSorted) +0.5), minTopElementValue * ones(1,length(0.5:0.5:(length(indxSubTypeSorted) +0.5))), 'LineStyle','-.','Color',[1 0 0]);
     
          title( [ uniqueTypes{j} ' - ' num2str( round(1000 * sum(processes(indxSubTypeSorted, i)))/10) '% of all mutations' ],'FontWeight','bold','FontSize', titleFont );

          % Plot error bars
          if ~verLessThan('matlab', '8.4')
               x =  handles.bars(1).XData + handles.bars(1).XOffset;
          else
               x =get(get(handles.bars(1),'children'), 'xdata');
               x = mean(x([1 3],:));
          end
          
          handles.errors(i) = errorbar(x,  processes(indxSubTypeSorted, i), errorsPlotType, 'k', 'linestyle', 'none', 'linewidth', 0.5);
          hold off;
          axis([ 0.5 (length(indxSubTypeSorted) +0.5) max(min(processes(:, i)-0.002),0) min(max(processes(:,i)+errors(:,i)/2+0.0003),1.03) ]);
          
          set(subplot1,'XTick', 1:length(indxSubTypeSorted) );
          set(subplot1,'XTickLabel',doubleSortedSubType);

  if ( length(uniqueTypes) == 6 )
                rotateticklabel(subplot1, 90, fontLabels);
          [sortedProcessDisplay sortedProcessDisplayIndex] = sort(processes(indxSubTypeSorted, i), 'descend');
          for iDisplay = 1 : totalDisplayElements
              displaylabel(iDisplay) = strcat(num2str(iDisplay), ' - ', doubleSortedSubType(sortedProcessDisplayIndex(iDisplay)), ':', num2str(round(10000*processes(indxSubTypeSorted(sortedProcessDisplayIndex(iDisplay)), i))/100), '%');
          end
   
          annotation(f1,'textbox',...
                texboxesLocations(j,:),...
                'String', displaylabel,...
                'LineStyle','none',...
                'Color',[1 0 0]);
            else
rotateticklabel(subplot1, 90, fontLabels);
     end
      end
      annotation(f1,'textarrow',[0.0585937500000001 0.12890625],...
    [0.535290355247701 0.582683720176611],'TextEdgeColor','none',...
    'TextRotation',90,...
    'FontWeight','bold',...
    'FontSize',20,...
    'String', processNames{i},...
    'HeadStyle','none',...
    'LineStyle','none');

annotation(f1,'textbox',...
    [0.003 0.96 0.11 0.043],...
    'String',{[num2str(size(input.sampleNames, 1)) ' ' cancerType ' samples'], [ThousandSep(sum(sum(input.originalGenomes))) ' mutations'], ['Stability ' num2str(round(processStabAvg(i)*100)/100)]},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor','none',...
    'Color',[1 0 0]);

  end
  
  
end