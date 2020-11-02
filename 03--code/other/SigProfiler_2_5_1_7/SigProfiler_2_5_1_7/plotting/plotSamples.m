function plotSamples(samples, types, subtypes, sampleNames)

  %samples = samples ./ sum(samples);
  screen_size = get(0, 'ScreenSize');
  totalProcesses = size(samples, 2);
  totalMutationTypes = size(samples, 1);
  topElements = 16;
  totalDisplayElements = 10;
  titleFont = 20;
  fontLabels = 4;
  [sortedType sortedIndex] = sort(types);
  uniqueTypes = unique(sortedType);
  sortedSubType = subtypes(sortedIndex);
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
  
  if ( size(types, 1) < 200 )
      fontLabels = 12;
      totalDisplayElements = 3;
      titleFont = 15;
  end
  
  
  if ( size(types, 1) == 96 )
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
  
  
  % Plot mutational samples and their features
  for i = 1 : totalProcesses
      tempFileName = ['temp/process_ ' num2str(i) '.txt'];
      fileID = fopen(tempFileName,'w');
      fprintf(fileID,'types\tSubtype\tContribution\n');
      
      f1 = figure('InvertHardcopy','off','Color',[1 1 1]);
      set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
      [sortedProcess sortedProcessIndex] = sort(samples(:, i)); % sortedProcess is unused
      minTopElementValue = samples(sortedProcessIndex(totalMutationTypes-topElements), i);
      
      for j = 1 : length(uniqueTypes)
          subplot1 = subplot(subplotLeft, subplotRight, j,'LineWidth',1,'FontWeight','bold', 'FontSize',20);
          indxSubType = find ( strcmp(sortedType, uniqueTypes{j}) );
          [doubleSortedSubType doubleSortedSubTypeIndex] = sort(sortedSubType(indxSubType));
          indxSubTypeSorted = indxSubType(doubleSortedSubTypeIndex);
          
          errorsPlotType = errors(indxSubTypeSorted, i);
          
          handles.bars = bar( 1:length(indxSubTypeSorted), samples(indxSubTypeSorted, i), 'FaceColor', colors(j, :), 'EdgeColor', edgeColor, 'BarWidth', 0.4);
   
          hold on;
          plot(0.5:0.5:(length(indxSubTypeSorted) +0.5), minTopElementValue * ones(1,length(0.5:0.5:(length(indxSubTypeSorted) +0.5))), 'LineStyle','-.','Color',[1 0 0]);
     
          title( [ uniqueTypes{j} ' - ' ThousandSep( sum(samples(indxSubTypeSorted, i))) ' mutations' ],'FontWeight','bold','FontSize', titleFont );

          % Plot error bars
          x =get(get(handles.bars(1),'children'), 'xdata');
          x = mean(x([1 3],:));
          handles.errors(i) = errorbar(x,  samples(indxSubTypeSorted, i), errorsPlotType, 'k', 'linestyle', 'none', 'linewidth', 0.5);
          hold off;
          
         axis([ 0.5 (length(indxSubTypeSorted) +0.5) 0 max(samples(:,i)+errors(:,i)/2+5) ]);
          
          set(subplot1,'XTick', 1:length(indxSubTypeSorted) );
          set(subplot1,'XTickLabel',doubleSortedSubType);
          for iPrint = 1 : length(indxSubTypeSorted)
              fprintf(fileID,'%s\t%s\t%f\n', uniqueTypes{j}, doubleSortedSubType{iPrint}, samples(indxSubTypeSorted(iPrint), i));
          end

  if ( length(uniqueTypes) == 6 )
                rotateticklabel(subplot1, 90, fontLabels);
          [sortedProcessDisplay sortedProcessDisplayIndex] = sort(samples(indxSubTypeSorted, i), 'descend');
          for iDisplay = 1 : totalDisplayElements
              displaylabel(iDisplay) = strcat(num2str(iDisplay), ' - ', doubleSortedSubType(sortedProcessDisplayIndex(iDisplay)), ':', num2str(samples(indxSubTypeSorted(sortedProcessDisplayIndex(iDisplay)))));
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

annotation(f1,'textbox',...
    [0.003 0.96 0.11 0.043],...
    'String',{sampleNames{i}, [ThousandSep(sum(samples(:,i))) ' mutations']},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor','none',...
    'Interpreter', 'none', ...
    'Color',[1 0 0]);
fclose(fileID);

  end
  
 
  
  
end