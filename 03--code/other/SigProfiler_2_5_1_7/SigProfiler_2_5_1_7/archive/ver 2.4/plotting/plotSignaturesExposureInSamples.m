%% This function plot genomes based the contributions of the mutational 
%% processes operative in them
function plotSignaturesExposureInSamples(exposure, input, processNames)
   
  % Defining constants
  if ( isfield(input, 'cancerType') == 0 )
    cancerType = 'Unknown';
  else
    cancerType = input.cancerType;  
  end
  
  totalProcesses = size(exposure, 1);
  totalGenomes = size(exposure, 2);
  maxPerPlot = 50;
  
  colors(1, :) = [ 0.00 0.00 0.50 ]; % navy blue
  colors(2, :) = [ 0.00 0.50 0.00 ]; % navy blue
  colors(3, :) = [ 0.50 0.00 0.00 ]; % navy blue
  colors(4, :) = [ 1.00    0.5000    0.0000 ]; %light coral
  colors(5, :) = [  0.2824    0.8196    0.8000 ];% medium turquoise
  colors(6, :) = [ 0  0 0 ];
  colors(7, :) = [ 1.0000    0.0784    0.5765 ];
  colors(8, :) = [ 0.50 0.25 0.60 ];
  
  [sortedTotalExp, sortedTotalExpInd] = sort(sum(exposure,1), 'descend'); % sortedProcess is unused);
  %sortedTotalExpInd = 1 : size(exposure,2);
  genomeNames = input.sampleNames(sortedTotalExpInd);
  genomes = exposure(:, sortedTotalExpInd);
  screen_size = get(0, 'ScreenSize');
  genomes = genomes + 10e-5;

  generateNames = 0;
  generateGenomeNames = 0;
  
  if ( exist('processNames', 'var') == 0 )
      generateNames = 1;
  end
  
  if ( exist('genomeNames', 'var') == 0 )
      generateGenomeNames = 1;
  end

  totalContribution = 0;
  for i = 1 : totalProcesses
      contribution = round(100*sum( genomes(i,:) )/sum(sum(genomes)));
     if ( i == totalProcesses)
         contribution = 100 - totalContribution;
     end
     
     if ( generateNames == 1 )
         processNames{i} = ['Signature ' num2str(i) ' (' num2str(contribution) '%)'];
     else
         processNames{i} = [ processNames{i} ' (' num2str(contribution) '%)'];
     end
     totalContribution = totalContribution + contribution;
     
  end
  
  for i = 1 : totalGenomes
      if ( generateGenomeNames == 1 )
          genomeNames{i} = ['S' num2str(i)];
      end
  end

  if ( exist('genomeNames', 'var') == 0 )
      for i = 1 : totalGenomes
          genomeNames{i} = [ 'Sample ' num2str(i)];
      end
  end
  
  % Calculating normalized genomes
  genomesNormalized = zeros( size(genomes) );
  for i = 1 : totalGenomes
     genomesNormalized(:, i) = genomes(:, i) ./ sum( genomes(:,i) );
  end
  genomesNormalized = genomesNormalized + 10e-5;
  
  for manyPlots = 1 : ceil(totalGenomes/maxPerPlot)
  minGenomeSet  = (manyPlots-1)*maxPerPlot + 1;
  maxGenomeSet  = min(maxPerPlot*manyPlots,totalGenomes);
  if ( manyPlots ~= ceil(totalGenomes/maxPerPlot) )
    maxPerPlotSet = maxPerPlot;
  else
    maxPerPlotSet = mod(totalGenomes, maxPerPlot);
    if ( maxPerPlotSet == 0 )
        maxPerPlotSet = maxPerPlot;
    end
  end
  
  genomesSet1 = genomes(:, minGenomeSet:maxGenomeSet );
  
  f1 = figure('InvertHardcopy','off','Color',[1 1 1]);
  set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
  
  subplot1 = subplot(1, 2, 1,'LineWidth',1,'FontWeight','bold', 'FontSize',20);
  bar1 = bar(genomesSet1', 'stacked', 'BarWidth', 0.4);
  for i = 1 : totalProcesses
    set(bar1(i), 'DisplayName', processNames{i});
    set(bar1(i), 'facecolor', colors(i,:), 'edgecolor', colors(i,:))
  end
  
  axis tight;
%   colormap(jet);
  title('Mutations in genomes');
  legend1 = legend(subplot1,'show');
  set(legend1,...
    'Position',[0.00335498252500543 0.52930013371481 0.123207517474995 0.181213841578805],...
    'FontWeight','bold',...
    'FontSize',15);
  set(subplot1,'XTick', 1:maxPerPlotSet );
  set(subplot1,'XTickLabel',genomeNames(minGenomeSet:maxGenomeSet));
  rotateticklabel(subplot1, 90, 11);
  
  genomesSet1Norm = genomesNormalized(:, minGenomeSet:maxGenomeSet);
  subplot2 = subplot(1, 2, 2,'LineWidth',1,'FontWeight','bold', 'FontSize',20);
  bar2 = bar(genomesSet1Norm', 'stacked', 'BarWidth', 0.4);
   for i = 1 : totalProcesses
    set(bar2(i), 'facecolor', colors(i,:), 'edgecolor', colors(i,:))
  end
  
  axis tight;
  colormap(jet);
  title('Normalized mutations in genomes');
    set(subplot2,'XTick', 1:maxPerPlotSet );
  set(subplot2,'XTickLabel',genomeNames(minGenomeSet:maxGenomeSet));
  rotateticklabel(subplot2, 90, 11);
  
  annotation(f1,'textbox',...
    [0.003 0.96 0.11 0.043],...
    'String',{[num2str(size(input.sampleNames, 1)) ' ' cancerType{1} ' samples'], [num2str(minGenomeSet) ' to ' num2str(maxGenomeSet)], [ThousandSep(sum(sum(input.originalGenomes))) ' mutations']},...
    'HorizontalAlignment','center',...
    'FontWeight','bold',...
    'FontSize',18,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'EdgeColor','none',...
    'Color',[1 0 0]);
  end
end