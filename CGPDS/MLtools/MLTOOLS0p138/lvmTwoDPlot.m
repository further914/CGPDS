function [returnVal, txtReturnVal] = lvmTwoDPlot(X, lbl, symbol)

% LVMTWODPLOT Helper function for plotting the labels in 2-D.
%
%	Description:
%
%	LVMTWODPLOT(X, LBL, SYMBOL) helper function for plotting an
%	embedding in 2-D with symbols.
%	 Arguments:
%	  X - the data to plot.
%	  LBL - the labels of the data point.
%	  SYMBOL - the symbols to use for the different labels.
%	
%
%	See also
%	LVMSCATTERPLOT, LVMVISUALISE


%	Copyright (c) 2004, 2005, 2006, 2008 Neil D. Lawrence
% 	lvmTwoDPlot.m CVS version 1.5
% 	lvmTwoDPlot.m SVN version 589
% 	last update 2009-11-08T13:06:39.000000Z

  if nargin < 2
    lbl = [];
  end
  if(strcmp(lbl, 'connect'))
    connect = true;
    lbl = [];
  else
    connect = false;
  end
  
  if iscell(lbl)
    lblStr = lbl{2};
    lbl = lbl{1};
    labelsString = true;
  else
    labelsString = false;
  end
  
  if nargin < 3
    if isempty(lbl)
      symbol = getSymbols(1);
    else
      symbol = getSymbols(size(lbl,2));
    end
  end
  axisHand = gca;
  returnVal = [];
  textReturnVal = [];
  nextPlot = get(axisHand, 'nextplot');
  for i = 1:size(X, 1)
    if i == 2
      set(axisHand, 'nextplot', 'add');
    end
    if ~isempty(lbl)
      labelNo = find(lbl(i, :));
    else
      labelNo = 1;
    end
    %try 
    returnVal = [returnVal; plot(X(i, 1), X(i, 2), symbol{labelNo})];
    if labelsString
      textReturnVal = [textReturnVal; text(X(i, 1), X(i, 2), ['   ' lblStr{i}])];
    end
    if connect
      if i>1
        line([X(i-1, 1) X(i, 1)], [X(i-1, 2) X(i, 2)]);
      end
    end
    %catch
    % if strcmp(lasterr, 'Index exceeds matrix dimensions.')
    %  error(['Only ' num2str(length(symbol)) ' labels supported (it''s easy to add more!)'])
    % end
  %end
  end
  set(axisHand, 'nextplot', nextPlot);
  set(returnVal, 'markersize', 10);
  set(returnVal, 'linewidth', 2);
end
