function model = modelLoadResult(type, dataSet, number, dataLoaderStr)

% MODELLOADRESULT Load a previously saved result.
%
%	Description:
%
%	MODEL = MODELLOADRESULT(TYPE, DATASET, NUMBER) loads a previously
%	saved model result.
%	 Returns:
%	  MODEL - the saved model.
%	 Arguments:
%	  TYPE - the type of model to load.
%	  DATASET - the name of the data set to load.
%	  NUMBER - the number of the model run to load.
%	
%
%	See also
%	MODELWRITERESULT


%	Copyright (c) 2009 Neil D. Lawrence
% 	modelLoadResult.m SVN version 1521
% 	last update 2011-07-21T10:11:54.000000Z
  if nargin < 4
    dataLoaderStr = 'lvmLoadData';
  end
  dataLoaderHandle = str2func(dataLoaderStr);
  origDataSet = dataSet;
  dataSet(1) = upper(dataSet(1));
  origType = type;

  type(1) = upper(type(1));
  fileName = ['dem' dataSet type num2str(number)];

  fhandle = [origType 'LoadResult'];
  if exist(fhandle)==2
    % There is load result code, use it.
    fhandle = str2func(fhandle);
    if nargin > 3
      model = fhandle(origDataSet, number, dataLoaderStr);
    else
      model = fhandle(origDataSet, number);
    end
  else
    fhandle = [origType 'Reconstruct'];
    if exist(fhandle)==2
      % There is reconstruct code, use it to reconstruct.
      [Y, lbls] = dataLoaderHandler(origDataSet);
      load(fileName);
      varName = [origType 'Info'];
      eval(['model = ' fhandle '(' varName ', Y);'])
    else
      % No code is provided, just load the file.
      load(fileName);
    end
  end
end