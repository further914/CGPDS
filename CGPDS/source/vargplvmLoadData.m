function [Y, lbls, Ytest, lblstest] = vargplvmLoadData(dataset, local, seedVal)

% VARGPLVMLOADDATA Load a latent variable model dataset from a local folder
%
%	Description:
%	or from the global repository. This function tries to load the file from
%	the following directories (in that order): the local directory where
%	the small datasets are kept, the local directory where the large files
%	are kept, the global directory as recognised in lvmLoadData, the current
%	directory. If the file is found while searching in any of the above
%	directories, the function returns so that no more searching is done.
%	Searching in local folders succeeds only if such a local folder is added in the
%	path and loadLocalData.m is also added in the path (e.g. inside one of
%	the local data folders). The loadLocalData.m function has exactly the
%	same form as lvmLoadData.m.
%
%	[Y, LBLS, YTEST, LBLSTEST] = VARGPLVMLOADDATA(DATASET, SEEDVAL,
%	LOCAL) loads a data set for a latent variable modelling problem.
%	 Returns:
%	  Y - the training data loaded in.
%	  LBLS - a set of labels for the data (if there are no labels it is
%	   empty).
%	  YTEST - the test data loaded in. If no test set is available it is
%	   empty.
%	  LBLSTEST - a set of labels for the test data (if there are no
%	   labels it is empty).
%	 Arguments:
%	  DATASET - the name of the data set to be loaded.
%	  SEEDVAL - set a value for the random seeds.
%	  LOCAL - set to true/false for the function to search for dataset
%	   in the local folders. Default is true.
%	Copyright: Andreas C. Damianou, 2011
%	
%	
%
%	See also
%	LVMLOADDATA
% 	vargplvmLoadData.m SVN version 1570
% 	last update 2011-08-30T14:57:48.000000Z

if nargin > 1
    searchLocally = local;
else
    searchLocally = 1;
end

if nargin > 2
    randn('seed', seedVal)
    rand('seed', seedVal)
end



% get directory
% if strcmp(fileCategory,'small')
%     baseDir = localDatasetsDirectorySmall;
% elseif strcmp(fileCategory,'large')
%     baseDir = localDatasetsDirectoryLarge;
% end

dirSep = filesep;

% Try the local directories
if searchLocally
    % First, try the local directory (if exists) where the small datasets
    % are kept
    try
        [Y, lbls, Ytest, lblstest]= loadLocalData(dataset, localDatasetsDirectorySmall,dirSep);
        fprintf('# The requested dataset was found in the local directory (for the small files).\n');
        return
    catch
        % do nothing
    end

    % If the above fails, try the local directory (if exists) where the
    % large datasets are kept
    try
        [Y, lbls, Ytest, lblstest]=loadLocalData(dataset, localDatasetsDirectoryLarge,dirSep);
        fprintf('# The requested dataset was found in the local directory (for the large files).\n');
        return
    catch
        % do nothing
    end
end

% If we reach here nothing of the above contains the dataset. Try the
% global directory.
try
    switch nargout
        case 1
           %% modified by ZhaoJing  lvmloadData->load
            load(dataset);
            % 2013-05-05
            [Y]=data;
        case 2
            [Y, lbls] = lvmLoadData(dataset);
        case 3
            [Y, lbls,Ytest] = lvmLoadData(dataset);
        case 4
            [Y, lbls, Ytest, lblstest] = lvmLoadData(dataset);
    end
    fprintf('# The requested dataset was loaded from the global DATASETS directory.\n');
    return
catch
    % do nothing
end

% If everything else fails, try to just load the dataset from the current
% dir.
load(dataset);
