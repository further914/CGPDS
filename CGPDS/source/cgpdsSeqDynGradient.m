function g = cgpdsSeqDynGradient(x, model, y,samd,samr)

% cgpdsSeqDynGradient Wrapper function for gradient of CGPDS for testing data.



%	Copyright (c) 2009 % COPYRIGHT Michalis K. Titsias and Neil D. Lawrence
% 	vargplvmSeqDynGradient.m SVN version 1441
% 	last update 2011-06-02T20:53:40.000000Z


vardistx = model.vardistx;
vardistx = vardistExpandParam(vardistx, x);
fhandle = str2func([model.type 'SeqDynLogLikeGradient']);
g = -fhandle(model, vardistx, y,samd,samr);

