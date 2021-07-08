function [K, P] = linard2VardistPsi1Compute(linard2kern, vardist, Z)

% LINARD2VARDISTPSI1COMPUTE description.
%
%	Description:
%	[K, P] = linard2VardistPsi1Compute(linard2kern, vardist, Z)
%% 	linard2VardistPsi1Compute.m SVN version 583
% 	last update 2011-06-14T09:16:40.000000Z


K = kernCompute(linard2kern,vardist.means,Z);

P = [];