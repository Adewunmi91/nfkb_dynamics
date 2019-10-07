function T = getFoldChanges(X, varargin)
% Calculates max fold change for time_series matrix
%INPUT
%   X:          n x t matix
%OUTPUT
%   Y= table; 
%   fold_change:  n x 1 vector
%   max_amp:      n x 1 vector
%   min_amp:      n x 1 vector

% Initialize fold_change 
p = inputParser; 
addRequired(p,'X', @isnumeric);
addOptional(p,'XMax', [],@isnumeric);
addOptional(p, 'Frame',[],@isnumeric);
parse(p, X, varargin{:});
% X = p.Results.X + min(

[gMin, gMinLoc] = min(X,[],2);


allFC = fillmissing(standardizeMissing(XMax.\X,[nan,Inf]) ,'constant',0);
[maxFrac,maxLoc] = nanmax(allFC,[],2);
assert(all(maxFrac==1),'Sanity check failed')
figure; histogram(rmoutliers(maxFC))
T = table('Size', [], 'VariableTypes','numeric', 'VariableName', ["MaxFoldChange", ]); 

end