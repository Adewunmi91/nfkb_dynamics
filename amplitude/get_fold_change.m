function [fold_change,max_amp,min_amp] = get_fold_change(X)
% Calculates max fold change for time_series matrix
%INPUT
%   X:          n x t matix
%OUTPUT
%   fold_change:  n x 1 vector
%   max_amp:      n x 1 vector
%   min_amp:      n x 1 vector

% Initialize fold_change 
fold_change = zeros(size(X,1),1);
max_amp=zeros(size(X,1),1);
min_amp=zeros(size(X,1),1);
for row = 1:size(X,1)
    [max_val, max_frame]=max(X(row,:));
    gt_zero= X(row,1:max_frame) >0;
     min_val = min(X(row,gt_zero)); 
     if isempty(min_val)
         continue
     end
    fold_change(row,:) =max_val/min_val;
    max_amp(row,:)=max_val;
    min_amp(row,:)=min_val;
end
end