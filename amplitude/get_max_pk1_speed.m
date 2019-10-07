function maxSpeed = get_max_pk1_speed(Metrics, varargin)
% -------------------------------------------------------------------------
% maxSpeed = get_max_pk1_speed(metrics)
% -------------------------------------------------------------------------
% INPUT
%       REQUIRED:
%             Metrics: struct array containg 'derivatives' and 'pk1_time'
%             fields
% OUTPUT
%     maxSpeed: maximum derivative preceding the first peak
p = inputParser; 
validStruct = @(x) isstruct(x) && (isfield(x, 'derivatives') && isfield(x, 'pk1_time'));
addRequired(p, 'Metrics', validStruct);
addParameter(p, 'Debug', false, @islogical);
parse (p, Metrics, varargin{:});
Metrics = p.Results.Metrics;
maxSpeed = zeros(size(Metrics.pk1_time));
plotProps = {'XLim', [0, size(Metrics.time_series, 2)], 'YGrid', 'on', 'XGrid', 'on'};
time2frame = @(x) (x *12) +1;
clr = setcolors;     
% smoothFun = {'none'}; 
smoothFun = {'wavelet', 'lowess'};
smoothed = smooth_trajectory( Metrics.derivatives,'Method',smoothFun );
    
for i =1 :numel(maxSpeed)
    pk1_frame= time2frame( Metrics.pk1_time(i,:));   
    [maxSpeed(i),loc] =max(smoothed(i,1:pk1_frame));
    if p.Results.Debug        
       figure('Position', get_scr_sz())
  
       X = chain(Metrics.time_series(i,:), @(x) smooth_trajectory(x, 'Method', smoothFun));
       
        
       plot(X, 'LineWidth', 3);
       hold on
       plot(loc,X(loc+1), 'o', 'LineWidth', 3, 'MarkerSize', 12, ...
           'MarkerEdgeColor', clr.orange, 'MarkerFaceColor', clr.orange);
          hold on
       plot(pk1_frame,X( pk1_frame), 'o', 'LineWidth', 3, 'MarkerSize', 12, ...
           'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
       set(gca, plotProps{:})
       keyboard
    end
    
end
end