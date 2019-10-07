function duration = computeDuration_v1(Trajectories, varargin)
% Compute duration
% -----------------------------------------------------------------------
% INPUT:
%   REQUIRED:
%       Trajectories:
%   OPTIONAL:
%       Deriv:  2
%       Verbose:  false
%       'Threshold': 0.5
p = inputParser;
addRequired(p, 'Trajectories', @isnumeric);
addOptional(p,'Deriv', 0, @isnumeric);
addOptional(p, 'Verbose',false ,@islogical);
addParameter(p, 'Threshold',  0.1, @isnumeric);
addParameter(p, 'UnitConversionFun', @frame2hr);

parse(p, Trajectories, varargin{:});
data = p.Results.Trajectories;
if isvector(data)
    data = toRow(data);
end
thresh = p.Results.Threshold;
if numel(data) <2
    info = read_id(data);
    data= get_trajectories(data, 'TrimFrame',max(info.TimeRange{:}),...
        'SortBy', {'off_times'} , 'Order', 'ascend', 'Cache', true);
end

clr = setcolors;
%2019-03-29: removing wavelet due to performance issues
% smoothFuns = {'wavelet','lowess', 'median'};

convertFun = p.Results.UnitConversionFun;
neg = @(j) j*-1;
deriv = p.Results.Deriv;
normFun = @(x) normalize(x, 'range');
duration = zeros(size(data,1),1);

for i =1: size(data,1)
    %%
    endpt =find(isfinite(data(i,:)),1,'last');
    X = fillmissing(data(i,1:endpt), 'linear')';
    t = transpose(1:endpt);
    f=fit(t,X , 'smoothingspline','SmoothingParam',0.2);
    x0 = feval(f, t);
    [x1,x2] = differentiate(f,t);
    grid on;
    
    %%  Find strongest decline rate
            funs = {normFun};
            z0 = chain(x0, funs{:});
    switch deriv
        case 0
             x = x0; pkArgs = {'MinPeakProminence',thresh, 'MinPeakWidth',6};
        case 1
            funs = {neg, normFun};
            x = x1;
            pkArgs = {'MinPeakProminence',thresh, 'MinPeakWidth',3};
        case 2
            x = x2;
            funs = {neg, @zscore,normFun};
            pkArgs = {'MinPeakProminence',thresh, 'MinPeakWidth',3};
    end
    
    Z = chain(x,funs{:});
    [~,zPkLoc] = findpeaks(Z,pkArgs{:});
    if isempty(zPkLoc) % if too stringent, find peaks
        [~,zPkLoc]= findpeaks(Z);
    end
    if isempty(zPkLoc)
        duration(i) = 0;
    else
        %% Find where trajectories taper off, i.e. derivative => 0
        pkIx = max(zPkLoc);
        % find where the signal crosses the base of the last peak prominence
        [~, maxDur]=getPeakBounds(z0,pkIx);
        
        if isempty(convertFun)
            duration(i) = maxDur;
        else
            duration(i) = convertFun(maxDur);
        end
    end
    %% Plots
    if p.Results.Verbose
        fig=figure('Position', get_scr_sz());
        subplot(2,2,1)
        plotProps = {'XLim', [0, numel(X)], 'YGrid', 'on', 'XGrid', 'on'};
        plot(X)
        title('Original')
        hold on
        
        ax1 = gca;
        plot(repelem(maxDur,2), ax1.YLim, ':', 'LineWidth', 3);
        
        
        
        set(gca, plotProps{:});
        subplot(2,2,2)
        title('Fit')
        
        plot(t,X, 'LineStyle', 'none', 'Marker','o','LineWidth', 3, 'MarkerSize', 3,...
            'MarkerEdgeColor', clr.orange, 'MarkerFaceColor', clr.orange ); hold on;
        plot( t,x0);     hold on; plot(t,Z); hold on;
        
        
        legend({'orig','fit','transform'});
        
        hold on;plot(repelem(maxDur,2), ax1.YLim, ':', 'LineWidth', 3);
        set(gca, plotProps{:});
        subplot(2,2,3)
        
        
        findpeaks(Z,'Annotate','extents', pkArgs{:})
        title('Peaks')
        hold on
        
        ax2= gca;
        plot(repelem(maxDur,2), ax1.YLim, ':', 'LineWidth', 3);
        
        
        
        set(gca, plotProps{:});
        keyboard
        delete(fig);
        
    end
end
end


function[leftBound, rightBound]= getPeakBounds(x,peakLoc)
%Returns the left boundary & right boundary of a peak
params = {'MinSeparation', 6, 'FlatSelection', 'first'};

leftBound=peakLoc - find(islocalmin(flip(x(1:peakLoc)), params{:}),1,'first')+1;
rightBound = find(islocalmin(x(peakLoc:end), params{:}),1, 'first') + peakLoc-1;
leftBound = max(leftBound, 1);
if isempty(rightBound)
    rightBound =  numel(x);
end
% base = x(peakLoc) - p;
% leftBound=peakLoc - find(flip(x(1:peakLoc))<=base,1,'first')+1;
% rightBound=find(x(peakLoc:end)<=base,1,'first') + peakLoc-1;
end






