function duration = computeDuration(Trajectories, varargin)
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
addParameter(p, 'Threshold',  0.5, @isnumeric);
addParameter(p, 'UnitConversionFun', @frame2hr);

parse(p, Trajectories, varargin{:});

data = p.Results.Trajectories;
if isvector(data)
    data = toRow(data);
end
thresh = p.Results.Threshold;
feats= [];
if numel(data) <2
    info = read_id(data);
    [data, feats ]= get_trajectories(data, 'TrimFrame',max(info.TimeRange{:}), 'SortBy', {'off_times'} , 'Order', 'ascend');
end

clr = setcolors;
%2019-03-29: removing wavelet due to performance issues
% smoothFuns = {'wavelet','lowess', 'median'};
smoothFuns = {'lowess', 'median'};
smooth=@(x) smooth_trajectory(x, 'Method',smoothFuns);
convertFun = p.Results.UnitConversionFun;
neg = @(j) j*-1;
deriv = p.Results.Deriv;
% normFun = @(x) normalize(x, 'zscore','std');
normFun = @(x) normalize(x, 'range');
duration = zeros(size(data,1),1);
for i =1: size(data,1)
    %%
    %3/29/2019: fillmissing
    %     X = data(i,:);
    X = fillmissing(data(i,:),'linear');
     a = [1:numel(X)]'; b =X';  f=fit(a,b , 'linearinterp'); figure;plot(f, a,b)
    endpt =find(isfinite(X),1,'last');
    X = smooth(X(1:endpt));
    
    %%  Find strongest decline rate
    switch deriv
        case 0
            funs = {@(x) medfilt1(x,6), normFun };
            pkArgs = {'MinPeakProminence',thresh/2, 'MinPeakWidth',6,'MinPeakHeight',thresh};
%             valleyArgs = {'MinPeakProminence', 0.1};
                valleyArgs = {'MinPeakProminence',thresh/4};
        case 1
            %             funs = {@diff,smooth, neg, @zscore};
            funs = {@diff,smooth, neg, @zscore, normFun};
            pkArgs = {'MinPeakProminence',thresh/2, 'MinPeakWidth',3,'MinPeakHeight',thresh};
            valleyArgs = {'MinPeakWidth',1};
        case 2
            funs = {@diff,smooth,@diff ,  @zscore, normFun};
            %             pkArgs = {'MinPeakProminence',0.15, 'MinPeakWidth',3,'MinPeakHeight', 0.15};
            pkArgs = {'MinPeakProminence',thresh/2, 'MinPeakWidth',1,'MinPeakHeight',thresh};
            valleyArgs = {'MinPeakProminence',0.1, 'MinPeakWidth',2,'MinPeakHeight', 0};
        case 3
            funs = {@diff, @diff,@diff ,@zscore};
            pkArgs = {'MinPeakProminence',thresh/2};
            valleyArgs = {'MinProminence', 0.4};
    end
    
    Z = chain(X,funs{:});
    warning('off');
    [~,zPkLoc, ~, pr] = findpeaks(Z,pkArgs{:});
    
    if isempty(zPkLoc) % if too stringent, find peaks
        [~,zPkLoc, ~, pr] = findpeaks(Z);
    end
    
    if isempty(zPkLoc)
        duration(i) = 0;
    else
        %% Find where trajectories taper off, i.e. derivative => 0
        [pkIx, maxIx] = max(zPkLoc);
        
        % find where the signal crosses the base of the last peak prominence
        
        tmp =  X(:,pkIx+1:end);
        postLastPeak = neg(tmp) + range(tmp);
        % find valley after the peak
        
        valleyLoc =[];
        if deriv > 0
            [localMinima,valleyLocs, ~,valleyProm] = findpeaks(postLastPeak,...
                valleyArgs{:}) ;
            if any(localMinima)
                [~, ix ] = max(valleyProm);% use for Z
                valleyLoc = valleyLocs(ix);
            end
        else
            [~, rB]=getPeakBounds(Z,pkIx,pr(maxIx));
            valleyLoc = rB - pkIx;
        end
        
        if isempty(valleyLoc)
            [~,valleyLoc] =max(postLastPeak);
        end
        dur = zPkLoc+valleyLoc;
        maxDur = max(dur);
        if maxDur > endpt
            keyboard;
        end
        warning ('on');
        if isempty(convertFun)
            duration(i) = maxDur;
        else
            duration(i) = convertFun(maxDur);
        end
    end
    %% Plots
    if p.Results.Verbose
        disp(dur)
        fig=figure('Position', get_scr_sz());
        subplot(2,2,1)
        plotProps = {'XLim', [0, numel(X)], 'YGrid', 'on', 'XGrid', 'on'};
        plot(data(i,:))
        title('NFkB')
        hold on
        if ~isempty(dur)
            ax1 = gca;
            plot(repelem(max(dur),2), ax1.YLim, ':', 'LineWidth', 3);
        end
        
        if ~isempty(feats)
            ax1 = gca;
            frameIx = round(feats.off_times(i)*12);
            plot(frameIx,X(frameIx), 'o', 'LineWidth', 3, 'MarkerSize', 12, 'MarkerEdgeColor', clr.orange, 'MarkerFaceColor', clr.orange);
        end
        set(gca, plotProps{:});
        
        
        subplot(2,2,2)
        plotProps = {'XLim', [0, numel(X)], 'YGrid', 'on', 'XGrid', 'on'};
        plot(X)
        title('smoothed')
        
        hold on
        if ~isempty(dur)
            ax1 = gca;
            plot(repelem(max(dur),2), ax1.YLim, ':', 'LineWidth', 3);
        end
        
        if ~isempty(feats)
            ax1 = gca;
            frameIx = round(feats.off_times(i)*12);
            plot(frameIx,X(frameIx), 'o', 'LineWidth', 3, 'MarkerSize', 12, 'MarkerEdgeColor', clr.orange, 'MarkerFaceColor', clr.orange);
        end
        set(gca, plotProps{:});
        
        subplot(2,2,3)
        links =chain(X,funs{1:end-1});
        findpeaks(links,'Annotate','extents', pkArgs{:})
        title('Transform')
        hold on
        if ~isempty(dur)
            ax2= gca;
            plot(repelem(max(dur),2), ax2.YLim, ':', 'LineWidth', 3);
        end
        set(gca, plotProps{:});
        subplot(2,2,4)
        findpeaks(postLastPeak,'Annotate','extents', valleyArgs{:})
        keyboard
        delete(fig);
        
    end
end
end


function[leftBound, rightBound]= getPeakBounds(x,peakLoc, p)
%Returns the left boundary & right boundary of a peak

base = x(peakLoc) - p;
leftBound=peakLoc - find(flip(x(1:peakLoc))<=base,1,'first')+1;
rightBound=find(x(peakLoc:end)<=base,1,'first') + peakLoc-1;
end
