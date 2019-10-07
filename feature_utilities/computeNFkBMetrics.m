function metrics = computeNFkBMetrics(Trajectories, varargin)
% Compute NFkB metrics: duration, envelope, 
% Adapted from nfkbmetrics
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NFKBMETRICS uses the see_nfkb_native function to filter and preprocess NFkB trajectories,
% then calculates related metrics regarding activation. Metric types include:
% 
% 1) time series (base NFkB dynamics, resampled to 12 frames/hr
% 2) integrated activity
% 3) differentiated activity
% 4) calculated metrics: measuring aspects of oscillation, duration, timing ,and amplitude
%
%
% INPUT PARAMETERS (optional; specify with name-value pairs)
% 'Baseline': 
% 'Interval':   spacing between each time point in hours    
% 'Duration':   Experiment duration  in hours    
% 'Trim':  Time to trime  
% OUTPUT: 
% metrics   structure with output fields
% aux       Extra data (e.g. fourier information (FFT, power, frequencies), thresholds used in envelope/duration)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%% INPUT PARSING
% Create input parser object, add required params from function input
p = inputParser;
addRequired(p, 'Trajectories',@isnumeric);
addParameter(p,'Baseline', 1.0, @isnumeric);
addParameter(p,'Interval', 1/12, @isnumeric);
addParameter(p,'Duration', 8, @isnumeric);
addParameter(p,'Trim', 6, @isnumeric);   
parse(p,Trajectories, varargin{:})

%% PARAMETETERS for finding off times - chosen using 'scan_off_params.m'
Baseline = p.Results.Baseline; % Minimum activity required for cell to register as 'on'
window_sz = 14; % ~1 hr windows (on either side of a given timepoint)
thresh = 0.9; % Pct of inactivity allowed in a given window
cutoff_time = 4; % time to look for cell activity before declaring it "off" (hrs)
off_pad = 12; % Signal time added to trajectory in  FFT calculation (keeps transients from being recorded as osc.)

%% INITIALIZATION. Load and process data. Interpolate time series, calculate deriv/integral approximations

%%
% 1) basic time series. Interpolate over "normal" interval (12 frames per hr) if required
[M,N] = size(Trajectories);
interval = p.Results.Interval; 
endTime = min(p.Results.Trim, p.Results.Duration);
t = 0:interval:endTime;
metrics.time_series = zeros([M, numel(t)]) ;
x =linspace(0, p.Results.Duration, N);
sampFq = round(1/interval); 
if length(t)~=N
    for i = 1:M
        metrics.time_series(i,:) = interp1(x,Trajectories(i,:),t);
    end
end


% 2) integrated activity
[M,N] = size(metrics.time_series);
metrics.integrals = nan([M,N]);
for i = 1:M
    metrics.integrals(i,:) = cumtrapz(t,metrics.time_series(i,:));
end

% 3) differentiated activity - use central finite difference
smoothed = medfilt1(metrics.time_series,3,[],2);
metrics.derivatives = (smoothed(:,3:end) - smoothed(:,1:end-2))/(1/6);

%% MISC METRICS

% 4) Integrals within one-hour windows (0-1, 1-2, 2-3) and three hour windows (0-3, 1-4, etc) of activity
max_hr = floor(max(t));
metrics.intwin1 = nan(size(metrics.time_series,1),max_hr);
metrics.intwin3 = nan(size(metrics.time_series,1),max_hr-2);
for i = 1:(max_hr)
    win = t>=(i-1) & t<(i);
    metrics.intwin1(:,i) = trapz(t(win),metrics.time_series(:,win),2);
    if i<= (max_hr-2)
        win = t>=(i-1) & t<(i+2);
        metrics.intwin3(:,i) = trapz(t(win),metrics.time_series(:,win),2);
    end
end

% 5) amplitude/peak/on-vs-off metrics
% MAX/MIN metrics
metrics.max_amplitude = nanmax(metrics.time_series,[],2);
metrics.max_integral = nanmax(metrics.integrals,[],2);
metrics.max_derivative = nanmax(metrics.derivatives,[],2);
metrics.min_derivative = nanmin(metrics.derivatives,[],2);

% ACTIVITY Compute an off-time for all cells
metrics.off_times = zeros(size(smoothed,1),1);
inactive = [repmat(nanmin(smoothed(:,1:7),[],2),1,window_sz*2+1),smoothed(:,:),...
    repmat(nanmedian(smoothed(:,(end-window_sz:end)),2),1,window_sz*2)];
inactive = smoothrows(inactive<(Baseline),(window_sz*2));
frontcrop = round(window_sz*2*(1-thresh))+window_sz+1;
inactive = inactive(:,frontcrop:end);
inactive = inactive(:,1:size(smoothed,2));
inactive(isnan(smoothed)) = nan;

% Find the final time each cell was active
for i = 1:length(metrics.off_times)
    active_times = find(inactive(i,:)<thresh);
    if ~isempty(active_times)
        if active_times(1) < (cutoff_time*sampFq) % ignore cells who only turned on after 6+ hrs.
            metrics.off_times(i) = active_times(end);
        end
    end    
end
metrics.off_times = (metrics.off_times-1)/sampFq;
metrics.off_times(metrics.off_times<0) = 0;

%% METRICS OF OSCILLATION
% Calculate fourier distribution (via FFT) & power
Fs = 1/300;
depth = max(metrics.off_times)*sampFq;
NFFT = 2^nextpow2(depth); % Next power of 2 from chosen depth
aux.fft = zeros(size(metrics.time_series,1),NFFT/2+1);
aux.freq = Fs/2*linspace(0,1,NFFT/2+1);
aux.power = zeros(size(aux.fft));


for i = 1:size(metrics.time_series,1)
    if(metrics.off_times(i)>0)
        y = metrics.time_series(i,1:depth);
        off_frame = min([length(y), metrics.off_times(i)*sampFq+1+off_pad]); % (Pad w/ 1 extra hr of content)
        y(off_frame:end) = nan;
        y(isnan(y)) = [];
        y = y-nanmean(y);
        if ~isempty(y)
            Y = fft(y,NFFT)/length(y);
            aux.fft(i,:) = abs(Y(1:NFFT/2+1));
            aux.power(i,:) = abs(Y(1:NFFT/2+1).^2);
        end
    end
end

% Find the point of peak (secondary) power
metrics.peakfreq = nan(size(aux.power,1),1);
for i =1:size(metrics.time_series,1)
    [pks,locs] = globalpeaks(aux.power(i,:),2);
    % Ensure we're not getting a totally spurious peak
    if min(pks) < (0.1*max(pks))
        locs(pks==min(pks)) = [];
    end
    if length(locs)>1
        idx = max(locs(1:2));
        metrics.peakfreq(i) = 3600*aux.freq(idx);
    elseif ~isempty(locs)
         metrics.peakfreq(i) = 3600*aux.freq(max([locs,3]));
    else
        metrics.peakfreq(i) = 3600*aux.freq(1);
    end
end
%%
% Find total oscillatory content of particular cells (using thresholds from 0.35 to 0.7 hrs^(-1))
freq_thresh = aux.freq( (aux.freq >= (0.35/3600)) & (aux.freq <= (0.7/3600)));
metrics.oscfrac = nan(size(aux.power,1),length(freq_thresh));
for j = 1:length(freq_thresh)
    for i =1:size(metrics.time_series,1)
        metrics.oscfrac(i,j) = nansum(aux.power(i,aux.freq >= freq_thresh(j))) /nansum(aux.power(i,:));
        if isnan(metrics.oscfrac(i,j))
            metrics.oscfrac(i,j) = 0;
        end
    end
end

%% METRICS OF AMPLITUDE AND TIMING
% 1st + 2nd peak time/amplitude
metrics.pk1_time = nan(size(metrics.time_series,1),1);
metrics.pk1_amp =  nan(size(metrics.time_series,1),1);
metrics.pk2_time = nan(size(metrics.time_series,1),1);
metrics.pk2_amp =  nan(size(metrics.time_series,1),1);

for i = 1:size(metrics.pk1_time,1)    
    [pks, locs] = globalpeaks(metrics.time_series(i,1:min([90,N])),5);
    % Supress any peaks that are within 6 frames of each other.
    [locs, order] = sort(locs,'ascend');
    pks = pks(order);
    while min(diff(locs))<6
        tmp = find(diff(locs)==min(diff(locs)),1,'first');
        tmp = tmp + (pks(tmp)>=pks(tmp+1));
        pks(tmp) = [];
        locs(tmp) = [];  
    end
    pks(locs<4) = [];
    locs(locs<4) = [];
    if ~isempty(locs)
        metrics.pk1_time(i) = locs(1);
        metrics.pk1_amp(i) = pks(1);
    end
    if length(locs)>1
        metrics.pk2_time(i) = locs(2);
        metrics.pk2_amp(i) = pks(2);
    end
end
metrics.pk1_time = (metrics.pk1_time-1)/sampFq;
metrics.pk2_time = (metrics.pk2_time-1)/sampFq;


%% METRICS OF DURATION
% Envelope width: maximum consecutive time above a threshold (envelope must begin within 1st 6 hrs)
smoothed2 = medfilt1(metrics.time_series,5,[],2);
aux.thresholds = linspace(0,Baseline*3,25);
metrics.envelope = zeros(size(metrics.time_series,1),length(aux.thresholds));
for j = 1:length(aux.thresholds)
    thresholded = smoothed2>aux.thresholds(j);
    for i = 1:size(thresholded,1)
        curr = 1;
        idx_start = 1;
        while (curr<size(thresholded,2)) && (idx_start< (6*sampFq))
            idx_start = find(thresholded(i,curr:end)==1,1,'first')+curr-1;
            if ~isempty(idx_start)
                idx_stop = find(thresholded(i,idx_start:end)==0,1,'first')+idx_start-1;
                if isempty(idx_stop)
                    idx_stop = find(~isnan(thresholded(i,:)),1,'last');
                end
                if (idx_stop-idx_start) > metrics.envelope(i,j)
                    metrics.envelope(i,j) = (idx_stop-idx_start);
                end
                curr = idx_stop;
            else
                break
            end
        end
    end
end
metrics.envelope = metrics.envelope/sampFq;


% Number of frames above a given threshold
metrics.duration = zeros(size(metrics.time_series,1),length(aux.thresholds));
for i = 1:length(aux.thresholds)
    metrics.duration(:,i) = nansum(smoothed>aux.thresholds(i),2)/sampFq;
end

end