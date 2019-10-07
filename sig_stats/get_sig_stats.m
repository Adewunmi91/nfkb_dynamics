function sig_stats = get_sig_stats(Data, varargin)
% Calculates signal statistics for NFkB trajectories
%-----------------------------------------------------------------------
% INPUTS
%   REQUIRED
%       Data:         nxt array 
%   OPTIONAL
%   'FqRange',  frequency range of interest (2-element vector)
%   'Fs'         sample rate
%   'FillMethod' method for filling missing Data 
%                     {'previous',next', 'linear', 'spline', 'pchip'}
%   'stats'      member of get_sig_stat_list. 
%                     default = {'peak2peak','peak2rms','obw'}
%   'SmoothMethod' {'lowess', 'sgolay', 'none', 'wavelet'}
% OUTPUT
%   stats:    struct of wavelet, and spectral statistics
% 12/29/2018: use smoothed data for calculations
%
extraFeats =convertCharsToStrings({'power', 'medfreq', 'meanfreq', 'psd', 'noise_est', 'max_pentropy'})';
is_valid_stat=@(x) all(ismember(x, [string(get_sig_stat_list()); extraFeats])); 
valid_FillMethod  = {'previous', 'next','linear', 'spline', 'pchip'};
isvalidfill =@(x)ismember (x, valid_FillMethod);
p=inputParser;
addRequired(p, 'Data',@isnumeric);
addParameter(p, 'FqRange',[0.33 1] , @isnumeric);
addParameter(p, 'Fs',12, @isnumeric);
addParameter(p, 'FillMethod','linear', isvalidfill);
addParameter(p, 'Stats', ["peak2peak", "peak2rms","obw","bandpower"],is_valid_stat);
addParameter(p, 'SmoothMethod',["sgolay", "lowess"], @istext);

parse(p,Data,varargin{:});

%% Time-Series Statistics
orig_data = Data;
sig_stats=struct; 
stats=convertCharsToStrings(p.Results.Stats);
Data=fillmissing(Data,p.Results.FillMethod, 2, 'EndValues','extrap');
% Data=fillmissing(Data,p.Results.FillMethod, 2, 'EndValues',0);
smoothData = smooth_trajectory(Data,'Method', p.Results.SmoothMethod); 
% smoothData = smooth_trajectory(Data,'Method', {'wavelet', 'lowess'}); 
% 3/29/2018: changed smoothing fun;
Fs=p.Results.Fs; freq_range =p.Results.FqRange;

if any(ismember(stats, ["bandpower","oscfreq", "oscpower"]))
   stats= unique(["psd";stats],'stable'); 
end
for i=1:numel(stats)
    stat = stats(i);
    switch(stat)
        case 'medfreq'
            sig_stats.(stat)=medfreq(smoothData', Fs,freq_range)';
        case 'meanfreq'            
            sig_stats.(stat)=meanfreq(smoothData', Fs, freq_range)'; 
        case 'peak2rms'
            sig_stats.(stat)=peak2rms(smoothData,2);
        case 'rms'            
            sig_stats.(stat)=rms(smoothData')';
        case 'peak2peak'
            sig_stats.(stat) = peak2peak(smoothData,2);
        case 'mean_movmad'
            sig_stats.(stat)=mean( movmad(smoothData,3,2),2);
        case 'mean_movstd'
            sig_stats.(stat)=mean( movstd(smoothData,3,[],2),2);
        case 'mean_movvar'
            sig_stats.(stat) =mean( movvar(smoothData,3,[],2),2);
        case {'psd', 'power'}    
%             [~, sig_stats.(stat)]=pwelch(Data',[],[],[],Fs,'one-sided','power');
%             [pwr, fq] = periodogram(smoothData',[], [], Fs, 'power');      
             n = size(Data,2); 
             [pwr,fq]=pwelch(smoothData',n,10,256,Fs,'one-sided',stat);
            sig_stats.fq=fq';
            %normalize power/psd to 1
            sig_stats.(stat)=transpose(pwr./sum(pwr,1));
        case {'bandpower', 'oscpower'}
            pwr = transpose(sig_stats.psd) ; fq = transpose(sig_stats.fq);%             
             bp= bandpower(pwr,fq,freq_range, 'psd')';
%              totalPower = bandpower(pwr,fq,[min(fq), max(fq)], 'psd')';
%              bpFrac = bp./totalPower;
%              sig_stats.(stat) = bpFrac;
             sig_stats.(stat) =bp;
        case { 'oscfreq'}
            %find peaks within the frequency range
           bandfilter= @(x) x<= max(freq_range) & x>= min(freq_range);normalize =@(x) x/sum(x);
           ix =bandfilter(sig_stats.fq);
            peakFun =@(a) arrayfun(@(j) findpeaks(a(j,:), sig_stats.fq(ix),...
                'SortStr', 'descend', 'MinPeakProminence', 0.0055), 1:size(a,1), 'UniformOutput',false);
            [peaks,locs] = peakFun(sig_stats.psd(:,ix)) ; %peaks = psd, locs = frequency
            freq =zeros(size(peaks)); 
            for j = 1:numel(peaks)
                if numel(peaks{j}) > 1
                    % more than one peak within the range, take weighted
                    % sum of frequency 
                    wgts = normalize(peaks{j}); 
                    freq(j) = sum(locs{j}.*wgts); 
                elseif ~isempty(peaks{j})
                    freq(j) = locs{j}; 
                end
            end
            sig_stats.(stat) = freq';            
            
        case {'obw', 'oscbandwidth'}
            sig_stats.(stat)=obw(smoothData',Fs)';
        case 'max_pentropy'
            max_entropy= zeros(size(smoothData,1),1);    
            time_pts = max_entropy;
            
            for j =1:numel(max_entropy)
              [sig_entropy, tp]= pentropy(smoothData(j,:), Fs/3600,...
                  'FrequencyLimit',freq_range./3600);   
               [max_entropy(j), ix] = max(sig_entropy);
               time_pts(j) = tp(ix)/3600;
               
            end
            sig_stats.(stat) = max_entropy; 
        case 'noise_est'
            sig_stats.(stat) =  wnoisest(Data)';
    end
end
end