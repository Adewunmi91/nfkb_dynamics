function output= get_peak_features (trajectories, varargin)

p=inputParser;
addRequired(p,'trajectories', @isnumeric);
addParameter (p,'min_pks',2,@isnumeric);
addParameter (p,'max_pk_diff',35,@isnumeric);
parse (p,trajectories,varargin{:});
trajectories = p.Results.trajectories;
min_pks =p.Results.min_pks;
max_pk_diff=p.Results.max_pk_diff;
%[metrics,~,~,info,~]=nfkbmetrics(id);

[output.peak_times,output.peak_amps, output.valley_times, output.valley_amps]=nfkbpeaks(trajectories);
ipt=diff(output.peak_times,1,2);

%%
ipt(ipt>max_pk_diff)=nan;

tot_pks = sum(~isnan(ipt),2)+1;
kept = ones([size(ipt,1),1]);

ipt(tot_pks<min_pks,:) = [];
kept(tot_pks<min_pks)=0;
ipt = ipt.*5;%convert to minutes
%Brooks' analysis excludes ipt > 35, excludes cells with <5peaks

output.mean_ipt=nanmean(ipt,2);
output.med_ipt=nanmedian(ipt,2);
output.std_ipt=nanstd(ipt,[],2);
output.var_ipt =nanvar(ipt,[],2);
output.max_ipt = nanmax(ipt,[],2);
output.min_ipt = nanmin(ipt,[],2);
output.cv_ipt= output.std_ipt./output.mean_ipt;
output.ipt=ipt;
output.num_peaks=tot_pks;
%%
output.peak2trough = [output.peak_amps(1:end-1)-output.valley_amps,...
    output.peak_amps(2:end)-output.valley_amps];
output.max_peak2trough=max(output.peak2trough, [],2);
output.kept = logical(kept);

end
