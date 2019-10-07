function output= get_ipt (trajectories, varargin)

p=inputParser;
addRequired(p,'trajectories', @isnumeric);
addParameter (p,'min_pks',2,@isnumeric);
addParameter (p,'max_pk_diff',35,@isnumeric);
parse (p,trajectories,varargin{:});
trajectories = p.Results.trajectories;
min_pks =p.Results.min_pks;
max_pk_diff=p.Results.max_pk_diff;
%[metrics,~,~,info,~]=nfkbmetrics(id);

[pk_times,pk_amps, valley_times, valley_amps]=nfkbpeaks(trajectories);
ipt=diff(pk_times,1,2);

%%
ipt(ipt>max_pk_diff)=nan;

tot_pks = sum(~isnan(ipt),2)+1;
kept = ones([size(ipt,1),1]);

ipt(tot_pks<min_pks,:) = [];
kept(tot_pks<min_pks)=0;
ipt = ipt.*5;%convert to minutes
%Brooks' analysis excludes ipt > 35, excludes cells with <5peaks
% mean_ipt=nan(size(ipt,1),1);
% std_ipt=nan(size(ipt,1),1);
% med_ipt=nan(size(ipt,1),1);
% mode_ipt=nan(size(ipt,1),1);

%% Summary statistics for each cell's interpeak time
% for i=1:size(ipt,1)% for each cell
%     mean_ipt(i)=nanmean(ipt(i,:));
%     med_ipt(i)=nanmedian(ipt(i,:));
%     std_ipt(i)=nanstd(ipt(i,:));    
%     mode_ipt (i)=mode(ipt(i,:));
%    
% end
    mean_ipt=nanmean(ipt,2);
    med_ipt=nanmedian(ipt,2);
    std_ipt=nanstd(ipt,[],2);    
    var_ipt =nanvar(ipt,[],2); 
    max_ipt = nanmax(ipt,[],2);
    min_ipt = nanmin(ipt,[],2); 
    
%% %Intracellular variance
%Ways of looking at variability
% 1. cv within a cell's trajectories
%       calculate mean,std of the c.v.,  
intra_cv= std_ipt./mean_ipt;
mean_cv=nanmean(intra_cv);
%med_cv=nanmedian(intra_cv);
%mode_cv=mode(intra_cv);
%std_var=nanstd(intra_cv);
%% Interceullular variance
% 2. mean, median, mode of cell's trajectory
%       cv of mean, median, mode     

cv_mean=nanstd(mean_ipt)/nanmean(mean_ipt);
%cv_median=nanstd(med_ipt)/nanmean(med_ipt);
%cv_mode=nanstd(mode_ipt)/nanmean(mode_ipt);
%%
cv =nanstd(ipt(:))/nanmean(ipt(:));
output.std_ipt =nanstd(ipt(:)); 

output.pks=tot_pks;
output.ipt=ipt;
output.mean_ipt=mean_ipt; 
output.median_ipt=med_ipt;
output.cv_all=cv;
output.cv_cell=intra_cv;
output.mean_cv=mean_cv;
output.cv_mean=cv_mean;
output.var_ipt=var_ipt; 
output.max_ipt = max_ipt; 
output.min_ipt = min_ipt; 
output.pk_amps=pk_amps;
output.valley_times=valley_times;
output.valley_amps = valley_amps;
output.kept = logical(kept); 

end
