function num_pks = get_num_pks(trajectories)
% Calculates number of peaks in trajectories
%INPUTs (optional)
%   'Data':     trajectories to use instead of generating new onex
%   
%OUTPUTs 
% number of peaks

% [metrics,~,~,~,~]=nfkbmetrics(id);

[pk_times,pk_amps,valley_times, valley_amps]=nfkbpeaks(trajectories);
ipt=diff(pk_times,1,2);

%%
ipt(ipt>35)=nan;
num_pks = sum(~isnan(ipt),2)+1;

end
