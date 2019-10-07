function peak2trough = get_peak2trough(trajectories)
% Calculates number of peaks in trajectories
%INPUT
%   'Data':     trajectories to use
%   
%OUTPUTs 
% number of peaks

% [metrics,~,~,~,~]=nfkbmetrics(id);
warning('defunct! use get_peak_stats'); 
% [~,pk_amps,~, valley_amps]=nfkbpeaks(trajectories);
% peak2trough =pk_amps(:,1:end-1) - valley_amps; 

end
