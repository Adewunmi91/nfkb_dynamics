function [list, scalars,vectors] = get_peak_stat_list
% vectors=sort({'peak_times', 'valley_times','valley_amps', 'peak_amps', 'ipt','peak2trough', 'trough2peak'});
% scalars = sort({'median_ipt', 'mean_ipt', 'std_ipt', 'max_ipt', 'min_ipt', 'cv_ipt',...
% 'num_peaks', 'mean_peak_amp', 'median_peak_amp', 'std_peak_amp', 'cv_peak_amp',...
%  'max_peak2trough','mean_peak2trough', 'median_peak2trough', 'std_peak2trough', 'cv_peak2trough', ...
%  'max_trough2peak', 'mean_trough2peak', 'median_trough2peak', 'std_trough2peak', 'cv_trough2peak', 'pk1_whm', 'pk2_whm', 'pk1_prom', 'pk2_prom'});
% list = sort([vectors, scalars]); 
[list,scalars] = getFeatNames('peak_stats'); 
vectors = setdiff(list, scalars);
end