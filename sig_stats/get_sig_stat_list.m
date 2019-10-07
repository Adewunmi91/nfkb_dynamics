    function  [sig_stats, scalars]=get_sig_stat_list 
%     [sig_stats, scalars]=get_sig_stat_list 
% ------------------------------------------------------------------------
%returns field names of signal statistics 

% stat_list =sort({'power','bandpower','rms', 'noise_est', ...
%     'obw','mean_movmad','mean_movstd','mean_movvar', 'peak2rms', 'medfreq',...
%     'meanfreq', 'peak2peak', 'oscfreq', 'oscpower'});

% 
% tbl = readtable('features.xlsx', 'Basic', true, 'Sheet', 'sig_stats');
% stat_list = tbl.Name; 
% scalarTbl = filterTable(tbl, @(t) ismember(t, 'scalar'),{'Category'});
% scalarList = scalarTbl.Name;

[sig_stats, scalars] = getFeatNames('sig_stats');
end