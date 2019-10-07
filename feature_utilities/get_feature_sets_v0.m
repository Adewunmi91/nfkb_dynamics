function  sets = get_feature_sets(varargin)
% Provides names of feature sets
%-------------------------------------------------------------------------
%  sets= get_feature_sets ('peaks');
%--------------------------------------------------------------------------
%   OUTPUT
%   sets:        iscell if nargin =1 else isstruct

switch nargin
    case 0
        X = {'all'};
    case 1
        if ischar(varargin{1})
            X = cellstr(varargin{:});
        else
            X = varargin{:};
        end
    otherwise
        X = varargin;
end
if isrow(X); X= X';end
Y =cell( size(X));

for i=1:numel(X)
    x =lower(X{i});
    switch x
        case {'codon', 'codons'}
            Y{i} = getCodonNames();
        case {'info', 'cc', 'mi', 'it'}
            Y{i} = [get_metrx_list(); get_peak_stat_list(); {'peak2peak'; 'peak2rms'; 'oscbandwidth';  'oscpower'}; get_derived_metrx_list()];
        case    {'msmt', 'msmts', 'measurements', 'measurement'}
            Y{i} = get_msmt_list();
        case {'morphology'}
            Y{i} = get_msmt_list({'morphology'});
        case {'signal', 'sig_stat', 'sig_stats'}
            Y{i} = get_sig_stat_list();
        case { 'stattests', 'stat','stattest', 'statistics','stats'}
%             Y{i} = get_scalar_metrics();
            Y{i}=setdiff(pick(@getFeatNames,2),{'osc_cat', 'OscCat'});
        case {'scalar'}
            [~, Y{i}] =getFeatNames();
        case {'violin', 'viol'}
            [~, metrx ]= get_metrx_list(); [~, derivedMetrx] = get_derived_metrx_list(); 
            [~, pkStats] = get_peak_stat_list(); 
            pkStats = pkStats([startsWith(pkStats, 'pk') | startsWith(pkStats,'median') | endsWith(pkStats, 'peaks')]);
            Y{i} =setdiff([pkStats; metrx; derivedMetrx; {'oscpower'; 'oscbandwidth'; 'oscfreq'}], 'osc_cat');
        case {'original_metrics'}
            Y{i} =setdiff(get_metrx_list(), {'intwin1'; 'intwin3'});
        case {'derived'}
            Y{i} = get_derived_metrx_list();
        case 'spots'
            Y{i} = get_spot_msmts();
        case {'grp', 'group'}
            Y{i} = {'OscCat'};
        case {'all'}
            Y{i} = [get_metrx_list(); get_msmt_list(); get_peak_stat_list();...
                get_derived_metrx_list(); get_sig_stat_list(); get_spot_msmts(); {'OscCat'}];
        case {'metrics', 'metric', 'metrx', 'metrix'}
            Y{i} = get_metrx_list();
        case {'vector'}
            Y{i} = setdiff([get_metrx_list(); get_msmt_list(); get_peak_stat_list(); ...
                get_derived_metrx_list(); get_sig_stat_list(); get_spot_msmts()], get_scalar_metrics());
        case {'peaks'}
           Y{i} = [get_peak_stat_list;{  'peak2peak'; 'peak2rms'; 'pk2_ratio'}];
        case {'test'}
            Y{i} = [{'osc_cat';'time_series'};...
                intersect(get_metrx_list(),get_scalar_metrics())];
        case {'legacy'}
            Y{i} = [get_metrx_list(); get_msmt_list()];
        case {'nfkb'}
            
             Y{i} = [get_metrx_list(); get_peak_stat_list();...
                get_derived_metrx_list(); get_sig_stat_list()];
        otherwise
            Y{i} = X{i};
    end
    
end
    
% sets = unique(vertcat(Y(:)), 'sorted');
sets = chain(Y,  @flatten,@unique, @flatten);
    
    %
    %
    % if nargin > 1
    %     name =varargin{1};
    %     sets=   get_feature_sets(varargin{2:end});
    %
    % elseif nargin ==1
    %     name = varargin{:};
    %     sets ={};
    % else
    %     name ='';
    %     sets=struct;
    % end
    %
    % s = cell(4,1);
    % msmts= get_msmt_list();
    % [metrx, s{1}]=get_metrx_list();
    % [derived_features, s{2}] = get_derived_metrx_list();
    % [sig_stats_list, s{3}] = get_sig_stat_list();
    % [peak_stats_list, s{4}]= get_peak_stat_list();
    % spot_msmts = get_spot_msmts();
    % codonNames = getCodonNames();
    % scalar_metrics = vertcat(s{:});
    % featNames.('group') = { 'OscCat'};
    % featNames.('codon') = codonNames;
    % featNames.('all')=[ metrx; derived_features; msmts; sig_stats_list; peak_stats_list; spot_msmts];
    % featNames.('peak_stats')= peak_stats_list;
    % featNames.('info')= [metrx; peak_stats_list; {'peak2peak'; 'peak2rms'; 'obw'}; derived_features];
    % featNames.('nfkb') = [metrx; peak_stats_list; derived_features; setdiff(sig_stats_list, {'noise_est', 'power'})];
    % featNames.('legacy') = [msmts; metrx];
    % featNames.('scalar') = scalar_metrics;
    % featNames.('peaks') =sort([ peak_stats_list; {'peak2peak'; 'peak2rms'; 'pk2_ratio'}]);
    % featNames.('core')=sort([ msmts; metrx;derived_features;peak_stats_list ;{'peak2peak'; 'peak2rms'}]);
    % featNames.('msmt')= msmts;
    % featNames.('spots')=spot_msmts;
    % featNames.('signal') =  sig_stats_list;
    % featNames.('metrics')= metrx;
    % featNames.('orig_metrics')= setdiff(metrx, {'intwin1'; 'intwin3'});
    % featNames.('derived')= derived_features;
    %
    % featNames.('test')=[{'osc_cat';'time_series'};intersect(metrx,scalar_metrics)];
    % featNames.('statTests') =setdiff(scalar_metrics,[ {'OscCat'; 'osc_cat'}; scalar_metrics(contains(scalar_metrics,{'min'; 'std'}))]);
    % %  featNames.('statTests')=setdiff(scalar_metrics,{ 'OscCat'});
    % featNames.('all')= sort([ metrx; derived_features; msmts; sig_stats_list; peak_stats_list; spot_msmts; featNames.group]);
    %
    % featNames.('vector')= setdiff(featNames.all, [featNames.group; featNames.scalar]);
    % featNames.('violin') = {'off_times', 'peakfreq', 'max_amplitude','max_integral',...
    %     'num_peaks', 'pk1_time', 'pk1_amp', 'median_derivative', 'pk1_prom',...
    %     'pk1_whm', 'median_ipt', 'median_derivative', 'obw', 'fold_change', ...
    %     'median_trough2peak', 'median_peak2trough', 'peak2rms', 'time2HalfMaxIntegral'}';
    % if isempty(name)
    %     sets=featNames;
    % else
    %     switch(lower(name)) %aliases for feature sets
    %         case lower(fieldnames(featNames))
    %             sets=[sets;featNames.(name)];
    %         case { 'mi', 'i', 'cc'}
    %             sets=[sets; featNames.('info')];
    %         case { 'metrx'}
    %             sets = [sets;featNames.('metrics')];
    %         case {'sig_stats'}
    %             sets=[sets;featNames.('signal')];
    %         case { 'derived_metrics', 'derived_metrx', 'derived_features'}
    %             sets=[sets;featNames.('derived')];
    %         case {'msmts', 'measurement', 'measurements'}
    %             sets=[sets;featNames.('msmt')];
    %         case {'spots', 'spot'}
    %             sets=[sets;featNames.('spots')];
    %         case {'stat', 'stats', 'statistics'}
    %             sets=[sets;featNames.('statTests')];
    %         case {'viol', 'violin'}
    %             sets=[sets;featNames.('violin'), 'time_series'];
    %         otherwise
    %             sets=[sets;{name}];
    %     end
    %
    % end
    
end