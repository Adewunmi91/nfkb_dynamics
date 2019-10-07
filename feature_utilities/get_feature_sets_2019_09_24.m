function  sets = get_feature_sets(varargin)
% Provides names of feature sets
%-------------------------------------------------------------------------
%  sets= get_feature_sets ('peaks');
%--------------------------------------------------------------------------
% INPUT
%   SETS : ["codons", "info","signal", "violin", "transform", "scalar",
%   "vector", "peaks", "stats", "nfkb", "legacy", "metric", "morphology"]
% OUTPUT
%   sets:        iscell if nargin =1 else isstruct
f=@(x) toCol(catElems( arrayfun(@get_feature_sets,x,'UniformOutput',false)));
switch nargin
    case 0
%         X = {'all'};
        X= "all";
    case 1
        X = convertCharsToStrings(varargin{:});
%         if ischar(varargin{1})
%             X = cellstr(varargin{:});
%         else
%             X = varargin{:};
%         end
    otherwise
%         X = varargin;
         X = convertCharsToStrings(varargin);
end
% if isrow(X); X= X';end
X= toCol(X);
Y =cell( size(X));
tfMeth=cellstr(lower(getTransformMethods())) ;
for i=1:numel(X)
    x =X{i};
    switch lower(x)
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
            Y{i} = "OscCat";
        case {'all'}
            Y{i} = [get_metrx_list(); get_msmt_list(); get_peak_stat_list();...
                get_derived_metrx_list(); get_sig_stat_list(); get_spot_msmts(); {'OscCat'}' ;...
                get_feature_sets('transform')];
        case "dynamics"
%             Y{i} = toCol(catElems( arrayfun(@get_feature_sets,...
%                 ["metrics","peaks","derived", "group","transform","signal"], 'UniformOutput',false)));
             Y{i} = f(["metrics","peaks","derived", "group", "transform", "signal"]);
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
        case {'transform'}
            
            Y{i} = mat2vec(getVectorMetrics()+ "_"+getTransformMethods());
        case tfMeth
            Y{i} = getVectorMetrics()+ "_" + x;
        otherwise
            Y{i} = X(i);
    end
    
end
    
% sets = unique(vertcat(Y(:)), 'sorted');
sets = unique(convertContainedStrings(Y),'stable');
% sets = chain(Y,  @flatten,@unique, @flatten, @string);
    
end

function vecMetrics = getVectorMetrics()
% [~,~,tbl] = getFeatNames(["metrics","peak_stats", "sig_stats"]);
% [~,~,tbl] = getFeatNames("metrics");
%             vecT = filterByVar(tbl, 'Category', "vector");
%             vecMetrics = vecT.Name;
vecMetrics = toCol(["time_series", "oscfrac", "derivatives", "duration", "envelope"]);
end