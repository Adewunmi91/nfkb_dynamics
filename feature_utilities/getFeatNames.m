function [FeatNames, Scalars, tbl] =getFeatNames(varargin)
% Returns the names of features
% -------------------------------------------------------------------------
% [FeatNames, Scalars] =getFeatNames(varargin)
% -------------------------------------------------------------------------
p = inputParser;
paths = loadpaths;

addOptional(p, 'Type', 'all', @istext);
parse(p, varargin{:});
featType = convertCharsToStrings(   p.Results.Type);
FeatNames = {}; Scalars = {};
fileName = fullfile(paths.analysis, 'dynamics','features.xlsx');
persistent sheets
persistent cache

tbl = table;
if ~isempty(featType)
    if isempty(sheets)
        
        [~,sheets] = xlsfinfo(fileName);
        
    end
    if strcmpi(featType, 'all')
        featType = convertCharsToStrings(sheets);
    end
    
    assert(all(ismember(featType, sheets)), 'feat type not found')
%     featType = cellstr(featType);
    if contains(lower(featType(1)), {'msmt', 'measurements', 'measurements', 'msmts'})
        FeatNames = get_msmt_list({'nfkbdim', 'morphology','nucintensity'});
    else
        if ~isempty(cache) && isfield(cache,featType(1))
           tbl= cache.(featType(1)).tbl;
           scalarTbl = cache.(featType(1)).scalarTbl;
        else
            
        
        tbl =sortrows( readtable(fileName, 'Basic', true, 'Sheet', featType(1), 'TextType', 'string'), 'Name', 'ascend');
        scalarTbl = filterByVar(tbl, "Category", "scalar");
        cache.(featType(1)).tbl = tbl;
        cache.(featType(1)).scalarTbl =scalarTbl;
        end
%         scalarTbl = filterTable(tbl, @(t) ismember(t, 'scalar'),{'Category'});
        Scalars = scalarTbl.Name;  FeatNames= tbl.Name;
    end
    
    if numel(featType) > 1
        [a, b,c] = getFeatNames(featType(2:end));
        FeatNames = sort(vertcat(FeatNames,a));
        Scalars = sort(vertcat(Scalars,b));
        tbl = vertcat(tbl, c);
    end
end
end