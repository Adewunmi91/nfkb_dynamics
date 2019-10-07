function oscFrac = getOscFrac(ID,varargin)
% Gets oscillatory fraction of experiment ID
%--------------------------------------------------------------------------
% oscFrac =getOscFrac(ID,GrpVar)
%--------------------------------------------------------------------------
%   INPUT
%      REQUIRED
%                ID:      experiment ID
%      OPTIONAL
%               'GrpVar'    {'ID', 'XYPos'}
%               'Data':     Data Table 
%   OUTPUT
%% Parse inputs

p = inputParser;
% [p,defaults] = addParamDefaults(p, 'nfkbmetrics'); 
[p,~,defaults] = addParamDefaults(p, 'nfkbmetrics'); 
[p,~,defaults2] = addParamDefaults(p, 'get_feature_tbl');
addRequired(p, 'ID', @isnumeric);
addParameter(p,'GrpVar', 'ID');
addParameter(p, 'Data', []);
parse(p, ID, varargin{:});
ID= p.Results.ID; GrpVar = p.Results.GrpVar;

%% Load table of NFkB oscillation categories {'osc', 'non_osc'} and count categories
if isempty(p.Results.Data)
    
    args= getNameValuePairs(p.Results, union(defaults, defaults2)); 
    args = updateProps(args, {'Features', "OscCat"});
    Data = get_feature_tbl(ID,  args{:});
else
    Data = p.Results.Data;
end

countCats = @(a) array2table(transpose(countcats(a)), 'VariableNames',transpose(categories(a)));

[g,gL]=findgroups(Data(:,GrpVar));
oscVar = intersect(Data.Properties.VariableNames, {'OscCat', 'osc_cat'});
    oscFrac =[gL, splitapply(countCats, Data.(oscVar{1}),g)];

% [g,~,gL] =grp2idx(Data{:,GrpVar});

% oscFrac =[array2table(gL, 'VariableNames', {GrpVar}),...
%     splitapply(countCats, Data.OscCat,g)];

oscFrac.on = oscFrac.osc+oscFrac.non_osc;
oscFrac.osc_ratio=round(oscFrac.osc./oscFrac.non_osc,1);
total=oscFrac.on+oscFrac.off;
oscFrac.osc= round(oscFrac.osc./oscFrac.on,2)*100;
oscFrac.non_osc = round(oscFrac.non_osc./oscFrac.on,2)*100;
oscFrac.on =round(oscFrac.on./total,2)*100;
oscFrac.off=round(oscFrac.off./total,2)*100;
oscFrac=movevars(oscFrac,{'on'}, 'Before',{'osc'});

%Ensures that oscFrac table maintains the same order as the ID
oscFrac = sortByVar(oscFrac,'ID',ID);

% cats=arrayfun(@(x) struct2array(get_Features(x,{catName})),ID, 'UniformOutput', false);
% counts = arrayfun(@(x) array2table(transpose(countcats(x{:})),....
%     'VariableNames',transpose(categories(x{:}))) ,...
%     cats,'UniformOutput',false);
%
% oscFrac=vertcat(counts{:});
% oscFrac=[array2table(ID','VariableNames',{'ID'}),oscFrac];
% %
% % counts =array2table(countcats(cats)', 'VariableNames', categories(cats));
% % counts.on = numel(cats) -counts.off;
% % counts.on = sum(~(off_times==0));
% % var_names = {'on', 'osc', 'non_osc'};
% % frac = counts(:,var_names);
% % frac = counts(:,2:end);
% % frac.on = round(counts.on/numel(cats),2);
% % frac{:,var_names(3:4)} = round(counts{:,var_names(3:4)}./counts.on,2);
end