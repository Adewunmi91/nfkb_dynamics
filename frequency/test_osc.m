%new metric to decid osc category

aKO =get_features(596, {'osc_cat','oscCat'}); 
wt = get_features(595, {'osc_cat','oscCat'}); 
%% 


wtoscFrac =sum(ismember(wt.oscCat, 'osc'))/numel(wt.oscCat)
aKOoscFrac =sum(ismember(aKO.oscCat, 'osc'))/numel(aKO.oscCat)

wt_osc_frac = sum(ismember(wt.osc_cat, 'osc'))/numel(wt.osc_cat)
aKO_osc_frac = sum(ismember(aKO.osc_cat, 'osc'))/numel(aKO.osc_cat)