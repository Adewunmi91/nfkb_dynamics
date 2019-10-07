function frac =get_frac(ID, varargin)
%Returns fraction of cells that are oscillatory
warning('Deprecated, use getOscFrac()');
% p=inputParser; 
% addRequired(p, 'ID', @isnumeric);
% addParameter(p, 'cutoff_fq', 0.42, @isnumeric); 
% addParameter(p, 'cutoff_time',2.75, @isnumeric); 
% addParameter(p, 'Params', []);
% parse(p,ID, varargin{:});
% ID = p.Results.ID; 
% 
% % Line plot of fraction of cells that are active
% all_fracs=cell(1,length(ID));
% cut_off_fq = p.Results.cutoff_fq;
% 
% for i=1:length(ID)
%     if isempty(p.Results.Params)
%         metrics=get_feature_tbl(ID(i));
%     else
%         metrics=get_feature_tbl(ID(i), p.Results.Params);
%     end
% [ ~,all_fracs{i}] = osc_cats(metrics.peakfreq,metrics.off_times, 'cutoff_fq', cut_off_fq);
% end
% frac = cat(1,all_fracs{:});
% end


% tPerst = p.Results.cutoff_time; % Cutoff time (envelope) for a cell to be declared "persistent"


%Baseline = 1.9; % Minimum activity required (99.9th pct confidence level over "off" set values)

% Baseline = info.Baseline;
% nfkb_lvl = find(abs(aux.thresholds - Baseline) ...
%                 == min(abs(aux.thresholds - Baseline)),1,'first');

%Initializations

% cell_frac = [];
% frac_trans = [];




%     cell_frac(1) = sum(metrics.off_times==0);
%     cell_frac(2) = sum((metrics.off_times>0)&(metrics.peakfreq<cut_off_fq));   
%     cell_frac(3) = sum(metrics.peakfreq>=cut_off_fq);    
%     frac_trans = sum((metrics.off_times>0)&(metrics.peakfreq<cut_off_fq)&(metrics.envelope (:, nfkb_lvl) <tPerst))/sum (cell_frac(:));
%     cell_frac(:)= cell_frac(:)/sum(cell_frac(:));
%          
% 
% frac_on = 1-cell_frac(:,1);
% % percent "oscillatory" = peak frequeny > cutoff 
% frac_osc =  cell_frac(:,3)./frac_on;
% %percent transient = on AND not persistent AND not oscillatory
% frac_perst= 1-frac_osc-frac_trans(:);

% all_fracs{i} =table(round(frac_on*100,1), round(frac_osc*100,1),round(frac_trans*100,1), round(frac_perst*100,1),'VariableNames',varNames);
%save('population_fraction','pop_frac');

