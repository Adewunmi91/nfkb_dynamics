function [frac_on, frac_osc,frac_trans, frac_perst]= get_osc_frac (metrx, cut_off_fq, tPerst, Baseline)
%%
%Calculates the fraction of cells that are on, the fractio of on cells that
%have oscillatory translocations, persistent translocations, and transient
%translocations. 
%[frac_on, frac_osc,frac_trans frac_perst]= get_osc_frac (metrx,  cut_off_fq,tOff, tPerst)
%INPUT: 
%       metrx:      struct of metrics & aux outputs from nfkbmetrics
%       cut_off_fq:     frequency above which a trajectory is considered
%       oscillatory
%       
%       tPerst      time above which a cell is dteremined to be persistent
%       Baseline    activity threshold to be considered activated   
%
% OUTPUTs:     
%       frac_on    fraction of cells above Baseline
%       frac_osc    fraction of cells with peak frequencies > cutt_off_fq
%       frac_trans   fraction of cells with 
%%











% %Initializations
% nCells=length(metrx);
% cell_frac = zeros(nCells,3);
% frac_trans = zeros([nCells 1]);
% %cut_off_fq = 0.42;
% %Baseline = 1.9; % Minimum activity required (99.9th pct confidence level over "off" set values)
% %tPerst = 2.75; % Cutoff time (envelope) for a cell to be declared "persistent"
% nfkb_lvl = find(abs(metrx(1).aux.thresholds - Baseline) ...
%                 == min(abs(metrx(1).aux.thresholds - Baseline)),1,'first');
%             
% for i = 1:numel(metrx)
%     cell_frac(i,1) = sum(metrx(i).metrics.off_times==0);
%     cell_frac(i,2) = sum((metrx(i).metrics.off_times>0)&(metrx(i).metrics.peakfreq<cut_off_fq));   
%     cell_frac(i,3) = sum(metrx(i).metrics.peakfreq>=cut_off_fq);    
%     frac_trans(i) = sum((metrx(i).metrics.off_times>0)&(metrx(i).metrics.peakfreq<cut_off_fq)&(metrx(i).metrics.envelope (:, nfkb_lvl) <tPerst))/sum (cell_frac(i,:));
%     cell_frac(i,:)= cell_frac(i,:)/sum(cell_frac(i,:));
%          
% end
% frac_on = 1-cell_frac(:,1);
% % percent "oscillatory" = peak frequeny > cutoff 
% frac_osc =  cell_frac(:,3)./frac_on;
% %percent transient = on AND not persistent AND not oscillatory
% frac_perst= 1-frac_osc-frac_trans(:);


%save('population_fraction','pop_frac');

%%
end
