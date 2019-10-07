function features=computeFeatures(featNames,metrics, msmt)
% Computes features from metrics and msmt outputs of nfkbmetrics 
%--------------------------------------------------------------------------
% features=computeFeatures(featNames,metrics, msmt)
%--------------------------------------------------------------------------
% INPUT
%     featNames             cell array containing substet of get_feature_sets()
%     metrics               struct output from nfkbmetrics()
%     msmt                  'measure' output of nfkbmetrics
% OUTPUT
%     features              struct with fields(features) == featNames
    
featSet= get_feature_sets;
features =struct; 
if any(contains(featNames, featSet.peak_stats)); pks =get_peak_stats(metrics);  end 

for j=1:length (featNames)
    featName = featNames{j};
    switch featName
        case {'derivative', 'integrals', 'time_series'}
            features.(featName) = metrics.(featName)(:,1: end);
        case {'median_derivative'}
            features.(featName) = nanmedian(metrics.derivatives,2);
        case {'mean_derivative'}
            features.(featName) = nanmean(metrics.derivatives,2);
            
        case [featSet.msmt; featSet.spots]
            if nargin>2
%                 features.(featName)=msmt(:,featName);
                if isfield(msmt, featName)
                    features.(featName)=msmt.(featName);
                else
                    features.(featName)=[];
                end
            end
          
        case featSet.peak_stats          
%             pks =get_peak_stats(metrics);            
            if contains(featName, 'ipt')
                all_cells= nan(size(pks.kept, 1),size(pks.(featName),2));
                all_cells(pks.kept,:)= pks.(featName);
               features.(featName) = all_cells;
%             elseif contains(feature,'oscCat')
%                 data.(feature)=get_osc_cats(pks.num_peaks, metrics.off_times);
            else
               features.(featName)=pks.(featName);
            end
        case {'pk2_ratio'}
            features.(featName)=metrics.pk2_amp./metrics.pk1_amp;
% %         case {'pk1_whm'}
% %             interpeak_frames = ((metrics.pk2_time-metrics.pk1_time)*12)+1;
% %             features.(featName)=half_max_width(metrics.time_series, metrics.pk1_time,interpeak_frames);
        case {'fold_change'}            
            features.(featName) = get_fold_change(metrics.time_series);              
        case {'time2HalfMaxIntegral'}
            endFrame = min(97, size(metrics.integrals,2)); 
           features.(featName) = get_time_to_half_max_integral(metrics.integrals(:,1:endFrame)); 
        case {'osc_cat', 'OscCat'}
            features.(featName)=osc_cats(metrics.peakfreq, metrics.off_times, 'cutoff_fq', 0.42);
        
        case 'last_falltime'               
             features.(featName) = computeDuration(metrics.time_series, 2); 
        case 'max_pk1_speed'
            features.(featName)=get_max_pk1_speed(metrics);
        case featSet.signal
            sig_stats = get_sig_stats(metrics.time_series,'stats', {featName});
            features.(featName)=sig_stats.(featName);
            
        otherwise
            features.(featName) =metrics.(featName);
            
    end

end 