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

features =struct;

sig_stat_list = get_sig_stat_list;
pk_stat_list = cellstr(get_peak_stat_list());


all_suffices = "_"+getTransformMethods();
hasSuffix = endsWith(featNames, all_suffices, 'IgnoreCase',true);
fullFeatNames = featNames;
suffices = strings(size(featNames));
tokens = chain(regexpi(featNames(hasSuffix), "^(\w+)_(\w+)$",'tokens'),...
    @unnest,@catElems, @collapseCellRows);
if any(hasSuffix)
[featNames(hasSuffix), suffices(hasSuffix)]=deal(tokens{:});
end
if any(ismember(featNames, pk_stat_list))
    pks =get_peak_stats(metrics); % avoids multiple calls inside the loop
end
for j=1:length (featNames)
    featName = featNames(j);
    suffix = suffices(j);
    % check to see if features struct has the features already calculated
    if ~isfield(features, featName)
        switch featName
            case {'derivative', 'integrals', 'time_series'}
                features.(featName) = metrics.(featName)(:,1: end);
            case {'median_derivative'}
                features.(featName) = nanmedian(metrics.derivatives,2);
            case {'mean_derivative'}
                features.(featName) = nanmean(metrics.derivatives,2);
            case pk_stat_list
                if contains(featName, 'ipt')
                    all_cells= nan(size(pks.kept, 1),size(pks.(featName),2));
                    all_cells(pks.kept,:)= pks.(featName);
                    features.(featName) = all_cells;
                else
                    features.(featName)=pks.(featName);
                end
            case {'pk2_ratio'}
                features.(featName)=metrics.pk2_amp./metrics.pk1_amp;
            case {'fold_change'}
                features.(featName) = get_fold_change(metrics.time_series);
            case {'time2HalfMaxIntegral'}
                endFrame = min(97, size(metrics.integrals,2));
                features.(featName) = get_time_to_half_max_integral(metrics.integrals(:,1:endFrame));
            case {'osc_cat', 'OscCat'}
                features.(featName)=osc_cats(metrics.peakfreq, metrics.off_times, 'cutoff_fq', 0.42);
            case 'last_falltime'
                features.(featName) = computeDuration(metrics.time_series, 1, false, 'Threshold',0.2);
            case 'max_pk1_speed'
                features.(featName)=get_max_pk1_speed(metrics);
            case cellstr(sig_stat_list)
                sig_stats = get_sig_stats(metrics.time_series,'stats', featName);
                features.(featName)=sig_stats.(featName);   
            otherwise
                data = copyFields(metrics, msmt);
                features.(featName) =data.(featName);    
        end
    end
    if ~isEmptyStr(suffices(j))
                 X = rmmissing(features.(featName),2);
%         X = features.(featName);
        N = size(X,2);
        M = size(X,1);
        
        d = min([ floor(sqrt(N)),5]);
        features.(fullFeatNames(j)) = nan(M,d);
        
        % check that the cov  matrix is positive definite since FA requires it
        if isempty(X)||(matches(suffix, ["FA", "factor", "factoran"])&& (min(svd(cov(rmmissing(X,2))))<= eps^(1/3)))
            continue
        end
        features.(fullFeatNames(j)) = transform_feats(X,suffix,d) ;
    end
    %         dExp = find(cumsum(normalize(s, 'norm',1))> 0.75,1);
    %         d= min([d,dExp]);   %         d = min([ceil((2-N)/((-2*N)+1)),5]);
end

