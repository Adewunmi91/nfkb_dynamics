%test peak prominence threshold
sig =struct; 

sig.ako = get_sig_stats(get_trajectories(596), 'Stats', {'power'});

sig.wt = get_sig_stats(get_trajectories(595), 'Stats', {'power'});
 
freq_range =[0.33, 1]; bandfilter= @(x) x<= max(freq_range) & x>= min(freq_range);normalize =@(x) x/sum(x);
zeroFrac = @(x) sum(x==0)/numel(x); 
%%
conds= fieldnames(sig);
thresh =linspace(0.001,0.02,100);
sep = zeros(size(thresh));

for k =1:1:numel(thresh)
    freqs = cell(2,1);
    for i =1:numel(conds)
        
      
        ix =bandfilter(sig.(conds{i}).fq);
        peakFun =@(a) arrayfun(@(j) findpeaks(a(j,:), sig.(conds{i}).fq(ix),...
            'SortStr', 'descend', 'MinPeakProminence', thresh(k)), 1:size(a,1), 'UniformOutput',false);
        [peaks,locs] = peakFun(sig.(conds{i}).power(:,ix)) ; %peaks = psd, locs = frequency
        freq =zeros(size(peaks));
        for j = 1:numel(peaks)
            if numel(peaks{j}) > 1
                % more than one peak within the range, take weighted
                % sum of frequency
                wgts = normalize(peaks{j});
                freq(j) = sum(locs{j}.*wgts);
            elseif ~isempty(peaks{j})
                freq(j) = locs{j};
            end
        end
        freqs{i} = freq';
    end
    sep(k) = abs(diff([zeroFrac(freqs{1}), zeroFrac(freqs{2})]));
end
results = table(thresh', sep', 'VariableNames', {'Threshold', 'Separation'}); 
results = sortrows(results, 'Separation', 'descend');
disp(results);