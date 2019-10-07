function [HMW,Prominence, Amp, Locs] = halfMaxWidth(time_series, PeakHeight)
%--------------------------------------------------------------------------
%[HMW,Prominence, Amp, Locs] = halfMaxWidth(time_series, PeakHeight)
%--------------------------------------------------------------------------

% if isrow(PeakHeight)
%    PeakHeight = transpose( PeakHeight);
% end 7/10/2019

eqDims= size(time_series) == size(PeakHeight); 
if ~eqDims(1) && eqDims(2)
    %check that dimensions of time_series and PeakHeights are aligned
     PeakHeight = transpose( PeakHeight);
end
nPeaks=size(PeakHeight,2);
HMW = nan(size(PeakHeight,1), nPeaks);

Amp=HMW; Locs = HMW; Prominence = HMW; 
lastFrame =size(time_series,2);
FPH = 12;
t=linspace(0,(lastFrame-1)/FPH, lastFrame);


for i =1:size(time_series,1)
    minPeakHgt=min(PeakHeight(i,:)*0.95);
    if isnan(minPeakHgt); continue; end
     [Amp(i,:), Locs(i,:), HMW(i,:), Prominence(i,:)] =findpeaks(time_series(i,:), t,'NPeaks',nPeaks,...
        'MinPeakHeight', minPeakHgt,'WidthReference', 'halfheight', 'Annotate', 'extents', 'MinPeakDistance', 6/FPH);
    
end
end