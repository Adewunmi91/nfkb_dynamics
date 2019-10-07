function width = half_max_width (time_series, pk_loc, wdw)
%INPUTS:
%time_series: N x T matrix, N = number of cells, T = number of frames
%pk_loc: N x 1 matrix, peak times
%wdw: how many frames forward or backward

idx= ~(isnan(pk_loc));
len = size (time_series, 1);
width = zeros (len, 1);
%Find index of the pk1

for i= 1:len 
    if idx(i)
    loc = uint16((pk_loc (i)*12)+1);
    %convert loc to an index
      
    half_max = time_series (i, loc)./2;
    win = uint16(wdw (i));
    start = max (1, loc - win);   
    left_diff = abs (time_series(i, start:loc) - half_max); %absolute diff
    [~,left]= min (left_diff);%index of the closes to half-max
    left = left+start-1;
    %right side
    the_end = min (loc+win, size (time_series,2));
    right_diff = abs (time_series (i,loc: the_end)- half_max);
    [~,right]= min (right_diff);
    right = right+loc-1;
    width (i) =(right-left)/12;
    end
    
end


end