function  time_to_half_max_integral=get_time_to_half_max_integral(integrals)
 halfMaxIntegral = max(integrals,2)/2;
 

 distances = abs(integrals- halfMaxIntegral);
[row,frame] = find(distances ==min(distances,[],2));
[~,ix] = sort(row, 'ascend'); 
frame = frame(ix); 
 time_to_half_max_integral = (frame-1)/12;
end