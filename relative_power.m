function [relative_power_per_window] = relative_power(power,freq_map,n_windows)
%This function takes the power matrix and converts the power to relative power
%according to bands in the settings.

n_freq_bands = length(freq_map);

total_power = sum(power);

%creating memory
relative_mat = zeros(n_freq_bands,n_windows);

%for each band, the loop sums only the relevent band width and then multiplying
%each value in 1/the total power that relvent for the band.
for band = 1:n_freq_bands
    
    %we take only the relvent spectrum using the index vectors in 'freq_map' and
    %sum it.
    relative_mat(band,:) = sum(power(freq_map{3,band},:));
    
    %then we split by the total power fot the selected band width.
    relative_mat(band,:) = relative_mat(band,:).*(1./total_power);
end

relative_power_per_window = relative_mat;

end

