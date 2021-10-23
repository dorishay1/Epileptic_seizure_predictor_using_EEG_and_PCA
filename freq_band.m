function [bands_map,f,n_freq_bands] = freq_band(Names,Band,resulotion)
%this function acceptes matching names and band vectors and erange them.
%the main cell will include the name and the range in a cell.
%the 3rd row in the cell contains the index in the f vector for the freq.
%this function also preduce frequancy vector according to resulution
%and the number of freq bands.

%first the function checks if there is no gaps between the bands.
for i = 1:(length(Band)-1)
    if Band(i,2) ~= Band(i+1,1)
        error('Bands illegal')
    end
end

%producing freq vector using the min and max freq and resulution.
max_freq = max(max(Band));
min_freq = min(min(Band));
f = min_freq:resulotion:max_freq;

%creating memory for the cell.
bands_map = cell(3, length(Names));

for N_Band = 1: length(Names)
    bands_map{1,N_Band} = Names(N_Band);
    bands_map{2,N_Band} = Band(N_Band, :);
    
    %matching the correct index from f for later usage.
    %each vector starts from the lower limit and gets until the upper limit
    %index in f so there is no overlap calculations.
    bands_map{3,N_Band} = (find(f == Band(N_Band,1))):((find(f == Band(N_Band,2))-1));
end

n_freq_bands = length(bands_map);

end

