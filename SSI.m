function [SSI_mat] = SSI(power,f)
%This function calculates the spectral slope and intercept for each window after ln
%transformation using the function polyfit.possible plot is avialble in comment
%the end of the function.

%transform to ln
ln_power = log(power);
ln_f = log(f);

[~,n_t_windows]= size(power);

%creating memory
SSI_mat = zeros(2,n_t_windows);

%polyfit for each window and storing it in the mat. the 1st row is for slope and
%the 2nd for intercept. each column is for a different window.
for i = 1:n_t_windows
    SSI_mat(:,i) = polyfit(ln_f',ln_power(:,i),1);
end

end

