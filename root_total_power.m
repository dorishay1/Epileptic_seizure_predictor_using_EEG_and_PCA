function [root_total_power_mat] = root_total_power(power)
%this function calculates the root total power for the power mat.

%summing up the power for each window.
total_power = sum(power);

%square root for the summing.
root_total_power_mat = sqrt(total_power);
end

