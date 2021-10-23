function [norm_power] = normalize_power(power)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
        total_power = sum(power);
        norm_power = power.*(1./total_power);

end

