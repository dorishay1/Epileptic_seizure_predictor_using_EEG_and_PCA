function [spec_entropy] = spectral_entropy(norm_power)
%this function computes the spectral entropy feature using the norm_power.

%simple calculation according to the formula.
log2_norm = log2(norm_power);
entropy_mat = norm_power.*log2_norm;
spec_entropy = -1 * sum(entropy_mat);


end

