function [Fedge] = Sedge_freq(norm_power,f,precentile)
%this function estimate the specific frequancy that is the edge
%between the inserted precentile (we have made here extra effort to be accurate
%and efficent - not using loops).
%*we assume that the small gap between each 2 freq (that is a result of the resulution)
%is linear.
%the general idea is to find the disterbutions that are the closest to the
%precentile and estimate the specific edge using them.

%extracting the resultion from f.
f_res = diff(f);
f_res = f_res(1);

%preparing cumulative distrbution for each window.
cum_dis = cumsum(norm_power);

%finds only the index relative to x & y that are above the precentile in each
%window (in all the matrix)
[xf, yf] = find(cum_dis>precentile);

%extracting the index in yf that is the first for each window.(i.e: the 2nd
%element in yf_cor represent the place of the first element in the 2nd window
%that is above the precentile).
[~, yf_cor] = unique(yf);

%extracting the index vector for the rows that relevent only to the first
%element that is above the precentile (and then matching vector for the columns).
x_cum_dis_idx = xf(yf_cor);
y_cum_dis_idx = yf(yf_cor);

%in order to get the disterbutions, we need to get the info from cum_dis and
%call it by matrix indexing. 'sun2ind' converts the places vectors we found to
%one index vector that fits the way matlab treats matrixs.
v_idx = sub2ind(size(cum_dis),x_cum_dis_idx,y_cum_dis_idx);

%because we want to be as accurate as we can and the disterbutions that we get
%are not exactly the wanted precentile (due to limited resulution) we will
%estimate the edge using the nearest freq.

%upper limit disterbution value will be the first element in each row that is
%above the precentile and the doun limit will be one element before it.
up_lim_dis = cum_dis(v_idx);
down_lim_dis = cum_dis(v_idx-1);

%we want to get the proportion of the step, how far is the down limit dis from
%the actual precentile and from the up limit. delta will be this proportion.
diff_dis = up_lim_dis-down_lim_dis;
dif_down_precentile = precentile - down_lim_dis;

delta = dif_down_precentile./diff_dis;

%getting the freq(Hz) of the down limit.
down_lim_f = f(x_cum_dis_idx-1);

%we are taking a step that is proportional to the disterbutions that calculated
%and the resulotion of f.
Fedge = down_lim_f+(delta*f_res)';


end

