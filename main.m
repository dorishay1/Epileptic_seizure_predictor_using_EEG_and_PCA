% Ex - 5

% Dor Ishay         312328339
% Liav Sommerfeld   312481468

%% Data handling
clear; close all; clc;

% cell array for handling the data.

%main_data is the primay cell
%nrow - 3 (Patient_n, Seisure_n, data_mat), ncol - index files (1:n_files).
%more details inside the function.
%n_measures - the number of different measures.
%n_elect - nunber of electrodes in total for each measure.
%t_bins - the number of the sampels in the data (length of the mat).
[main_data,n_measures,n_elect,t_bins] = data2cell('Data.zip');
%% General Settings

fs = 250;                       %sampling frequency [Hz]
t_window = 40*fs;               %t = time, length of time windows [sec]
t_overlap = t_window/2;         %overlap between windows [sec]/[as part of window]
f_resulotion = 0.1;             %frequancy vector resulution [Hz].

%frequancy bands

%this 2 vectors can be modified. each index in the name vectors is connected
%later with the same range index. more info inside 'freq_band'.
Bands_Name = ["delta" , "theta" , "low_alpha", "high_alpha", "beta", "gamma"];
Bands_range = [ 1, 4.5;  4.5,8;   8, 11.5;   11.5, 15;  15,30 ;  30, 40];

[freq_map,f,n_freq_bands] = freq_band(Bands_Name,Bands_range,f_resulotion);

%features settings

%pwelch
p_window = 2*fs;                    %pwelch window size[sec]
p_overlap = p_window/2;             %pwelch overlap between windows [sec]/[as part of window]

%spectral edge
S_edge_precntile = 0.9;             %the precentile for spectral edge calculations.

%extra analyze (bonus hunters)
Extra_analysis = 1;        %for extra analysis 'Extra_analysis' need to be equal to 1.
PCA_max_dim = 5;

%% The main loop
%the main loop is working per measure. for each measure the code analyse the
%data from all the elctrodes, calculate PCA & and ploting the results.

for n_file = 1:n_measures
    
    current_data = main_data{3,n_file};         %loading the current measurement.
    
    for idx_elect = 1:n_elect
        
        cur_elect = current_data(idx_elect,:);  %loading the currnet electrode.
        
        %spliting the original data into windos in the length of 't_window' with
        %regarding the 't_overlap' settings.for making sure that we take only
        %the full windows (without zero padding) - 'nodelay' makes sure that the zero
        %padding is in the end and taking only the full matrix without the
        %leftovers (including the zeros) is with saving only the first argument.
        %*each column will represent a window.
        [windows_mat, ~] = buffer(cur_elect,t_window,t_overlap,'nodelay');
        
        %preparing memory - there will be matrix for each electrode, each row
        %represent different featue and each column is different time window. We
        %will do it only once (only if it is the first electrode). we do it here
        %because we must use buffer in order to know the number of windows.
        if idx_elect == 1
            [~, n_t_windows] = size(windows_mat);
            data = zeros((2*n_freq_bands + 6),n_t_windows,n_elect);
        end
        
        
        %calculate power spectra for each window using pwelch.
        power = pwelch(windows_mat,p_window,p_overlap,f,fs);
        
        norm_power = normalize_power(power);        %normalizing the power matrix.
        
        %next we calculate the different features for each elctrode usually
        %using a uniqe function (more info inside the functions).
        
        %relative power(1-6)
        data((1:n_freq_bands),:,idx_elect) = relative_power(power,freq_map,n_t_windows);
        
        %log relative power (7-12)
        
        %we convert the power to log power using the following trick(more
        %explantion in the report) and use the 'relative_power' function to
        %calculate the relative log power.
        log_power = log(exp(1).*power./min(power));
        data((n_freq_bands+1:2*n_freq_bands),:,idx_elect) = ...
            relative_power(log_power,freq_map,n_t_windows);
        
        %root total power (13)
        data((2*n_freq_bands+1),:,idx_elect) = root_total_power(power);
        
        %spectral slope & intercept (14-15)
        data(((2*n_freq_bands+2):(2*n_freq_bands+3)),:,idx_elect) = SSI(power,f);
        
        %spectral moment (16)
        %simple matrix multiplication between the f vector and norm_power fits
        %the formula for spectral moment.
        data((2*n_freq_bands+4),:,idx_elect) = f*norm_power;
        
        %spectral edge (17)
        data((2*n_freq_bands+5),:,idx_elect) = Sedge_freq(norm_power,f,S_edge_precntile);
        
        %spectral entropy (18)
        data((2*n_freq_bands+6),:,idx_elect) = spectral_entropy(norm_power);
        
    end
    
    
    %% reshpae & z-score
    %after the computation for all the features for all the electrodes per
    %measurment, we reerange the data so we can easily do it PCA.
    
    
    data = permute(data,[1 3 2]);
    data= reshape(data, [], n_t_windows);
    data = zscore(data,0,2);
    
    
    %% PCA
    %standert PCA calculations. we skipped the part of substracing the mean
    %example because the PCA is for z scores and such as the mean is 0.
    C = (data*data')./(n_t_windows-1);
    
    %taking only the first 3 eigenvectors (most covariance).
    %usings eigs we can choose to compute and present only the number of
    %eigenvectors that we want and by that saving unneccesery computaion time.
    [U, D] = eigs(C,3);
    
    PCA_dim = (U')*data;
    
    %% Plot
    
    %calculations of the experiment time according to the sample freq:
    %total experiment time in minutes and time vec thet goes down from the
    %total time of the exp until zero.
    exp_time = (t_bins/fs)/60;                                  %[minutes]
    Time_vec = (-exp_time:exp_time/(n_t_windows-1):0);          %[minutes]
    
    
    %main plots
    fig = figure('Units','normalized','Position', [0 0 1 1]);
    hold on;
    
    sgtitle(['Patient - ',main_data{1,n_file},'  ,Seizure-', main_data{2,n_file}], 'FontSize', 20);
    
    % 2D plot
    subplot(1,2,1);
    scatter(PCA_dim(1,:),PCA_dim(2,:),8,Time_vec,'filled');
    xlabel('PC1','FontSize', 16);ylabel('PC2','FontSize', 16);
    title('2D projection', 'FontSize', 16);
    pbaspect ([1,1,1])
    set(gca,'FontSize',14)
    
    % 3D plot
    subplot(1,2,2);
    scatter3(PCA_dim(1,:),PCA_dim(2,:),PCA_dim(3,:),8,Time_vec,'filled');
    xlabel('PC1','FontSize', 16);ylabel('PC2','FontSize', 16); zlabel('PC3','FontSize', 16);
    title('3D projection', 'FontSize', 16);
    
    % settings
    cb = colorbar;
    cb.Label.String = 'Time to seizure [min]';
    cb.Label.FontSize = 14;
    cb.FontSize = 14;
    cb.Units = 'normalized';
    cb.TicksMode = 'manual';
    cb.Position = [0.93 0.1400 0.0167 0.7150];
    cb.AxisLocation = 'out';
    colormap(winter);
    pbaspect ([1,1,1]);
    set(gca,'FontSize',14)
    
    %% Extra analysis - explained variance
    % if the user want the code to calculate the explained variance
    %(Extra_analysis == 1), the code will calculate the sum of the biggest eigenvalues of matrix C
    %(in accordance to the variable 'PCA_max_dim'). then it will divide it in the
    %sum of all the eigenvalues of matrix C to get explained variance.
    
    if Extra_analysis == 1
        if n_file == 1
            explained_v = zeros (PCA_max_dim, n_measures);
        end
        
        Eigen_v = eig(C);
        n_eigen_vec = length (Eigen_v);
        Eigen_v_sum = sum(Eigen_v);
        
        for N_dim = 1: PCA_max_dim
            Eigen_v_part = sum(Eigen_v (n_eigen_vec-N_dim +1 : n_eigen_vec));
            explained_v(N_dim ,n_file) =   Eigen_v_part / Eigen_v_sum;
        end
        
        
        % plot of change in explained variance as function of PCA dimennsions.
        if n_file == n_measures
            precentage_explained_v = 100*explained_v;
            figure;
            bars = bar(precentage_explained_v);
            
            hold on
            mean_exp = mean(precentage_explained_v,2);
            plot_a = plot(1:PCA_max_dim,mean_exp);
            mean_sub_var = scatter(1:PCA_max_dim,mean_exp);
            barslgd = legend([bars, mean_sub_var,plot_a],...
                [main_data(1,:),'mean explained variance','trend'],'FontSize',12);
            barslgd.Location = 'northwest';
            title(barslgd,'patient number');
            xlabel('number of PCA dimenissions','FontSize',14)
            ylabel('explained variance [%]','FontSize',14)
            set(gca,'FontSize',14)
            title('Explained variance as function of PCA dimenssions')
        end
    end
end
