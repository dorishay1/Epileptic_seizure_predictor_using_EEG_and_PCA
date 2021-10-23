function [datacell, n_sub, n_elect,time_bins] = data2cell(zip_file_name)
%This function takes a zipfile with mat files and converts it to a easy to handle data cell.
%the function also saves for each mesurement the original patient number and the seasiure number
%as writen in the file name.
%number of total measures, electrodes and time sampls also preduce.
unzip(zip_file_name,'..\DATA_DIR\');        %unziping the data.

%first filter: files will contain only files in 'DATA_DIR' folder and subfolders
%that end with '.mat'.
files = dir('..\DATA_DIR\**/*.mat');

%This is the allowed pattern. The first digits after 'p' that shown will be stored as
%patient number, next, the digits after 's' will be stored as seizure number.
look_str = 'p(?<patient_number>\d+)_s(?<seizure>\d+)';


%Generally, for each file in files (that contain only the mat files), this loop
%will check and extract the details from the file name. Using 'load' will read
%the file and store the relevant information.

%The genral idea of the loop -  making sure that the file name is according to the
%requasted pattern and storing the relvent info.

%preparing memory - 1st row is for patient id, 2nd is for seizure number and 3rd is for
%the data mat.
%the length is number of valid measures.

datacell = cell(3,length(files));


for i = 1:length(files)
    current_str = string(files(i).name);
    
    %checks if file name(current_str) is fit for the allowed pattern (look_str).
    %if it's not, error msg is shown. if it is, it stores it in temp variable.
    current_data = regexp(current_str,look_str,'names');
    if isempty(current_data) == 1
        error('unvalid file name: %s' ,current_str)
    end
    
    datacell{1,i} = char(current_data.patient_number);
    
    datacell{2,i} = char(current_data.seizure);
    
    %we use 'fullfile' in order to get the path for the files (they are not
    %stored in the same directory as the main code.
    current_file = fullfile(files(i).folder,files(i).name);
    current_file = load(current_file);
    
    datacell{3,i} = current_file.data;
    
    
end


%setting number of subjects, number of electrodes & time bins.
n_sub = length(datacell);

[n_elect, time_bins] = size(datacell{3,1});




end

