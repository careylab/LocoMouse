% Track collections of videos as they come from the setup.

%% User input:  check paths
locomouse_path = 'path\to\LocoMouse';
session_path = fullfile(locomouse_path,'movies'); % path where videos for tracking are stored
output_path = fullfile(locomouse_path,'output'); % output path for tracking files
model_path = fullfile(locomouse_path,'tracking code','model_01112013.mat');  % path to model file
dist_image_path = fullfile(locomouse_path,'output','Distortion images'); % output path for images
distortion_path =fullfile(locomouse_path,'tracking code','IDX_pen.mat');  % path to calibration file
error_path = fullfile(locomouse_path,'output','error');  % output path for tracking errors

skip_analyzed_videos = false; % skip files 
skip_undistortion = true; % skip undistortion if images are available
specific = false;

%make input file structure is compatible (see example movies folder)
specific_data_list = {{'S1_1_21_2014',{'G6AK4'},{'9'}}};%{'01_11_13_S5',{'G6AE9'},'9'}


%% Browsing the data:
browsing_data_structure_dfolder
