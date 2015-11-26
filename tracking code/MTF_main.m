% Track collections of videos as they come from the setup.
%
% The naming convention of the system is:
% SessionID/MouseID/MouseID_md1_..._mdn_type_SxTy.avi
%
% where: 
% SessionID: usually sesson number and date.
% MouseID: a code that identifies the animal.
% mdx: are info about the mouse (size, weight, etc.);
% type: control/pcd.
% Sx: x is session number.
% Ty: y is trial number. 

%% User input:  check paths
locomouse_path = 'path\to\LocoMouse'; % path for the LocoMouse root folder

session_path = fullfile(locomouse_path,'movies'); % path where videos for tracking are stored
output_path = fullfile(locomouse_path,'output'); % output path for tracking files
model_path = fullfile(locomouse_path,'tracking code','model_01112013.mat');  % path to model file
dist_image_path = fullfile(locomouse_path,'output','Distortion images'); % output path for images
distortion_path =fullfile(locomouse_path,'tracking code','IDX_pen.mat');  % path to calibration file
error_path = fullfile(locomouse_path,'output','error');  % output path for tracking errors

skip_analyzed_videos = false; % skips if data is found for that video. 
skip_undistortion = true; % skip undistortion if images are available.

%% Browsing the data:
browsing_data_structure_dfolder(session_path,output_path,model_path,dist_image_path,distortion_path,error_path,skip_analyzed_videos,skip_undistortion);
