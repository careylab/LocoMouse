function file_list = lsOSIndependent(string)
% LSOSINDEPENDENT OS independent file_list consistent with MATLAB ls on
% windows.
%
% Input:
%
% string: a string with the regular expression to list (e.g. '*.png' to
% list all png files.).
%
% Ouput:
% file_list: a NxM matrix where each row is the name of a file matching the
% search string.

if ispc
    file_list = ls(string);
else
    file_list = dir(string);
    file_list = cell2mat({file_list(:).name}');
end
    
