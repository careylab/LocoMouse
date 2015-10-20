function data = getDataList(path)
data = lsOSIndependent(path);
[path,~,~] = fileparts(path);
data = [repmat([path filesep],size(data,1),1) data];