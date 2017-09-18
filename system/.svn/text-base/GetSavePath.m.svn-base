function [savepath savedir] = GetSavePath(filepath,filename)
poz = find(filepath =='.',1,'last');
SizeExt = length(filepath) - poz;
savedir = filepath(1:end-SizeExt-1);
savepath = [savedir '/' filename];