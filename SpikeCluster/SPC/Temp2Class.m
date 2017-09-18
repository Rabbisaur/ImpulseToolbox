function [Class] = Temp2Class(Class, Temp, MinSize)
% Usage: [Class] = Temp2Class(Class, Temp, MinSize)
Class = Class(Temp,3:end)+1;
Class = Class';
ClassPool = unique(Class);
for i = 1:numel(ClassPool)
    if sum(Class==ClassPool(i))<MinSize
        Class(Class==ClassPool(i)) = 0;
    end
end