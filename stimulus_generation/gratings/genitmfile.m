dirpath= './ctx';
itmpath = './NCOR1.ITM';
ctxpath = GetFilepath(dirpath,'*.ctx');

if ispc
    slash = '\';
else
    slash = '/';
end

fp = fopen(itmpath,'w');
% write the first lines
fprintf(fp,'ITEM TYPE FILLED CENTERX CENTERY BITPAN WIN_WIDE WIN_TALL HEIGHT WIDTH ANGLE OUTER INNER -R- -G- -B- C A ------FILENAME------ ---PATTERN_FILE----- FLOAT1 FLOAT2 FLOAT3 FLOAT4 INT1 INT2 INT3 INT4 CHAR1 CHAR2');
fprintf(fp,'\n');
fprintf(fp,'  -4    1      1    0.00    0.00      0                   1.00   1.00  0.00              126 126 126 x');
fprintf(fp,'\n');
fprintf(fp,'  -3    2      0    0.00    0.00      0                                      0.04  0.04  126 126 126 x');
fprintf(fp,'\n');
fprintf(fp,'  -2    2      0    0.00    0.00      0                                      0     0       0   0   0 x');
fprintf(fp,'\n');
fprintf(fp,'  -1    2      1    0.00    0.00      0                                      0.2   0     220   0   0 x');
fprintf(fp,'\n');
fprintf(fp,'   0    1      1    0.00    0.00      0                   1.00   1.00  0.00              126 126 126 x');
fprintf(fp,'\n');

% write the rest of the file
% write the rest of the file
for i = 1:numel(ctxpath)
    name = ctxpath(i).name;
    idx1 = find(name == slash,1,'last');
%     idx2 = find(name == '.', 1,'last');
    name = name(idx1+1:end);
    numzeros = 4 - (floor(log10(i))+1);
    
    fprintf(fp,[repmat(' ',1,numzeros),'%d    8           0.00    0.00      0                                                    0   0   0 n   %s'],i,name);
    fprintf(fp,'\n');
end
fclose(fp);