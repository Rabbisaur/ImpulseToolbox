numcond = 32;
itm = 1:32;

cndfilepath = './NCOR1.CND';
fp = fopen(cndfilepath,'w');


fprintf(fp,'COND# BCKGND TIMING FIX_ID TEST0 TEST1 TEST2 ---COLOR-PALETTE---');
fprintf(fp,'\n\r\n');

for thiscond = 1:numcond
    numzeros = (floor(log10(thiscond))+1);
fprintf(fp,[repmat(' ',1,5-numzeros),'%d     -4      1     -3    -1 ',repmat(' ',1,5-numzeros), '%d     0'],thiscond,itm(thiscond));
fprintf(fp,'\n');
end

fclose(fp);