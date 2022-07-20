% combines sac files into a saf file
%
% example:
% sac2saf('01.G014.HHZ.sac','01.G014.HH1.sac','01.G014.HH2.sac','01.G04.saf')
%
% Fatih Bulut, Bogazici University, 2020


function [] = sac2saf(input1,input2,input3,output)
x1 = sacread(input1);
x2 = sacread(input2);
x3 = sacread(input3);

SAYAC = min([ length(x1.data) length(x2.data) length(x3.data)]);


fid = fopen(output,'w');

out  = fprintf(fid,'SESAME ASCII data format (saf) v. 1   (this line must not  be modified)\n');
out  = fprintf(fid,'# Data converted from SEISAN with program SAFSEI\n');
out  = fprintf(fid,'# Original file name: 2008.SEI\n',x1.aheader(1:4));
out  = fprintf(fid,'STA_CODE = %4s\n',x1.aheader(1:4));
out  = fprintf(fid,'SAMP_FREQ =%9.3f\n',1/x1.fheader(1));
out  = fprintf(fid,'NDAT =%9.0f\n',SAYAC);

out  = fprintf(fid,'START_TIME =  %4.0f  0 00 00 00  00.000\n',x1.iheader(1));
out  = fprintf(fid,'UNITS = count\n');
out  = fprintf(fid,'CH0_ID = S Z\n');
out  = fprintf(fid,'CH1_ID = S N\n');
out  = fprintf(fid,'CH2_ID = S E\n');
out  = fprintf(fid,'####-------------------------------------------\n');


for i = 1:SAYAC
 out  = fprintf(fid,'%12.3f %12.3f %12.3f\n',x1.data(i),x2.data(i),x3.data(i));
end

fclose(fid);


function sac = sacread(filename)
sac.filename = filename;
sacfile = fopen(sac.filename,'rb');
sac.fheader = fread(sacfile,70,'float32');
sac.iheader = fread(sacfile,40,'int32');
sac.aheader = (  char(  fread(sacfile,192,'char')  )  );
sac.data = fread(sacfile,'float32');
fclose(sacfile);
