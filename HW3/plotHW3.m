fileID = fopen('filtered0.01.txt','r');
%Define the format of the data to read. Use '%f' to specify floating-point numbers.

formatSpec = '%f';
%Read the file data, filling output array, A, in column order. fscanf reapplies the format, formatSpec, throughout the file.

XM2 = fscanf(fileID,formatSpec);

fclose(fileID);

fileID = fopen('filtered1.txt','r');
XM0 = fscanf(fileID,formatSpec);
fclose(fileID);

fileID = fopen('filtered100.txt','r');
XP2 = fscanf(fileID,formatSpec);
fclose(fileID);

fileID = fopen('filtered10000.txt','r');
XP4 = fscanf(fileID,formatSpec);
fclose(fileID);

fileID = fopen('denoise.txt','r');
XCOR = fscanf(fileID,formatSpec);
fclose(fileID);

figure(1)
hold
plot([1:1000],XCOR)
plot([1:1000],XM2)
plot([1:1000],XM0)
plot([1:1000],XP2)
plot([1:1000],XP4)
title('filtered')
legend('XCOR','0.01','1','100','10000')

figure(2)
hold
plot([1:1000],XCOR-XM0)
plot([1:1000],XCOR-XM2)
plot([1:1000],XCOR-XP2)
plot([1:1000],XCOR-XP4)
title('noise')
legend('0.01','1','100','10000')

