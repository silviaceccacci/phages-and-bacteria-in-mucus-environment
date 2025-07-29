function [u] = readBAMGbb(fileName)

if(strcmp(fileName(end-2:end),'.bb'))
    fileNameDat = fileName;
else
    fileNameDat = [fileName '.bb'];
end

fid = fopen(fileNameDat, 'r');

dim = [fscanf(fid, '%f %f %f %f', [4 1])];
%kk = fgetl(fid);
%kk = fgetl(fid);
%kk = fgetl(fid);
%dim = fscanf(fid,'%d',1)
%kk = fgetl(fid);

u = [fscanf(fid, '\n %f', [1 dim(3)])];

fclose(fid);