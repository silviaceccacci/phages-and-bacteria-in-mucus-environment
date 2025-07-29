function [ demCell] = readASC( file, numDesiredPixels )

% ncols 1024
% nrows 1024
% xllcorner 423236.49798656208441
% yllcorner 4574236.5009865621105
% cellsize 2.2016611099243164062
% NODATA_value -9999

ni = fscanf(file, 'ncols %d\n\n',1);
nj = fscanf(file, 'nrows %d\n\n',1);
x0 = fscanf(file, 'xllcorner %e\n\n',1);
y0 = fscanf(file, 'yllcorner %e\n\n',1);
cellSize = fscanf(file, 'cellsize %e\n\n',1);
noDataValue = fscanf(file, 'NODATA_value %e\n\n',1);

%elevation_i_jInv = fscanf(file, '%e', [nj, ni]);
elevation_i_jInv = fscanf(file, '%e', [ni, nj]);

elevation_i_j = elevation_i_jInv(:, end:-1:1);

%cellSize = sqrt(cellSize);
if(numDesiredPixels==0 || numDesiredPixels==ni)
    elevationProxy_i_j = elevation_i_j;
else
    dividePix = ni/numDesiredPixels;
    %di = int32(ni / double(dividePix));
    %dj = int32(nj / double(dividePix));
    di = dividePix;
    dj = dividePix;
    elevationProxy_i_j = elevation_i_j(1:di:end, 1:dj:end);
    if(di~=dj)
        error('cannot do it')
    end
    ni = size(elevationProxy_i_j,1);
    nj = size(elevationProxy_i_j,2);
    cellSize = cellSize*double(di);
end

demCell ={};
% dem.elevation_x_y = elevation_x_y;
demCell.elevationProxy_i_j = elevationProxy_i_j;
demCell.ni = ni;
demCell.nj = nj;
demCell.x0 = x0;
demCell.y0 = y0;
demCell.cellSize = cellSize;
demCell.noDataValue = noDataValue;



end

