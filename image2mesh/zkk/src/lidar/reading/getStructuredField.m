function [ field ] = getStructuredField( dem )

demProxy_i_j = proxyDEM(dem);
demProxy_j_i = demProxy_i_j';
demProxy_jInv_i = demProxy_j_i(end:-1:1,:);

demOutput = demProxy_i_j;

kkindex = find(demOutput(:)<-1000);
demOutput(kkindex) = 0.0;

field.structured = true;
field.z = demOutput;
field.hx = dem.cellSize;
field.hy = dem.cellSize;
field.x0 = [dem.x0 dem.y0];
field.nx = dem.ni*(dem.NI+1);
field.ny = dem.nj*(dem.NJ+1);
%field.width = width;
%field.height = height;
%field.x = px;
%field.y = py;
%field.points = [Xg(:) Yg(:)];

end
