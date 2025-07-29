function [f_out]=interpolatePoints(F_interp,x_out)

nfields = size(F_interp([0,0]),2);
nout = size(x_out,1);

f_out = zeros(nout,nfields);
for i=1:nout
    f_out(i,:) = F_interp(x_out(i,:));
end
