function [F_interp]=set_interpolant(mesh_in,f_in,adim,domain)

nfields = size(f_in,2);

isStructured = false;
if(isfield(mesh_in,"structured"))
    isStructured = mesh_in.structured;
end

L   = domain.Lchar;
x0  = domain.origin;
Lx  = domain.Lx;
Ly  = domain.Ly;

if(isStructured)
    error('implement structured interpolant')
    % 1. need meshgrid structure
    [x, y] = meshgrid(0:0.1:1, 0:0.1:1);  
    F =@(point)( interp2(x, y, f_in_grid, point(1), point(2), 'linear') );
else
    if(    nfields==1)

        F = scatteredInterpolant(mesh_in.X(:,1), mesh_in.X(:,2), f_in, 'linear', 'nearest');
        F_interp = @(x)( F( y1(x,adim,L,x0,Lx),y2(x,adim,L,x0,Ly) ) );

    elseif(nfields==2)
        F1 = scatteredInterpolant(mesh_in.X(:,1), mesh_in.X(:,2), f_in(:,1), 'linear', 'nearest');
        F2 = scatteredInterpolant(mesh_in.X(:,1), mesh_in.X(:,2), f_in(:,2), 'linear', 'nearest');
        F1_interp = @(x)( F1(y1(x,adim,L,x0,Lx),y2(x,adim,L,x0,Ly)) );
        F2_interp = @(x)( F2(y1(x,adim,L,x0,Lx),y2(x,adim,L,x0,Ly)) );
        F_interp = @(x)([F1_interp(x),F2_interp(x)]);
    else
        error('not implemented for more than 2 fields')
    end
end

end

function y=point_adim(x,adim,Lchar,x0)
    y = x;
    if(adim)
        y(1) = y(1)-x0(1);
        y(2) = y(2)-x0(2);
        y(1) = y(1)/Lchar;
        y(2) = y(2)/Lchar;
    end
end
function [y]=y1(x,adim,Lchar,x0,Lmod)
    y=point_adim(x,adim,Lchar,x0);
    y = y(1);
    if(adim)
        y = mod(y,Lmod);
    end
end
function [y]=y2(x,adim,Lchar,x0,Lmod)
    y=point_adim(x,adim,Lchar,x0);
    y = y(2);
    if(adim)
        y = mod(y,Lmod);
    end
end