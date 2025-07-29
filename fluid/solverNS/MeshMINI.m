% input: X, T 

% output: X, T mesh for velocity and Xp, Tp, for pressure

function [mesh] = MeshMINI(mesh)
    X =mesh.X;
    Xp=X;
    Tp=mesh.T;
    
    T=zeros(size(Tp,1),4);
    T(:,1:3)=mesh.T;
    for i=1:size(T,1)
        %x1 = [X(T(i,1),1) X(T(i,2),1) X(T(i,3),1)];
        %y1 = [X(T(i,1),2) X(T(i,2),2) X(T(i,3),2)];
        %polyin = polyshape({x1},{y1});
        %[x,y] = centroid(polyin);
        %X=[X; x y];
        %T(i,4) = size(X,1);
        T(i,4) = size(X,1)+i;
    end

    mesh.T = T;
    mesh.Tp= Tp;
    mesh.Xp= Xp;
   
    %xy = ( X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:) )/3.0;
    %fileName = './temps/z_temp_mini';
    %X_mini = [X;xy];
    %DT = delaunay(X_mini);
    %mesh.X = X_mini';
    %mesh.T = DT;
    %msh = printMSH(fileName,mesh);
end