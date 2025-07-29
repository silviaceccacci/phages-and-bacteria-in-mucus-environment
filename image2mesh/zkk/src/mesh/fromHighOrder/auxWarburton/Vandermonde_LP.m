function V = Vandermonde_LP(nDeg,coord,element)
%
% V = Vandermonde_LP(nDeg,coord)
%
% Function to compute the Vandermonde matrix: V(i,j) = p_j(xi_i)
%
% Input:
% nDeg:  degree of polynomials
% coord: nodal coordinates coord(i,:) = [xi,eta] (or [xi,eta,zeta]) of node i
%
% Output:
% V:     Vandemonde matrix
%

    nsd = size(coord,2);
    if nsd==1
        V = Vandermonde_LP1D(nDeg,coord);
    elseif nsd==2
        V = Vandermonde_LP2D(nDeg,coord,element);
    elseif nsd ==3
        V = Vandermonde_LP3D(nDeg,coord,element);
    else
        error('Vandermonde_LP requires coordinates in 1D, 2D or 3D')
    end
end

%
% 1D matrix
%
function V = Vandermonde_LP1D(nDeg,coord)

    N = length(coord);%number of nodes/polynomials
    if ( N ~= nDeg+1 ) 
        error('The number of polynomials does not coincide with the number of nodes')
    end

    % V  = zeros(N,N);
    % 
    % for i = 1:N
    %  x = coord(i);
    %  p = orthopoly1D(x,nDeg);
    %  V(i,:) = p';
    % end
    p = orthopoly1D(coord,nDeg);
    V = p';
end

%
% 2D matrix
%
function V = Vandermonde_LP2D(nDeg,coord,element)
%     disp('Cuidao que a Vandermonde_LP2D he desactivat el checking de trisHO')
%     N = size(coord,1);%number of nodes/polynomials
%     if ( N ~= (nDeg+1)*(nDeg+2)/2 ) 
%         error('The number of polynomials does not coincide with the number of nodes')
%     end

%     p = orthopoly2D(coord,nDeg,element);
    p = orthopoly2D_general(coord,nDeg,element);
    V = p';
end

%
% 3D matrix
%
function V = Vandermonde_LP3D(nDeg,coord,element)

%     N = size(coord,1);%number of nodes/polynomials
%     if ( N ~= (nDeg+1)*(nDeg+2)*(nDeg+3)/6 ) 
%         error('The number of polynomials does not coincide with the number of nodes')
%     end


    p = orthopoly3D_general(coord,element);
    V = p';
end

