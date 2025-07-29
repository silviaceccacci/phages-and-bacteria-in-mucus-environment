%input: mesh
%output: vector with the Dirichlet BC prescribed (first all x components,
%then y components).

function val = bc_cos(X,parameters)

x = X(:,1); y = X(:,2);
val = zeros(length(x)+length(y),1);
nodesone = 1:length(x);

val(nodesone)        = 0.1 + cos(y*pi)/10;
val(length(x)+1:end) = 0.0;  
