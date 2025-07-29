function [faceNodes]=getFace(Telem,iface)

%Telem = [ 1 2 3 4]
%Faces pointing outside
%Face1 = [ 2 3 4 ];
%Face2 = [ 1 4 3 ];
%Face3 = [ 1 2 4 ];
%Face4 = [ 1 3 2 ];

faceNodes = [ 2 3 4 
              1 4 3
              1 2 4
              1 3 2 ];

faceNodes = Telem(:,faceNodes(iface,:));

end