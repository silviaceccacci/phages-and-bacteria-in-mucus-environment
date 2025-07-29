function [  ] = showDEMCell(demCell)

elevationProxy_i_jInv = demCell.elevationProxy_i_j(:,end:-1:1);

elevationProxy_jInv_i = elevationProxy_i_jInv';

image(elevationProxy_jInv_i);

end

