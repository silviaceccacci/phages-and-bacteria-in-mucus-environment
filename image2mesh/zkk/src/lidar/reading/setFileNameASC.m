function [ fileName ] = setFileNameASC( filePrefix, I, J)

fileName = [filePrefix,'_',int2str(I),'_',int2str(J),'.asc'];

end

