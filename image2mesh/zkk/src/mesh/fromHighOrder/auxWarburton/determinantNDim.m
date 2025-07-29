function [det] = determinantNDim(A)

    dimVec = size(A);
    numDim = length(dimVec);

    if(      numDim == 2)
        if(          dimVec(numDim)==2)
            det =  A(1,1).*A(2,2)-A(2,1).*A(1,2)  ;
        elseif(      dimVec(numDim)==3)
            det =  A(1,1).*(A(2,2).*A(3,3) - A(2,3).*A(3,2) ) ...
              - A(2,1).*(A(1,2).*A(3,3) - A(3,2).*A(1,3) ) ...
              +A(3,1).*(A(1,2).*A(2,3) - A(2,2).*A(1,3) ) ;
        end
    elseif( numDim ==3)
        if(           dimVec(numDim)==2)
            det =  A(:,1,1).*A(:,2,2)-A(:,2,1).*A(:,1,2)  ;
        elseif(      dimVec(numDim)==3)
            det =  A(:,1,1).*(A(:,2,2).*A(:,3,3) - A(:,2,3).*A(:,3,2) ) ...
              -A(:,2,1).*(A(:,1,2).*A(:,3,3) - A(:,3,2).*A(:,1,3) ) ...
              +A(:,3,1).*(A(:,1,2).*A(:,2,3) - A(:,2,2).*A(:,1,3) ) ;
        end
    elseif( numDim ==4)
        if(           dimVec(numDim)==2)
            det =  A(:,:,1,1).*A(:,:,2,2)-A(:,:,2,1).*A(:,:,1,2)  ;
        elseif(      dimVec(numDim)==3)
            det =  A(:,:,1,1).*(A(:,:,2,2).*A(:,:,3,3) - A(:,:,2,3).*A(:,:,3,2) ) ...
              - A(:,:,2,1).*(A(:,:,1,2).*A(:,:,3,3) - A(:,:,3,2).*A(:,:,1,3) ) ...
              +A(:,:,3,1).*(A(:,:,1,2).*A(:,:,2,3) - A(:,:,2,2).*A(:,:,1,3) ) ;
        elseif(      dimVec(numDim)==4)
            det = A(:,:,1,1).*determinantNDim(A(:,:,2:4,2:4))...
                 -A(:,:,2,1).*determinantNDim(A(:,:,[1 3 4],2:4))...
                 +A(:,:,3,1).*determinantNDim(A(:,:,[1 2 4],2:4))...
                 -A(:,:,4,1).*determinantNDim(A(:,:,[1 2 3],2:4));
        end
	elseif( numDim ==5)
        if(           dimVec(numDim)==2)
            det =  A(:,:,:,1,1).*A(:,:,:,2,2)-A(:,:,:,2,1).*A(:,:,:,1,2)  ;
        elseif(      dimVec(numDim)==3)
            det =  A(:,:,:,1,1).*(A(:,:,:,2,2).*A(:,:,:,3,3) - A(:,:,:,2,3).*A(:,:,:,3,2) ) ...
              - A(:,:,:,2,1).*(A(:,:,:,1,2).*A(:,:,:,3,3) - A(:,:,:,3,2).*A(:,:,:,1,3) ) ...
              +A(:,:,:,3,1).*(A(:,:,:,1,2).*A(:,:,:,2,3) - A(:,:,:,2,2).*A(:,:,:,1,3) ) ;
        end
    else
        error('not programmed case')
    end

end