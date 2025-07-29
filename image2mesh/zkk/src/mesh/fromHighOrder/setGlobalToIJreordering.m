function [element] = setGlobalToIJreordering(element)


    switch(element.type)

        case 'tri'
            [element] = setGlobalToIJreorderingTri(element);

        case 'quad'
            [element] = setGlobalToIJreorderingQuad(element);
            
        case 'tet'
            [element] = setGlobalToIJreorderingTet(element);
    end


end

function [element] = setGlobalToIJreorderingTri(element)

    p=element.order;
    
    IJToGlobal = zeros(p+1,p+1);
    globalToIJ = zeros(element.numNod,2);

    IJToGlobal(1,1) = 1 ;
    IJToGlobal(1,p+1) = 2 ;
    IJToGlobal(p+1,1) = 3 ;
    globalToIJ(1,:) = [1 1];
    globalToIJ(2,:) = [1 p+1];
    globalToIJ(3,:) = [p+1 1];

    
    for i=2:p
        IJToGlobal(1,i)     = 2+i;
        IJToGlobal(i,p+2-i) = 2+1*(p-1)+i;
        IJToGlobal(p+2-i,1) = 2+2*(p-1)+i;

        globalToIJ(2+i,:)         = [1 i];
        globalToIJ(2+(p-1)+i,:)   = [p+2-i  i];
        globalToIJ(2+2*(p-1)+i,:) = [p+2-i  1];
    end

    counter = element.numNod;% esta invertit respecte ez4u
    for i=2:p
        for j=2:(p+1-i)
            IJToGlobal(i,j)       = counter;
            globalToIJ(counter,:) = [i j];
            counter = counter-1;
        end
    end
        
    element.IJToGlobal=IJToGlobal;
    element.globalToIJ=globalToIJ;
    
end


function [element] = setGlobalToIJreorderingQuad(element)

    p=element.order;
    
    IJToGlobal = zeros(element.order+1,element.order+1);
    globalToIJ = zeros(element.numNod,2);

    IJToGlobal(1,1) = 1 ;
    IJToGlobal(1,p+1) = 2 ;
    IJToGlobal(p+1,p+1) = 3 ;
    IJToGlobal(p+1,1) = 4 ;
    globalToIJ(1,:) = [1 1];
    globalToIJ(2,:) = [1 p+1];
    globalToIJ(3,:) = [p+1 p+1];
    globalToIJ(4,:) = [p+1 1];

    
    for i=2:p
        IJToGlobal(1,i)     = 3+i;
        IJToGlobal(i,p+1)   = 3+1*(p-1)+i;
        IJToGlobal(p+1,p+2-i)= 3+2*(p-1)+i;
        IJToGlobal(p+2-i,1) = 3+3*(p-1)+i;

        globalToIJ(3+i,:)         = [1 i];
        globalToIJ(3+(p-1)+i,:)   = [i  p+1];
        globalToIJ(3+2*(p-1)+i,:)   = [p+1  p+2-i];
        globalToIJ(3+3*(p-1)+i,:) = [p+2-i  1];
    end
     
    counter = 4*p;
    for i=2:p
        for j=2:p
            counter = counter +1;
            IJToGlobal(i,j)       = counter;
            globalToIJ(counter,:) = [i j];
        end
    end
    
    element.IJToGlobal=IJToGlobal;
    element.globalToIJ=globalToIJ;

end


function [element] = setGlobalToIJreorderingTet(element)

    p=element.order;
    
    IJToGlobal = zeros(element.order+1,element.order+1,element.order+1);
    globalToIJ = zeros(element.numNod,3);

    %only nodes under x+y+z plane are taken
    count = 0;
    for k=1:p
        for j=1:p
            for i=1:p
                if( i+j+k<= (p+1) )
                    count = count+1;
                    IJToGlobal(i,j,k) = count;
                    globalToIJ(count,:) = [i j k];
                end
            end
        end
    end
    
    element.IJToGlobal=IJToGlobal;
    element.globalToIJ=globalToIJ;
    
end




