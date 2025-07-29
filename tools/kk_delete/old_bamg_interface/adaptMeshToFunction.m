function [] = adaptMeshToFunction()

    maxIt = 10;

    load('test');   
    
    countFig = 1;
    figure(2)
    subplot(2,4,countFig)
    plotMesh(X,T,0)
    
    it=0;
    isConverged = false;
    numElemsOld = size(T,1);
    while(it>maxIt || isConverged==false)
        it
        
        U = myFun(X(:,1),X(:,2));
        
        [X,T,U] = adapt(X,T,U);

        it = it+1;
        numElems = size(T,1)
        isConverged =  abs( numElemsOld -numElems ) < numElems*0.05;
        numElemsOld = numElems;

        figure(2)
        countFig = countFig+1;
        subplot(2,4,countFig)
        plotMesh(X,T,0)
    end
    
    figure;
    plotMesh(X,T,0)

end

function [u]=myFun(x,y)
    r = radius(x,y);
    z = r - (0.5);
	u=tanh( 30*z );
end

function [r]=radius(x,y)
    r=sqrt(x.^2 + y.^2);
end




