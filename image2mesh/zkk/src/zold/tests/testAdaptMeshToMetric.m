function [] = testAdaptMeshToMetric(testCase)

    addpath('../.')

    switch testCase
        case 1
            test1();
        case 2
            test2();
        case 3
            test3();
        case 4
            test4();
    end

end


function [] = test1()

    load('test0');

    [Xh,Th] = adapt(X,T,h);

    figure(1)
    plot(Xh(:,1),Xh(:,2),'*')
    hold on
    plotMesh(X,T,0)

end

function [] = test2()

    maxIt = 100;

    X = [0 0; 1 0; 1 1; 0 1];
    T = [1 2 4; 2 3 4];

    it=0;
    converged = false;
    numElemsOld = size(T,1);
    while(it>maxIt || converged==false)
    %while(it<maxIt)
        it
        
        h = hTest2(X);
        [X,T] = adapt(X,T,h);

        it = it+1;
        numElems = size(T,1)
        converged = (numElemsOld + numElemsOld*0.05) > numElems;
        numElemsOld = numElems;
        
        figure(1)
        clf
        plot(X(:,1),X(:,2),'r*')
        hold on
        plotMesh(X,T,0)
    end
    
    it

end

function [h] = hTest2(x)
minH = 0.000001;
    %h = 0.0001 + abs(x(:,1)).^2 + abs(x(:,2)).^2 ;
    %h = [ 1./(0.0001+abs(x(:,1)).^4) ,  (0.*x(:,2))  , (10.*ones(size((x(:,2)))))]  ;
    h = [ 1./(minH+abs(tanh((x(:,1)-0.5)/10))).^2 ,  (0.*x(:,2))  , (10.*ones(size((x(:,2)))))]  ;
    %h = [ 1./(minH+abs(((x(:,1)-0.5)/1).^2)) ,  (0.*x(:,2))  , (10.*ones(size((x(:,2)))))]  ;
end

function [] = test3()

    load('test3');
    
    [Xh,Th] = adapt(X,T,h);

    figure(1)
    plot(X(:,1),X(:,2),'r*')
    hold on
    plotMesh(Xh,Th,0)
end

function [] = test4()

    load('test4');
    
    [Xh,Th] = adapt(X,T,h);

    figure(1)
    plot(X(:,1),X(:,2),'r*')
    hold on
    plotMesh(Xh,Th,0)
end





