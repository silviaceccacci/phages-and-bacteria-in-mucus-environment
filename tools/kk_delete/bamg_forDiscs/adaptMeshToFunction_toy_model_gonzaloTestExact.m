
% n is number of figure where the solution is
function [X,T,nodesDISC,elementsDISC,n,exacta,model] = adaptMeshToFunction_toy_model_gonzalo()
    
    exacta = 1; %1 geometria exacta, 0 no
    maxIt = 25;
    
    %parametres del bamg
    beta = .1;
    
    data.coef   = '1'; %default: 1
    data.ratio  = '0'; %default: 0
    data.hmin   = '1'; %0 without hmin
    data.hmax   = '20'; %0 without hmax
    data.errg   = '0.001'; %never bigger than 0.7071...
    data.err    = '0.001';
    data.iso    = '0'; % 1 si volem que sigui isomètric, 0 vol dir anisometric
    data.anisomax = '10'; %0 without anisomax
    data.nbv = '3000';
    
    
    ne1dx = 1000;
    ne1dy = 100;
    [X,T] = generateRectangleMesh(ne1dx,ne1dy);
    
    
    x_LL = 100;
    x_LR = 300;
    y_L = 125;
    
    countFig = 1;
    figure(2)
    subplot(4,2,countFig)
    plotMesh(X,T,0)
    
    %nou domini [-x_LL,x_LR]x[-y_L,y_L]
    for i=1:size(X,1)
        X(i,1)=(x_LL+x_LR)*X(i,1)+100;
        X(i,2)=y_L*X(i,2);
    end
    countFig = countFig+1;
    subplot(4,2,countFig)
    plotMesh(X,T,0)
    
    
    it=0;
    isConverged = false;
    numElemsOld = size(T,1);
    
    if exacta == 1
        
        model = 2;  %model = 1 (linia {0}x[-25,25]), 2 (disc [-1,1]x[-25,25], 3 (2 discos [-40,-38]x[-25,25] i [40,42]x[-25,25])
                    %els he fet rapid, no estan optimitzats (sobretot el 3) pero funcionen
        
        while(it<=maxIt && isConverged==false)
            it
        
            switch model
                case 1
                        U = Fun2(X(:,1),X(:,2),beta);
                case 2
                        U = Fun3(X(:,1),X(:,2),beta);
                case 3
                        U = Fun7(X(:,1),X(:,2),beta);
                otherwise
                    error('I do not know this model');
            end                    
        
        
            [X,T,U] = adapt_toy_model_gonzalo_exact(X,T,U,data,model);
        
            it = it+1;
            numElems = size(T,1);
            isConverged =  abs( numElemsOld -numElems ) < numElems*0.05; %why?
            numElemsOld = numElems;
        
            figur = floor(countFig/8);
            figure(2+figur)
            countFig = countFig+1;
            p=mod(countFig,8);
            if p==0
                p=8;
            end
            subplot(4,2,p)
            plotMesh(X,T,0)
        end
        
        %plot de la funció tanh(beta*d(x,y)) en la malla final
        PlotDistanceFinalMesh(X,T,model,beta);
    
    
        if(it>maxIt && isConverged==false)
            error('there is no convergence');
        end
        
        %plot the final disc
        figure;
        n = get(gcf,'Number');
        triplot(T, X(:,1), X(:,2))
        set(gca,'Fontsize',16); 
        colormap(jet)
        lighting flat
        axis equal
        title('Final Mesh with Exact Geometry')
        dim = [0.8 0.4 0.3 0.3];
        str = {['model = ' num2str(model)],['coef = ' data.coef],['ratio = ' data.ratio],['hmin = ' data.hmin],['hmax = ' data.hmax],['errg = ' data.errg],['erg = ' data.err],['iso = ' data.iso],['anisomax = ' data.anisomax],['nbv = ' data.nbv],['beta = ',num2str(beta)]};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        
        switch model
            case 1
                nodesDISC = FindNodesDisc(X,8);
                x1 = linspace(-25,25);
                y1 = zeros(size(x1));
                hold on, plot(y1,x1,'red','LineWidth',2)
                hold on, plot(X(nodesDISC,1),X(nodesDISC,2),'bo','MarkerSize',10,'MarkerEdgeColor','k')
                

            case 2
                nodesDISC = FindNodesDisc(X,3);
                x = [-1, 1, 1, -1, -1];
                y = [-25, -25, 25, 25, -25];
                hold on, plot(x, y, 'r','LineWidth', 3)
                hold on, plot(X(nodesDISC,1),X(nodesDISC,2),'bo','MarkerSize',10,'MarkerEdgeColor','k')
                
            case 3
                nodesDISC = FindNodesDisc(X,7);
                x = [-40, -38, -38, -40, -40];
                y = [-25, -25, 25, 25, -25];
                hold on, plot(x, y, 'r','LineWidth', 3)
                x = [40, 42, 42, 40, 40];
                y = [-25, -25, 25, 25, -25];
                hold on, plot(x, y, 'r','LineWidth', 3)
                hold on, plot(X(nodesDISC,1),X(nodesDISC,2),'bo','MarkerSize',10,'MarkerEdgeColor','k')           
                
            otherwise
                error('Not implemented yet');
        end
        
    else
        
        
        
        
        
        model = 7; %between 1 and 8
    
        %plot de la funció tanh(beta*d(x,y)) (aquesta comanda acostuma a ser cara d'evaluar, així que només quan volem veure la distancia)
        %PlotDistance(model,beta);
    
        %load('test');
    
        while(it<=maxIt && isConverged==false)
            it
        
            switch model
                case 1
                        U = Fun1(X(:,1),X(:,2),beta);
                case 2
                        U = Fun2(X(:,1),X(:,2),beta);
                case 3
                        U = Fun3(X(:,1),X(:,2),beta);
                case 4
                        U = Fun4(X(:,1),X(:,2),beta);
                case 5
                        U = Fun5(X(:,1),X(:,2),beta);
                case 6
                        U = Fun6(X(:,1),X(:,2),beta);
                case 7
                        U = Fun7(X(:,1),X(:,2),beta);
                case 8
                        U = Fun8(X(:,1),X(:,2),beta);
                otherwise
                    error('I do not know this model');
            end
        
            [X,T,U] = adapt_toy_model_gonzalo(X,T,U,data);

            it = it+1;
            numElems = size(T,1);
            isConverged =  abs( numElemsOld -numElems ) < numElems*0.05; %why?
            numElemsOld = numElems;
        
            figur = floor(countFig/8);
            figure(2+figur)
            countFig = countFig+1;
            p=mod(countFig,8);
            if p==0
                p=8;
            end
            subplot(4,2,p)
            plotMesh(X,T,0)
        end
    
    %find nodes inside the final DISC
    nodesDISC = FindNodesDisc(X,model);
    
    %plot the final disc
    figure;
    n = get(gcf,'Number');
    plotMeshFinal(X,T,model,0,nodesDISC)
    title('Final Mesh')
    dim = [0.8 0.4 0.3 0.3];
    str = {['model = ' num2str(model)],['coef = ' data.coef],['ratio = ' data.ratio],['hmin = ' data.hmin],['hmax = ' data.hmax],['errg = ' data.errg],['erg = ' data.err],['iso = ' data.iso],['anisomax = ' data.anisomax],['nbv = ' data.nbv],['beta = ',num2str(beta)]};
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
    %plot de la funció tanh(beta*d(x,y)) en la malla final
    PlotDistanceFinalMesh(X,T,model,beta);
    
    
    if(it>maxIt && isConverged==false)
        error('there is no convergence');
    end
    
    end
    
        %look for the elements inside the disc
        elementsDISC = [];
        for i=1:size(T,1)
             AreDisc = ismember(T(i,:),nodesDISC);
             %disp(T(i,:));
             %disp(nodesDISC);
             if sum(AreDisc) == 3
                  elementsDISC=[elementsDISC, i];
             end
         end
        

end