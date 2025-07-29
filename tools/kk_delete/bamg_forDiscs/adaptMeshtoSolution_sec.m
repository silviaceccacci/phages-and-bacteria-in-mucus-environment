
% n is number of figure where the solution is
function [X,T,nodesDISC,elementsDISC,n] = adaptMeshtoSolution(...
    ifplot,exacta,data,X,T,Xp,Tp,adapt_variable,ifpreviousmesh,omega,mesh_solution)
    
    if exacta == 1
        error('Still to be done. I do not know if it can be done');
    else
        
    %PARAMETERS  
    maxIt = 1; %maximum interations in Bamg (DO NOT CHANGE TO AVOID PRESSURE INTERPOLATION PROBLEMS)
    n = 0;
    
    %SETING THE INITIAL MESH
    if ifpreviousmesh == 1 %from the previous mesh
%         countFig = 1;
%         if ifplot == 1
%             figure;
%             n = get(gcf,'Number');
%             subplot(4,2,countFig)
%             %plotMesh(Xp,Tp,0)
%         end
%     
%         countFig = countFig+1;
%         if ifplot == 1
%             subplot(4,2,countFig)
%             %plotMesh(Xp,Tp,0)
%         end
    else %from zero
        error('should not enter here...')
%         [X,T] = generateRectangleMesh(data.ne1dx,data.ne1dy);
%     
%         countFig = 1;
%         if ifplot == 1
%             figure;
%             n = get(gcf,'Number');
%             subplot(4,2,countFig)
%             %plotMesh(X,T,0)
%         end
%     
%         %nou domini [-x_LL,x_LR]x[-y_L,y_L]
%         for i=1:size(X,1)
%             X(i,1)=((-omega.x_LL+omega.x_LR)/2)*X(i,1)+((-omega.x_LL+omega.x_LR)/2)+omega.x_LL;
%             X(i,2)=((-omega.y_LD+omega.y_LU)/2)*X(i,2)+((-omega.y_LD+omega.y_LU)/2)+omega.y_LD;
%         end
%         countFig = countFig+1;
%         if ifplot == 1
%             subplot(4,2,countFig)
%             %plotMesh(X,T,0)
%         end
%         
%         %INTERPOLATION IN THE NEW MESH
%         %save the mesh
%         fileName = './temps/z_temp2';%mesh goes here
%         mesh.X = X';
%         mesh.T = T;
%         msh = printMSH(fileName,mesh);
%         
%         fileName = './temps/z_temp_sol';%interpolation goes here
%         %in order to do the interpolation, we need the Bamg
%         [bamg_exe]=get_bamg_exe();
%         %interpolation in the mesh
%         adaptativeStep = [ bamg_exe ' -b ./temps/z_temp.msh' ' -r ./temps/z_temp2.msh '  ' -rbb ./temps/z_temp_sol.bb ' ' -wbb ./temps/z_temp_sol.bb' ];
%         eval(adaptativeStep)
%         interpolated = readBAMGbb(fileName);
%         clear interpolated;
%         
%         %change it to the mesh
%         fileName = './temps/z_temp';%mesh goes here
%         msh = printMSH(fileName,mesh);

    end
    
    
    %print the solution in the file
    fileName_metric_input      = './temps/temp_metric_input';
    fileName_solution_input    = './temps/temp_solution_input';
    fileName_pressure_input    = './temps/temp_pressure_input';
    fileName_combination_input = './temps/temp_combination_input';
    
    it=0;
    isConverged = false;
    numElemsOld = size(T,1);
    while(it < maxIt && isConverged==false)
        %it
        
        if adapt_variable.comb == 1
            metric_input = ' '; solution_input = ' '; pressure_input = ' ';
            
            if adapt_variable.d == true && adapt_variable.u == true && adapt_variable.p == true
                U = Fun(X(:,1),X(:,2),data.beta,omega);
                mod_u = mesh_solution.u;
                p = mesh_solution.p;
                U2 = adapt_variable.alfa*U+adapt_variable.beta*mod_u'+adapt_variable.gamma*p';
                combination_input = printSOL(fileName_combination_input,U2');
                
            elseif adapt_variable.d == true && adapt_variable.u == true && adapt_variable.p == false
                U = Fun(X(:,1),X(:,2),data.beta,omega);
                mod_u = mesh_solution.u;
                U2 = adapt_variable.alfa*U+adapt_variable.beta*mod_u';
                combination_input = printSOL(fileName_combination_input,U2');
                              
            elseif adapt_variable.d == true && adapt_variable.u == false && adapt_variable.p == true
                U = Fun(X(:,1),X(:,2),data.beta,omega);
                p = mesh_solution.p;
                U2 = adapt_variable.alfa*U+adapt_variable.gamma*p';
                combination_input = printSOL(fileName_combination_input,U2');
                
            elseif adapt_variable.d == true && adapt_variable.u == false && adapt_variable.p == false
                U = Fun(X(:,1),X(:,2),data.beta,omega);
                U2 = adapt_variable.alfa*U;
                combination_input = printSOL(fileName_combination_input,U2');
                
            elseif adapt_variable.d == false && adapt_variable.u == true && adapt_variable.p == true
                mod_u = mesh_solution.u;
                p = mesh_solution.p;
                U2 = adapt_variable.beta*mod_u'+adapt_variable.gamma*p';
                combination_input = printSOL(fileName_combination_input,U2');
                
            elseif adapt_variable.d == false && adapt_variable.u == true && adapt_variable.p == false
                mod_u = mesh_solution.u;
                U2 = adapt_variable.beta*mod_u';
                combination_input = printSOL(fileName_combination_input,U2');
                
            elseif adapt_variable.d == false && adapt_variable.u == false && adapt_variable.p == true
                p = mesh_solution.p;
                U2 = adapt_variable.gamma*p';
                combination_input = printSOL(fileName_combination_input,U2');
                
            else
                error('wrong combination');
            end
            
        else
            combination_input = ' ';
            solution_input = printSOL(fileName_solution_input,mesh_solution.u');
            pressure_input = printSOL(fileName_pressure_input,mesh_solution.p');
            U2 = Fun(X(:,1),X(:,2),data.beta,omega);
            metric_input = printSOL(fileName_metric_input,U2');
        end
        
        %save in files
        [X,T,u_inter,p_inter] = adapt_nonexact_after_adaptation(Xp,Tp,data,...
            adapt_variable,metric_input,solution_input, pressure_input,combination_input);
        
        %delete files
        delete(metric_input);
        delete(solution_input);
        delete(pressure_input);
        
        it = it+1;
        numElems = size(T,1);
        isConverged =  abs( numElemsOld -numElems ) < numElems*0.05;
        numElemsOld = numElems;

    end
    
    
    %if(it>maxIt && isConverged==false)
    %    error('there is no convergence');
    %end
    
    nodesDISC    = FindNodesDisc(X,omega);
    [elementsDISC,nodesDISC] = computeElementsDisc(X,T,omega,nodesDISC);
    
    fprintf('  numElem: %d   numNodes: %d \n',size(T,1),size(X,1))
    
    %plot the final disc
    n = 0; %to return something
    if ifplot == 1
        figure;
        n = get(gcf,'Number');
        plotMeshFinal(X,T,omega,0,nodesDISC(:,1))
        title('Final Adapted Mesh')
        dim = [0.8 0.4 0.3 0.3];
        str = {['model = ' num2str(omega.model)],['coef = ' data.coef],['ratio = ' data.ratio],['hmin = ' data.hmin],['hmax = ' data.hmax],['errg = ' data.errg],['erg = ' data.err],['iso = ' data.iso],['anisomax = ' data.anisomax],['nbv = ' data.nbv],['beta = ',num2str(data.beta)]};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
    end
    
    
    %plot de la funció beta*tanh(beta*d(x,y))+solucio en la malla final
    % U = Fun(X(:,1),X(:,2),data.beta,omega);
    if ifplot == 1    
        if adapt_variable.d == true && adapt_variable.u == true && adapt_variable.p == true
            U = Fun(X(:,1),X(:,2),data.beta,omega);
            U2 = adapt_variable.alfa*U+adapt_variable.beta*u_inter'+adapt_variable.gamma*p_inter';
            PlotDistanceFinalMesh2(X,T,omega,data.beta,U2,adapt_variable);

        elseif adapt_variable.d == true && adapt_variable.u == true && adapt_variable.p == false
            U = Fun(X(:,1),X(:,2),data.beta,omega);
            U2 = adapt_variable.alfa*U+adapt_variable.beta*u_inter';
            PlotDistanceFinalMesh2(X,T,omega,data.beta,U2,adapt_variable);

        elseif adapt_variable.d == true && adapt_variable.u == false && adapt_variable.p == true
            U = Fun(X(:,1),X(:,2),data.beta,omega);
            U2 = adapt_variable.alfa*U+adapt_variable.gamma*p_inter';
            PlotDistanceFinalMesh2(X,T,omega,data.beta,U2,adapt_variable);

        elseif adapt_variable.d == true && adapt_variable.u == false && adapt_variable.p == false
            U = Fun(X(:,1),X(:,2),data.beta,omega);
            U2 = adapt_variable.alfa*U;
            PlotDistanceFinalMesh2(X,T,omega,data.beta,U2,adapt_variable);

        elseif adapt_variable.d == false && adapt_variable.u == true && adapt_variable.p == true
            U2 = adapt_variable.beta*u_inter'+adapt_variable.gamma*p_inter';
            PlotDistanceFinalMesh2(X,T,omega,data.beta,U2,adapt_variable);

        elseif adapt_variable.d == false && adapt_variable.u == true && adapt_variable.p == false
            U2 = adapt_variable.beta*u_inter';
            PlotDistanceFinalMesh2(X,T,omega,data.beta,U2,adapt_variable);

        elseif adapt_variable.d == false && adapt_variable.u == false && adapt_variable.p == true
            U2 = adapt_variable.gamma*p_inter';
            PlotDistanceFinalMesh2(X,T,omega,data.beta,U2,adapt_variable);

        else
            error('wrong combination');
        end
          
    end
    
    
    end
        

end