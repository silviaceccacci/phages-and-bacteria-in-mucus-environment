function [polylines_out,doCollapse,doCloseRegions,countLoops] =...
            collapseBoundaryOfRegions_buildings(mesh,tolCollapse,tolClose,doCollapse,doCloseRegions,...
            countLoops,maxLoops,doCheckArea,minElementAreaAllowed)


polylines_out.X = [];
polylines_out.T = [];

numRegions = size(mesh.fieldElements,1);

%% create polylines that define the regions
region_polylines = cell(numRegions,1);
polylines = zeros(0,2);
polylines_to_region = zeros(0,1);

for iregion = (mesh.groundRegion+1):numRegions
    regionElements = mesh.fieldElements{iregion};

    adjElem = mesh.matrixAdjacentElement(regionElements,:);
    linEdgeIds = find(mesh.elementField(adjElem)~=iregion);
    [bel,bed] = ind2sub(size(adjElem),linEdgeIds);

    edgeNodes = [1 2; 2 3; 3 1];
    bEdgeNodes = zeros(length(bel),2);
    for iaux = 1:length(bel)
        bEdgeNodes((iaux),1) = mesh.T(regionElements(bel(iaux)),edgeNodes(bed(iaux),1));    
        bEdgeNodes((iaux),2) = mesh.T(regionElements(bel(iaux)),edgeNodes(bed(iaux),2)); 
    end

    numBEdges = size(bEdgeNodes,1);
    
    sortedEdges = sort(bEdgeNodes,2);
    
    region_polylines{iregion} = size(polylines,1) + (1:numBEdges);%sortedEdges;%bEdgeNodes
    
    polylines = [polylines;  sortedEdges];
    polylines_to_region = [polylines_to_region; iregion*ones(numBEdges,1)];
    
end
[polylines_unique,polylines_mantained,polylines_map_to_unique] = unique(polylines,'rows');
%polylines_unique = polylines(polylines_mantained);
numPolylines_u = size(polylines_unique,1);

%% split the polylines into patches that are shared by one or two regions
numCurves = 0;
curves = cell(0,1);
%curves_numContainers = zeros(0,1);
%curves_containers = zeros(0,2);

pol_u_checked = false(numPolylines_u,1);
for iregion = (mesh.groundRegion+1):numRegions
    
    %region_poly_ids_i = find(polylines_to_region==iregion);%region_polylines{iregion};%
    region_poly_ids_i = region_polylines{iregion};%find(polylines_to_region==iregion);%
    region_poly_ids_i_u = polylines_map_to_unique(region_poly_ids_i);
    region_poly_ids_i_u = region_poly_ids_i_u(find(pol_u_checked(region_poly_ids_i_u)==0));
    
    for jregion = (iregion+1):numRegions
        
        %region_poly_ids_j = find(polylines_to_region==jregion);%region_polylines{jregion};%
        region_poly_ids_j = region_polylines{jregion};%find(polylines_to_region==jregion);%
        region_poly_ids_j_u = polylines_map_to_unique(region_poly_ids_j);
        region_poly_ids_j_u = region_poly_ids_j_u(find(pol_u_checked(region_poly_ids_j_u)==0));
        
        diffLines_ij_u = setdiff(region_poly_ids_i_u,region_poly_ids_j_u);
        eqLines_ij_u = setdiff(region_poly_ids_i_u,diffLines_ij_u);
        
        if(not(isempty(eqLines_ij_u)))
            numCurves = numCurves+1;
            curves{numCurves} = eqLines_ij_u';
            %curves_numContainers(numCurves) = 2;
            %curves_containers(numCurves,:) = [iregion,jregion];
            
            pol_u_checked(eqLines_ij_u) = true; 
            
            region_poly_ids_i_u = region_poly_ids_i_u(find(pol_u_checked(region_poly_ids_i_u)==0));
        end
    end 
    
    if(not(isempty(region_poly_ids_i_u)))
        % a unique curve that can be non-connected
        numCurves = numCurves+1;
        curves{numCurves} = region_poly_ids_i_u;
        %curves_numContainers(numCurves) = 1;
        %curves_containers(numCurves,:) = [iregion];

        pol_u_checked(region_poly_ids_i_u) = 1; 
     end     
end

%% Split curves into straight pieces
angleSplit = 10;
angleSplit_rad = angleSplit*(pi/180);
splitFactor = cos(angleSplit_rad);
num_aligned_curves = 0; %at least there will be numCurves
aligned_curves = cell(numCurves,1);
%aligned_curves_numContainers = zeros(numCurves,1);
%aligned_curves_containers = zeros(numCurves,2);
for icurve = 1:numCurves

    % separe the curves into sets that are linearly dependent
    poly_ids_i_u = curves{icurve};

    edges = polylines_unique(poly_ids_i_u,:);
    v_e = mesh.X(edges(:,2),:)-mesh.X(edges(:,1),:);
    v_e_n = v_e./sqrt(sum(v_e.*v_e,2));
    while(not(isempty(v_e_n)))
        v1 = v_e_n(1,:);
        
        v_dot_edges = abs(v_e_n(:,1).*v1(1) + v_e_n(:,2).*v1(2));
        edges_angle = find(v_dot_edges>splitFactor);

        num_aligned_curves = num_aligned_curves+1;
        aligned_curves{num_aligned_curves} = poly_ids_i_u(edges_angle);
        %aligned_curves_numContainers(num_aligned_curves) = curves_numContainers(icurve);
        %aligned_curves_containers(num_aligned_curves,:) = curves_containers(icurve,:);
        
        if(isempty(edges_angle))
            angle
            angleSplit
            error('not possible.. at least should find himself')
        end
        
        poly_ids_i_u(edges_angle) = [];
        edges(edges_angle,:) = [];
        v_e_n(edges_angle,:) = [];
    end    

    % build output polylines structure
    %polylines_out.T = [polylines_out.T ; polylines_unique(connected_curves{icurve},:)]
end
numCurves = num_aligned_curves;
curves = aligned_curves;
%curves_numContainers = aligned_curves_numContainers;
%curves_containers = aligned_curves_containers;

%% Corbes en components connexes
num_connected_curves = 0; %at least there will be numCurves
connected_curves = cell(numCurves,1);
%connected_curves_numContainers = zeros(numCurves,1);
%connected_curves_containers = zeros(numCurves,2);
for icurve =1:numCurves
    
    poly_ids_i_u = curves{icurve};
    polylines_in_curve = polylines_unique(poly_ids_i_u,:);

    % find endpoints
    points = polylines_in_curve(:);
    point_count = sparse(points,ones(length(points),1),ones(length(points),1));
    endPoints = find(point_count==1);
    
    checkPoints = find(point_count>2);
    if(not(isempty(checkPoints)))
        polylines_in_curve
        checkPoints
        error('not possible.. a point should appear just once or twice in a polycurve')
    end

    count = 0;
    maxCount = 1000;
    if(not(isempty(endPoints))) % it is not a simply connected curve

        while(not(isempty(endPoints)))

            connectedCurve_ids_u = [];

            prevPoint = endPoints(1);
            endPoints(1) = [];
            doContinue = true;
            while(doContinue && maxCount>count)
                count = count+1;

                [edge,vertex]=find(polylines_in_curve==prevPoint);
                line = poly_ids_i_u(edge);

                if(isempty(edge))
                    endPoint_idToRemove = find(endPoints==prevPoint);
                    endPoints(endPoint_idToRemove) = [];
                    doContinue = false;
                else
                    connectedCurve_ids_u = [connectedCurve_ids_u; line];
                    if(vertex==1)
                        prevPoint = polylines_in_curve(edge,2);
                    else
                        prevPoint = polylines_in_curve(edge,1);
                    end
                    polylines_in_curve(edge,:) = [];
                    poly_ids_i_u(edge) = [];
                end
            end
            num_connected_curves = num_connected_curves+1;
            connected_curves{num_connected_curves} = connectedCurve_ids_u;
            theRegions = unique(polylines_to_region(polylines_mantained(connectedCurve_ids_u)));
            numReg = length(theRegions);
            if(numReg>2) 
                theRegions
                error('the connected curve should be only in one or two buildings')
                % think if it can even be more than one...
            end
            %connected_curves_numContainers(num_connected_curves) = numReg;
            %connected_curves_containers(num_connected_curves,1:numReg) = theRegions;

            pol_u_checked(connectedCurve_ids_u) = 1; 
        end
        
    else % it was already a connected curve
        num_connected_curves = num_connected_curves+1;
        connected_curves{num_connected_curves} = curves{icurve};
        %connected_curves_numContainers(num_connected_curves) = curves_numContainers(icurve);
        %connected_curves_containers(num_connected_curves,:) = curves_containers(icurve,:);
    end
    
end
numCurves = num_connected_curves;
curves = connected_curves;
%curves_numContainers = connected_curves_numContainers;
%curves_containers = connected_curves_containers;

%% Simplify straight curves (all) by two points
doSimplifyStraightCurves = true;
if(doSimplifyStraightCurves)
%     num_straight_curves = 0; %at least there will be numCurves
%     straight_curves = cell(numCurves,1);
%     straight_curves_numContainers = zeros(numCurves,1);
%     straight_curves_containers = zeros(numCurves,2);

    polylines_out.T = zeros(numCurves,2);

    for icurve = 1:numCurves 

        poly_ids_i_u = curves{icurve};
        polylines_in_curve = polylines_unique(poly_ids_i_u,:);

        % find endpoints
        points = polylines_in_curve(:);
        point_count = sparse(points,ones(length(points),1),ones(length(points),1));
        endPoints = find(point_count==1); 
        
        if(length(endPoints)>2)
            error('should have been already splitted')
        end
        
        polylines_out.T(icurve,:) = endPoints;
    end
else
    listCurves = cell2mat(curves);
    if(length(listCurves)==numPolylines_u)
        %length(listCurves)
        %numPolylines_u
    else
        error('something wrong')
    end
    polylines_out.T = polylines_unique(listCurves,:);

end

%% build output polylines structure
finalPoints = unique(polylines_out.T);
numPoints = length(finalPoints);

mapToRemaining = zeros(size(mesh.X,1),1);
mapToRemaining(finalPoints) = 1:numPoints;

polylines_out.T(:,:) = mapToRemaining(polylines_out.T);
polylines_out.X = mesh.X(finalPoints,1:2);

doCollapse=false;
doCloseRegions=false;
countLoops = countLoops+1;

% curvesPlot = connected_curves;
% curvesPlot_num = num_connected_curves;
% % curvesPlot = aligned_curves;
% % curvesPlot_num = num_aligned_curves;
% figure(15); clf
% plot(polylines_out.X(:,1),polylines_out.X(:,2),'*')
% for icurve = 1:curvesPlot_num
%     %curves{icurve}
%     plot_points = polylines_unique(curvesPlot{icurve},:);
%     figure(15)
%     hold on;
%     if(size(plot_points,1)==1)
%         plot(mesh.X(plot_points(:),1),mesh.X(plot_points(:),2),'--')
%     else
%         theRandColor = rand(1,3);
%         for iline=1:size(plot_points,1)
%             plot(mesh.X(plot_points(iline,:),1),mesh.X(plot_points(iline,:),2),'-.','color',theRandColor)
%         end
%     end
% pause()
% end
% axis equal
% pause()
% disp('---->>>> EL DIBUIX ESTA A COLLAPSEBOUNDARYOFREGIONS_BUILGINDS')
%% basura

% edgesSharedByRegions = cell(numRegions,numRegions);
% 
% for iregion = (mesh.groundRegion+1):numRegions
%   
%     for jregion = (iregion+1):numRegions
%         
%         
%     
%     end
% end

%%


%             region_poly_ids_i_u = diffLines_ij_u;
%             
%             %diffLines = setdiff(region_poly_ids_j_u,eqLines);
%             %eqLines = setdiff(1:length(region_poly_ids_i),region_poly_ids_i);
%             
%             notFound_j_u = setdiff(region_poly_ids_j_u,eqLines_ij_u);
%             found_j_u = setdiff(region_poly_ids_j_u,notFound_j_u);
%             %region_poly_ids_j = setdiff(region_poly_ids_j_u,found_j_u
%             polylines_to_region() = -1;

%%

%     if(curves_numContainers(icurve) == 2)
%        
%         polylines_in_curve = polylines_unique(curves{icurve},:);
%         % find endpoints
%         points = polylines_in_curve(:);
%         endPoints = find(sparse(points,ones(length(points),1),ones(length(points),1))==1);
%         % compose line between endpoints
%         if(length(endPoints)==2)
%             polylines_unique = [polylines_unique; endPoints'];
%             curves{icurve} = size(polylines_unique,1);
%             
% %             figure(15)
% %             hold on;
% %             plot(mesh.X(endPoints,1),mesh.X(endPoints,2),'-')
%         elseif(length(endPoints)==size(polylines_in_curve,1))
%             % do nothing, the building is isolated and is a block itself
%         else
%             count = 0;
%             maxCount = 1000;
%             while(not(isempty(endPoints)))
%                 prevPoint = endPoints(1);
%                 endPoints(1) = [];
%                 firstVertex = prevPoint;
%                 doContinue = true;
%                 while(doContinue && maxCount>count)
%                     count = count+1;
%                     
%                     [line,vertex]=find(polylines_in_curve==prevPoint);
%                     
%                     if(isempty(line))
%                         doContinue = false;
%                     else
%                         if(vertex==1)
%                             prevPoint = polylines_in_curve(line,2);
%                         else
%                             prevPoint = polylines_in_curve(line,1);
%                         end
%                         polylines_in_curve(line,:) = [];
%                     end
%                 end
%                 secondVertex = prevPoint;
%                 
%                 polylines_unique = [polylines_unique; [firstVertex,secondVertex]];
%                 curves{icurve} = size(polylines_unique,1);
%                 
%             end
%             
%             if(count>=maxCount)
%                 polylines_in_curve            
%                 endPoints
% 
%                 figure(18)
%                 hold on;
%                 plot(mesh.X(endPoints,1),mesh.X(endPoints,2),'*')
%                 for iline = 1:size(polylines_in_curve,1)
%                     plot(mesh.X(polylines_in_curve(iline,:),1),mesh.X(polylines_in_curve(iline,:),2),'-')
%                 end
%                 error('some problem')
%             end
%         end
%     end







