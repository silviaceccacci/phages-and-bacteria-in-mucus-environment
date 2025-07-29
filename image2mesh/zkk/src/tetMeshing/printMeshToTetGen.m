function [inpName]=printMeshToTetGen(mesh,fileName,parameters)

    regionsInGeometry = false;
    
    if(~regionsInGeometry)
        [inpName]=printMeshAsSeparateFacets(mesh,fileName);
    else
        [inpName]=printMeshWithGeometryFacets(mesh,fileName);
    end
    
    if(isfield(parameters,'size') && parameters.size)
        %metric description: tetgen manual page 70
        printMetric(mesh,fileName);
    end
    

end

function [inpName]=printMeshWithGeometryFacets(mesh,fileName)
    inpName = [fileName '.poly'];

    numNodes = size(mesh.X,1);
    
    stream = fopen(inpName,'w');
    
    % Nodes
    fprintf(stream,'%d 3 0 0\n',numNodes);
    fprintf(stream,'%d %f %f %f\n', [1:numNodes ;  mesh.X']);
    
    % Elements
    numFacets = length(mesh.facets);
    fprintf(stream,'%d 1\n',numFacets);
    for ifacet = 1:numFacets
        numPoligons = length(mesh.facets{ifacet});
        facetInfo = zeros(1,3);
        facetInfo(1) = numPoligons; % number of polygons
        facetInfo(2) = 0; % number of holes
        facetInfo(3) = mesh.facetBoundId(ifacet); % boundary markers
        fprintf(stream,'%d %d %d\n', facetInfo );
        fprintf(stream,' %d %d %d %d\n',...
            [3*ones(numPoligons,1) mesh.T(mesh.facets{ifacet},:)]');
    end
    
%     fprintf(stream,'%d 0\n',numElements);
%     fprintf(stream,'%d 0 0\n',numElements);
%     fprintf(stream,'%d %d %d %d\n', [3*ones(1,numElements); mesh.T']);

%     fprintf(stream,'%d 1\n',numElements);
%     extraInfo = zeros(numElements,3);
%     extraInfo(:,1) = 1; % number of polygons
%     extraInfo(:,2) = 0; % number of holes
%     extraInfo(:,3) = mesh.elementBoundId; % boundary markers
%     fprintf(stream,'%d %d %d\n %d %d %d %d\n', [extraInfo 3*ones(numElements,1) mesh.T]');

    % Holes
    %fprintf(stream,'%d\n',0);
    numHoles = size(mesh.holeBuildings,1)-1;
    fprintf(stream,'%d\n',numHoles);
    fprintf(stream,'%d %f %f %f\n', [1:numHoles ;  mesh.holeBuildings(2:end,:)']);
       
    fclose(stream);
   
end

function [inpName]=printMeshAsSeparateFacets(mesh,fileName)
    inpName = [fileName '.poly'];
    numNodes = size(mesh.X,1);
    numElements = size(mesh.T,1);
    
    stream = fopen(inpName,'w');
    % Nodes
    fprintf(stream,'%d 3 0 0\n',numNodes);
    fprintf(stream,'%d %f %f %f\n', [1:numNodes ;  mesh.X']);
    % Elements
    fprintf(stream,'%d 1\n',numElements);
    extraInfo = zeros(numElements,3);
    extraInfo(:,1) = 1; % number of polygons
    extraInfo(:,2) = 0; % number of holes
    extraInfo(:,3) = mesh.elementBoundId; % boundary markers
    fprintf(stream,'%d %d %d\n %d %d %d %d\n', [extraInfo 3*ones(numElements,1) mesh.T]');
    % Holes
    fprintf(stream,'0\n');
    
    fclose(stream);
end

function [inpName]=printMetric(mesh,fileName)
    inpName = [fileName '.mtr'];

    numNodes = size(mesh.X,1);
    
    stream = fopen(inpName,'w');
    
    fprintf(stream,'%d 1\n',numNodes);
    fprintf(stream,'%f\n', mesh.hnodes);

    fclose(stream);
end


%   gapName = 'Gap'
%   gapFile = './temp/'//TRIM(gapName)//'.poly'
% 	f1 = TRIM(path)//TRIM(gapFile)
% 	open(90,file=TRIM(f1),status='unknown')
% 	write(90,'(1(i8,1x),3(1x,i1))') ,numBoundPoints+numBoundPointsDisc,3,0,0 !the last 0 could be 1 to set different boundary markers
%   !**Part 1 - node list
%   countNode = 0
% 	do inode = 1,numBoundPoints
%     countNode = countNode+1
%     x(:) = boundTriCoords(:,inode) - translation(:)
%     !write(90,'(i8,3(1x,e30.20),1x,i1)') countNode,x(1),x(2),x(3),1
%     write(90,'(i8,3(1x,e30.20))') countNode,x(1),x(2),x(3)
% 	end do
% 	do inode = 1,numBoundPointsDisc
%     countNode = countNode+1
%     x(:) = boundTriCoordsDisc(:,inode) - translation(:)
%     !write(90,'(i8,3(1x,e30.20),1x,i1)') countNode,x(1),x(2),x(3),2
%     write(90,'(i8,3(1x,e30.20))') countNode,x(1),x(2),x(3)
% 	end do
%   !**Part 2 - facet list
%   countElem = 0
% 	write(90,'(1(i8,1x),3(1x,i1))') ,numBoundTri+numBoundTriDisc,0!the last 0 could be 1 to set different boundary markers
% 	do ielem = 1,numBoundTri
%     countElem = countElem+1
%     triFace = boundTriMesh(:,ielem)
%     write(90,'(i1,1x,i1,1x,i1)') 1,0,0
% 		!write(90,'(i1,3(1x,i8),1x,i1)') countElem,(triFace(inode),inode=1,3),2
% 		write(90,'(i1,3(1x,i8),1x,i1)') 3,(triFace(inode),inode=1,3)
% 	end do      
% 	do ielem = 1,numBoundTriDisc
%     countElem = countElem+1
%     triFace = boundTriMeshDisc(:,ielem)+numBoundPoints
%     write(90,'(i1,1x,i1,1x,i1)') 1,0,0
% 		!write(90,'(i1,3(1x,i8),1x,i1)') countElem,(triFace(inode),inode=1,3),2
% 		write(90,'(i1,3(1x,i8),1x,i1)') 3,(triFace(inode),inode=1,3)
% 	end do        
% 	!**Part 3 - hole list
% 	write(90,'(1(i8,1x))') ,numInsertions
% 	do ibox = 1,numInsertions
%     x(:) = cmInsertedDiscs(:,ibox) - translation(:)
%     write(90,'(i8,3(1x,e30.20))') ibox,x(1),x(2),x(3)
% 	end do
%   !**Part 4 - region attributes list (OPTIONAL)
%   !One line: <# of region>
%   !Following lines list # of region attributes:
%   !<region #> <x> <y> <z> <region number> <region attribute>
% 	close(90)		