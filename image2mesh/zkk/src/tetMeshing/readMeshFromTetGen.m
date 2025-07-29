function [mesh]=readMeshFromTetGen(fileName,faceSwitch)
%% Nodes
nodeDat = [fileName '.1.node'];
fidnode = fopen(nodeDat, 'r');

outline = fscanf(fidnode,'%d %d %d %d\n',4);
numNodes = outline(1);
X = fscanf(fidnode, '%f %f %f %f\n', [4 numNodes])';
X = X(:,2:4);

fclose(fidnode);

%% Elems
elemDat = [fileName '.1.ele'];
fidelem = fopen(elemDat, 'r');

outline = fscanf(fidelem,'%d %d %d\n',3);
numElems = outline(1);
T = fscanf(fidelem, '%d %d %d %d %d\n', [5 numElems])';
T = T(:,2:5);

fclose(fidelem);

%% Set mesh structure
mesh.T = T;
mesh.X = X;

element.dim = 3;
element.numCoord = 3;
element.type = 'tet';
element.order = 1;
element.distribution = 'lineal';
mesh.element = defineElement(element);

%% Set mesh boundary
nodeDat = [fileName '.1.face'];
fidnode = fopen(nodeDat, 'r');

if(isempty(faceSwitch))
    outline = fscanf(fidnode,'%d %d\n',2);
    numBfaces = outline(1);
    isThereBMarker = outline(2);
    if(isThereBMarker~=1) 
        error('No boundary markers... need them to generate boundary conditions')
    end
    boundInfo = fscanf(fidnode, '%d %d %d %d %d %d %d\n', [7 numBfaces])';
else
    outline = fscanf(fidnode,'%d %d\n',2);
    numFaces = outline(1);
    isThereBMarker = outline(2);
    if(isThereBMarker~=1) 
        error('No boundary markers... need them to generate boundary conditions')
    end
    faceInfo = fscanf(fidnode, '%d %d %d %d %d %d %d\n', [7 numFaces])';
    boundFaces = find(faceInfo(:,5)>0);
    boundInfo = faceInfo(boundFaces,:);
    
    faceInfo(boundFaces,:) = [];
    faceNodes = sort(faceInfo(:,2:4),2);
    faceElements = faceInfo(:,6:7);
    mesh.innerFaces.elements = faceElements;
    mesh.innerFaces.nodes = faceNodes;
end
boundaries.faces = boundInfo(:,2:4);
boundaries.marks = boundInfo(:,5);
boundaries.elems = max(boundInfo(:,6:7),[],2);

mesh.boundaries = boundaries;
mesh.boundaryNodes = unique(mesh.boundaries.faces(:));

fclose(fidnode);

% meshBoundary.name = 'boundary';
% meshBoundary.X = mesh.X;
% meshBoundary.T = boundaries.faces;
% meshBoundary.elementField = boundaries.marks;
% exportTriMeshToParaview(meshBoundary);

end

% 	open(91,FILE=TRIM(path)//'temp/'//TRIM(gapName)//'.1.node',STATUS='old',ERR=101)
% 	open(92,FILE=TRIM(path)//'temp/'//TRIM(gapName)//'.1.ele',STATUS='old',ERR=101)
%   !** Read nodes
% 	read(91,*) numNodeTetgen,kk1,kk2,kk3
% 	allocate(nodeDelTemp(3,numNodeTetgen))	
% 	do inode = 1,numNodeTetgen
% 		read(91,*)  kk1,nodeDelTemp(1,inode),nodeDelTemp(2,inode),nodeDelTemp(3,inode)
%     nodeDelTemp(:,inode) = nodeDelTemp(:,inode) + translation
% 	end do	
%   !** Read elements
% 	read(92,*) numElemTetgen,kk1,kk2
% 	numTetDel = numElemTetgen
% 	allocate(tetrahedraDel(4,numTetDel))
% 	badMesh = .false.
% 	tetId = 0
% 	elemLoop: do ielem= 1,numTetDel
% 		read(92,*) tetId,tetrahedraDel(1,ielem),tetrahedraDel(2,ielem),tetrahedraDel(3,ielem),tetrahedraDel(4,ielem)
% 		tetNodes(:,1) = nodeDelTemp(:,tetrahedraDel(1,ielem))
% 		tetNodes(:,2) = nodeDelTemp(:,tetrahedraDel(2,ielem))
% 		tetNodes(:,3) = nodeDelTemp(:,tetrahedraDel(3,ielem))
% 		tetNodes(:,4) = nodeDelTemp(:,tetrahedraDel(4,ielem))
% 		determinant = computeDet(tetNodes)
% 		!*** Check mesh
% 		minDet = min(minDet,determinant)
% 		maxDet = max(maxDet,determinant)
% 		meanDet = meanDet+determinant
% 		if(determinant<tolDet) then
% 			badMesh = .true.
% 			tetId = ielem
% 			deallocate(tetrahedraDel)
% 			deallocate(nodeDelTemp)
% 			exit elemLoop
% 		end if 
% 	end do elemLoop
% 	close(91)	
% 	close(92)	



