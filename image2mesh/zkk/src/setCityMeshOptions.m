function [domainLimits,topHeightOverTerrain,edgeLengthGround,edgeLengthCeil,...
    edgeLengthGround_coarse,meshSize2D_ground,meshAngle2D,...
    meshSize2D_ground_coarse,parameters_tet,elementOptions,...
    tolHeight,minElevFromGround,tolFacadePoints,tolCloseRegions,minElementAreaAllowed,...
    doSplitBoundaryTets,doOptimizeMesh,doTreatSpecialBuildings,typeTreatSpecialBuildings,...
    idealizationLevel,edgeMarks,meshFacadesWithQuads,meshNotClassificableRegionsWithLidar]...
 =setCityMeshOptions(...
  topHeightOverTerrain,edgeLengthGround,edgeLengthCeil,edgeLengthGround_coarse)
%% Code options
%-%-%-%-%-%-%-%-%-%-%-%-Meshing domain-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-
domainLimits.type = 'cadastre';%'cadastre';'dem';'lidar';'intersection';'reset';
% domainLimits.xmin = ;domainLimits.xmax = ;domainLimits.ymin = ;domainLimits.ymax = ;

domainLimits.margin      =50%100;% 10; %500 %0==detault, margin not zero is number of meters of margin
% domainLimits.marginLeft  = 5000%200; %100; 
% domainLimits.marginRight = 9000%200; %200; 
% domainLimits.marginUp    = 1000%200; %75; 
% domainLimits.marginDown  = 6000%400; %75;

domainLimits.sizeBox = false; % surface discretized using fine mesh in city, coarse outside the city
domainLimits.typeOut = 'lidar';%'cadastre';'dem';'lidar';'intersection';'reset';
% domainLimits.xmin = 2556;
% domainLimits.ymin = 1210;
% domainLimits.xmax = 3540;
% domainLimits.ymax = 2278;

%-%-%-%-%-%-%-%-%-%-%-%-Mesh sizing-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% %edgeLengthBuildings = XXX ==> i fer que els nodes de facana tinguin aquesta mida mes fina

factor_srf_volume = 1%2;%4%2 %num of divisions of the surface elements to generate volume mesh
meshSize2D_ground = edgeLengthGround*edgeLengthGround/2.0;
meshSize2D_ground = meshSize2D_ground*factor_srf_volume^2; % the volume mesh will be the surface/2.0

meshAngle2D = 30; %desired minimum angle in the 2D mesh

doSplitBoundaryTets = false; % if a tet has all nodes in boundary, split it to avoid this
doOptimizeMesh = false;
%-%-%-%-%-%-%-%-%-%-%-%-Meshing tolerance-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
tolHeight         = 5; %absolute tolerance between lidar and topo to consider a catastral region a building
minElevFromGround = 2.5; %minimum elevationpar to the ground at any point in the building to consider it building
tolFacadePoints   = min([5,edgeLengthGround/2.0]); %points in a facade that are closer than this distance will be removed (abans estava edge/2)
tolCloseRegions   = edgeLengthGround*2.0; %minimum size to not remove an entrance to a block
minElementAreaAllowed = meshSize2D_ground/10000; % a building with a smaller element will be removed
%% -%-%-%-%-%-%-%-%-%-%-%-Testing stuff-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
doTreatSpecialBuildings = false;
typeTreatSpecialBuildings = 'fineMesh';%'fineMesh','anisMesh','visualMesh'
%% -%-%-%-%-%-%-%-%-%-%-%-Fixed stuff-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
parameters_tet.size = true;
parameters_tet.edgeLengthCeil   = edgeLengthCeil;
parameters_tet.areaCeil   = 0.5*edgeLengthCeil*edgeLengthCeil*factor_srf_volume^2;%edgeLengthCeil;%meshSize2D_ceil;%100;
parameters_tet.factor_srf_volume = factor_srf_volume;
%parameters_tet.sizeGround = 10;%5;
%-%-%-%-%-%-%-%-%-%-%-%-Different sizings-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
meshSize2D_ground_coarse = 0;
if(domainLimits.sizeBox)
    meshSize2D_ground_coarse = 0.5*edgeLengthGround_coarse*edgeLengthGround_coarse;
    %meshSize2D_ground_coarse = meshSize2D_ground_coarse*factor_srf_volume^2;
end
%-%-%-%-%-%-%-%-%-%-%-%-Other options-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
idealizationLevel = 'block'; % "block" level , or "building" level? ==> right now only illa (block)
idealizationLevel = 'building'; % "block" level , or "building" level? ==> right now only illa (block)

elementOptions.orderQuadrature = 3;
elementOptions.quadrature = 'gaussStandard';
elementOptions.shapeFunctions = 'monomial';

global fieldIntegrationType; % DO NOT CHANGE GAUSSPOINTS
fieldIntegrationType = 'gaussPoints';%'pixels','gaussPoints';

edgeMarks = true;

meshFacadesWithQuads = false;
meshNotClassificableRegionsWithLidar = false;




