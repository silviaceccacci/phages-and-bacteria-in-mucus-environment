function [mesh_facade3D]=projectingMeshToField(type,mesh,field,fakeBuildingHeight)
if(nargin==3)
    fakeBuildingHeight=0.0;
end

if(strcmp(type,'nodal'))
    [mesh_facade3D]=projectingMeshToField_nodal(mesh,field,fakeBuildingHeight);
elseif(strcmp(type,'elemental'))
    [mesh_facade3D]=projectingMeshToField_elemental(mesh,field,fakeBuildingHeight);     
elseif(strcmp(type,'elementalStructured'))
    [mesh_facade3D]=projectingMeshToField_elementalStructured(mesh,field,fakeBuildingHeight);        
else
    error('not implemented')        
end