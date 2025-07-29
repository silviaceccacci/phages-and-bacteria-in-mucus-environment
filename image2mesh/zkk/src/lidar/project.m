function [z_on_x] = project(field, target, options)  

    if(field.structured)
        [z_on_x] = project_structured(field, target, options);
    else
        [z_on_x] = project_toMesh(field, target, options);
    end
    
end


