function [logic] = checkLinearFace(element,X)
% input: X is the coordinates of an element
    logic=true;
    for iEdge = 1:element.numEdges
        Xedge = X(:,getEdge(1:element.numNod,element,iEdge));
        logic = min([logic checkLinearEdge(Xedge)]);
        if(logic==false)
            break;
        end
    end
end