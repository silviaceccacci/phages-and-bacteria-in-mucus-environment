function [polyBox,polyBox2D,domainLimits]=addBoxTo1DGeometry(...
    polylines,fields,translationVector,domainLimits,outScreen,edgeLength)

    if(nargin==4)
        outScreen = true;
    end

    polyBox = polylines;

    X = polylines.X;
    
    if(nargin==1)
    %    if(true)
        xmin = min(X(:,1));
        xmax = max(X(:,1));
        ymin = min(X(:,2));
        ymax = max(X(:,2));

        equiSpacing = sqrt( (xmax-xmin)*(ymax-ymin) / size(polylines.T,1) );

        %spacingFactor = 10;%0.01;%0.01;%10;
        %spacingX = spacingFactor*equiSpacing;
        %spacingY = spacingFactor*equiSpacing;
        spacingX = edgeLength;
        spacingY = edgeLength;

        
        if(isfield(domainLimits,'margin') && domainLimits.margin>0)
            widthX = domainLimits.margin;
            widthY = domainLimits.margin;
        elseif(isfield(domainLimits,'marginLeft') && domainLimits.marginLeft>0)
            widthLeft  = domainLimits.marginLeft;
            widthRight = domainLimits.marginRight;
            widthUp    = domainLimits.marginUp;
            widthDown  = domainLimits.marginDown;
        else
            widthFactor = 2;
            widthX = widthFactor*spacingX;
            widthY = widthFactor*spacingY;
%             widthX = edgeLength;
%             widthY = edgeLength;
        end
        
        if(isfield(domainLimits,'marginLeft'))
            xmin = min(X(:,1))-widthLeft;
            xmax = max(X(:,1))+widthRight;
            ymin = min(X(:,2))-widthDown;
            ymax = max(X(:,2))+widthUp;
        else
            xmin = min(X(:,1))-widthX;
            xmax = max(X(:,1))+widthX;
            ymin = min(X(:,2))-widthY;
            ymax = max(X(:,2))+widthY;
        end
    else
        switch domainLimits.type
            case 'cadastre'
                xmin = min(X(polylines.T(:),1));
                xmax = max(X(polylines.T(:),1));
                ymin = min(X(polylines.T(:),2));
                ymax = max(X(polylines.T(:),2));

                equiSpacing = sqrt( (xmax-xmin)*(ymax-ymin) / size(polylines.T,1) );

               % disp('NO MARGINNNNNNSNSNSNSNS')
                spacingFactor = 1;%10;%0.01;%0.01;%10;
                spacingX = spacingFactor*equiSpacing;
                spacingY = spacingFactor*equiSpacing;

        if(isfield(domainLimits,'margin') && domainLimits.margin>0)
            widthX = domainLimits.margin;
            widthY = domainLimits.margin;
        elseif(isfield(domainLimits,'marginLeft') && domainLimits.marginLeft>0)
            widthLeft  = domainLimits.marginLeft;
            widthRight = domainLimits.marginRight;
            widthUp    = domainLimits.marginUp;
            widthDown  = domainLimits.marginDown;
        else
            widthFactor = 2;
            widthX = widthFactor*spacingX;
            widthY = widthFactor*spacingY;
%             widthX = edgeLength;
%             widthY = edgeLength;
        end
        
        if(isfield(domainLimits,'marginLeft'))
            xmin = min(X(:,1))-widthLeft;
            xmax = max(X(:,1))+widthRight;
            ymin = min(X(:,2))-widthDown;
            ymax = max(X(:,2))+widthUp;
        else
            xmin = min(X(:,1))-widthX;
            xmax = max(X(:,1))+widthX;
            ymin = min(X(:,2))-widthY;
            ymax = max(X(:,2))+widthY;
        end
%                 if(isfield(domainLimits,'margin') && domainLimits.margin>0)
%                     widthX = domainLimits.margin;
%                     widthY = domainLimits.margin;
%                 else
%                     widthFactor = 2;
%                     widthX = widthFactor*spacingX;
%                     widthY = widthFactor*spacingY;
%                 end
%                 
%                 xmin = min(X(:,1))-widthX;
%                 xmax = max(X(:,1))+widthX;
%                 ymin = min(X(:,2))-widthY;
%                 ymax = max(X(:,2))+widthY;
%                 xmin = min(polylines.X(:,1));
%                 ymin = min(polylines.X(:,2));
%                 xmax = max(polylines.X(:,1));
%                 ymax = max(polylines.X(:,2));
% 
%                 marginFactor = 1000.0;
%                 
%                 lx = xmax-xmin;
%                 ly = ymax-ymin;
% %                 xmin = xmin - lx/marginFactor;
% %                 ymin = ymin - ly/marginFactor;
% %                 xmax = xmax + lx/marginFactor;
% %                 ymax = ymax + ly/marginFactor;
%                 xmin = xmin - 20;
%                 ymin = ymin - 20;
%                 xmax = xmax + 20;
%                 ymax = ymax + 20;
% 
%                 spacingX = lx/5;
%                 spacingY = ly/5; 
            case 'dem'
                field = fields.ground;
                
                if(isfield(domainLimits,'margin') && domainLimits.margin>0)
                    myMargin = domainLimits.margin;
                else
                    myMargin = edgeLength;
                end
                
                xmin = field.x0(1)-translationVector(1);
                ymin = field.x0(2)-translationVector(2);
                xmax = xmin + field.nx*field.hx - myMargin;
                ymax = ymin + field.ny*field.hy - myMargin;
                xmin = xmin + myMargin;
                ymin = ymin + myMargin;

                spacingX = edgeLength;%200;
                spacingY = edgeLength;%200; 
                
                
            case 'lidar'
                field = fields.roof;
                
                if(isfield(domainLimits,'margin') && domainLimits.margin>0)
                    myMargin = domainLimits.margin;
                else
                    myMargin = edgeLength;
                end
                
                xmin = field.x0(1)-translationVector(1);
                ymin = field.x0(2)-translationVector(2);
                xmax = xmin + field.nx*field.hx - myMargin;
                ymax = ymin + field.ny*field.hy - myMargin;
                xmin = xmin + myMargin;
                ymin = ymin + myMargin;

                spacingX = edgeLength;%200;
                spacingY = edgeLength;%200;  
            case 'intersection'
                %polyline limits
                xmin = min(X(polylines.T(:),1));
                xmax = max(X(polylines.T(:),1));
                ymin = min(X(polylines.T(:),2));
                ymax = max(X(polylines.T(:),2));
                equiSpacing = sqrt( (xmax-xmin)*(ymax-ymin) / size(polylines.T,1) );
                spacingFactor = 2;%10     %0.01;%0.01;%10;
                spacingX = spacingFactor*equiSpacing;
                spacingY = spacingFactor*equiSpacing;
                if(isfield(domainLimits,'margin') && domainLimits.margin>0)
                widthX = domainLimits.margin;
                widthY = domainLimits.margin;
                else
                widthFactor = 2;
                widthX = widthFactor*spacingX;
                widthY = widthFactor*spacingY;
                end
                xmin_c = min(X(:,1))-widthX;
                xmax_c = max(X(:,1))+widthX;
                ymin_c = min(X(:,2))-widthY;
                ymax_c = max(X(:,2))+widthY;
                %lidar limits
                field = fields.roof;
                myMargin = edgeLength;%10;
                xmin_l = field.x0(1)-translationVector(1);
                ymin_l = field.x0(2)-translationVector(2);
                xmax_l = xmin_l + field.nx*field.hx - myMargin;
                ymax_l = ymin_l + field.ny*field.hy - myMargin;
                xmin_l = xmin_l + myMargin;
                ymin_l = ymin_l + myMargin;
                % intersection
                xmin = max([xmin_l xmin_c]);
                ymin = max([ymin_l ymin_c]);
                xmax = min([xmax_l xmax_c]);
                ymax = min([ymax_l ymax_c]);                
                spacingX = edgeLength;%200;
                spacingY = edgeLength;%200;  
            case 'reset'
                xmin = domainLimits.xmin;
                xmax = domainLimits.xmax;
                ymin = domainLimits.ymin;
                ymax = domainLimits.ymax;

                spacingX = edgeLength;%domainLimits.hbound;
                spacingY = edgeLength;%domainLimits.hbound;  
                
            otherwise
                error('domain not set')
        end
    end
    
    if(outScreen)
        fprintf('   ...Domain: [%f,%f]x[%f,%f]\n',xmin,xmax,ymin,ymax)
    end
    
    nx = ceil((xmax-xmin)/spacingX);
    spacingX = (xmax-xmin)/nx - 1e-12;
    ny = ceil((ymax-ymin)/spacingY);
    spacingY = (ymax-ymin)/ny - 1e-12;
    
    Xpoints = xmin:spacingX:xmax; Xpoints = Xpoints';
    Ypoints = ymin:spacingY:ymax; Ypoints = Ypoints';
    
    mx = length(Xpoints)-1;
    my = length(Ypoints)-1;
    Xbox = [Xpoints(1:mx)              Ypoints(1)*ones(mx,1)
            Xpoints(mx+1)*ones(my,1)   Ypoints(1:my)
            Xpoints((mx+1):-1:2)       Ypoints(my+1)*ones(mx,1)
            Xpoints(1)*ones(my,1)      Ypoints((my+1):-1:2)];
        
	polyBox2D.X = Xbox;
         
	T = [ [1:size(Xbox,1)]' [2:size(Xbox,1) 1]'];
    
    polyBox2D.T = T;
        
    T = T + size(X,1);
    
    polyBox2D.globIdsXbox = (size(X,1)+1):(size(X,1)+size(Xbox,1));
    
    polyBox.X = [ X; Xbox];
    polyBox.T = [polylines.T; T];
    
    buildingMark = 2;
    boundaryMark = 3;
    edgeMarks = [buildingMark*ones(size(polylines.T,1),1) ; boundaryMark*ones(size(T,1),1)];
    
    polyBox.buildingMark = buildingMark;
    polyBox.boundaryMark = boundaryMark;
    polyBox.mark = edgeMarks;
    polyBox2D.buildingMark = buildingMark;
    polyBox2D.boundaryMark = boundaryMark;
    polyBox2D.mark = edgeMarks;
    
    
    domainLimits.type = 'reset';
    domainLimits.xmin = xmin;
    domainLimits.xmax = xmax;
    domainLimits.ymin = ymin;
    domainLimits.ymax = ymax;
    domainLimits.hbound = (spacingX+spacingY)/2.0;
    
    
end





