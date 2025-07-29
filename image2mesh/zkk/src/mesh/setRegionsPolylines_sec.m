function [polylines]=setRegionsPolylines(polylines)
    
    segments = polylines.T;
    
    numSegments = size(segments,1);
    
    
    segmentRegion =zeros(numSegments,1);
    numRegionSegments = zeros(numSegments,1);
    regionSegments = zeros(numSegments,numSegments);
    
	iseg = 1;
    firstVertexFacade = segments(iseg,1);
    previousPoint = segments(iseg,2); 
    countClassifiedSegments = 0;
    
    iregion = 1;
%     numRegionSegments(iregion) = 1;
%     segmentRegion(iseg) = iregion;
%     regionSegments(iregion,numRegionSegments(iregion)) = iseg;

    while(iseg<=numSegments)
       
        %iseg = iseg+1;
        
        while(iseg<=numSegments && segments(iseg,2)~=firstVertexFacade)
            
            segmentRegion(iseg) = iregion;
            numRegionSegments(iregion) = numRegionSegments(iregion)+1;
            regionSegments(iregion,numRegionSegments(iregion)) = iseg; 
            
            %figure(1); hold on;
            %plot(polylines.X(segments(iseg,1:2),1),polylines.X(segments(iseg,1:2),2))
            
            iseg = iseg+1;
            
            %[iregion iseg]
            %[segments(iseg,:) segments(iseg+1,:)]
            %if(iregion>1) pause(1)
        end
        iregion = iregion+1;

        if(iseg<=numSegments) 
            firstVertexFacade = segments(iseg,1);
        end
        
%         segmentRegion(iseg) = iregion;
%         numRegionSegments(iregion) = numRegionSegments(iregion)+1;
%         regionSegments(iregion,numRegionSegments(iregion)) = iseg; 
       
        %iregion
        %iseg
        
        
    end
    iregion
    
%     iseg = 1;
%     firstVertexFacade = segments(iseg,1);
%     previousPoint = segments(iseg,2); 
%     countClassifiedSegments = 0;
%     
%     iregion = 1;
% %     numRegionSegments(iregion) = 1;
% %     segmentRegion(iseg) = iregion;
% %     regionSegments(iregion,numRegionSegments(iregion)) = iseg;
% 
%     while(iseg~=numSegments)
%        
%         %iseg = iseg+1;
%         
%         while(segments(iseg,2)~=firstVertexFacade)
%             
%             segmentRegion(iseg) = iregion;
%             numRegionSegments(iregion) = numRegionSegments(iregion)+1;
%             regionSegments(iregion,numRegionSegments(iregion)) = iseg; 
%             
%             %figure(1); hold on;
%             %plot(polylines.X(segments(iseg,1:2),1),polylines.X(segments(iseg,1:2),2))
%             
%             iseg = iseg+1;
%             
%             %[segments(iseg,:) segments(iseg+1,:)]
%             %if(iregion>1) pause(1)
%         end
%         pause()
%         iregion = iregion+1;
%         firstVertexFacade = segments(iseg,1);
%         
% %         segmentRegion(iseg) = iregion;
% %         numRegionSegments(iregion) = numRegionSegments(iregion)+1;
% %         regionSegments(iregion,numRegionSegments(iregion)) = iseg; 
%        
%         iregion
%         iseg
%         
%     end

    polylines.numRegions = iregion-1;
    polylines.segmentRegion = segmentRegion;
    polylines.numRegionSegments = numRegionSegments;
    polylines.regionSegments = regionSegments;

end


% % %     numRegionSegments(iregion) = 1;
% % %     segmentRegion(iseg) = iregion;
% % %     regionSegments(iregion,numRegionSegments(iregion)) = iseg;
% % 
% %     while(iseg~=numSegments)
% %        
% %         %iseg = iseg+1;
%         
% %         while(segments(iseg,2)~=firstVertexFacade)
% %             
% %             segmentRegion(iseg) = iregion;
% %             numRegionSegments(iregion) = numRegionSegments(iregion)+1;
% %             regionSegments(iregion,numRegionSegments(iregion)) = iseg; 
% %             
% %             %figure(1); hold on;
% %             %plot(polylines.X(segments(iseg,1:2),1),polylines.X(segments(iseg,1:2),2))
% %             
% %             iseg = iseg+1
% %         end
%     while(countClassifiedSegments~=numSegments)
%     
%         while(segments(iseg,2)~=firstVertexFacade)
%             if(iseg==1 || segments(iseg,1)==previousPoint)
%                 segmentRegion(iseg) = iregion;
%                 numRegionSegments(iregion) = numRegionSegments(iregion)+1;
%                 regionSegments(iregion,numRegionSegments(iregion)) = iseg; 
%                 previousPoint = segments(iseg,2);
%                 countClassifiedSegments = countClassifiedSegments+1;
%             end
%             
%             %figure(1); hold on;
%             %plot(polylines.X(segments(iseg,1:2),1),polylines.X(segments(iseg,1:2),2))
%             
%             iseg = iseg+1
%         end
%         iregion = iregion+1;
%         firstVertexFacade = segments(iseg,1);
%         
% %         segmentRegion(iseg) = iregion;
% %         numRegionSegments(iregion) = numRegionSegments(iregion)+1;
% %         regionSegments(iregion,numRegionSegments(iregion)) = iseg; 
%        
%         iregion
%         iseg
%         
%         iseg = 2
%     end
