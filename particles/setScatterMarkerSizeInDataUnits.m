function setScatterMarkerSizeInDataUnits(hScatter, desiredDiameterInDataUnits)
    % hScatter: handle to the scatter plot
    % desiredDiameterInDataUnits: marker diameter in axis units (assumes square aspect)

    ax = ancestor(hScatter, 'axes');

    % Save and temporarily set axes units to pixels
    originalUnits = get(ax, 'Units');
    set(ax, 'Units', 'pixels');
    axPos = get(ax, 'Position');  % [left bottom width height]
    set(ax, 'Units', originalUnits);  % restore

    % Get axis limits and ranges
    xLimits = xlim(ax);
    yLimits = ylim(ax);
    xRange = diff(xLimits);
    yRange = diff(yLimits);

    % Pixels per data unit
    xScale = axPos(3) / xRange;  % pixels per x data unit
    yScale = axPos(4) / yRange;  % pixels per y data unit

    % Use average scale for isotropic marker size
    avgScale = mean([xScale, yScale]);

    % Convert diameter in data units → diameter in pixels → points → points^2
    markerDiameterPixels = desiredDiameterInDataUnits * avgScale;
    markerDiameterPoints = markerDiameterPixels * 72 / get(gcf, 'ScreenPixelsPerInch');  % px → pt
    markerSizePointsSquared = markerDiameterPoints^2;

    % Set the scatter size
    hScatter.SizeData = markerSizePointsSquared;
end
