function setPlotMarkerSizeInDataUnits(hPlot, desiredDiameterInDataUnits)
    % hPlot: handle to the plot (line) object
    % desiredDiameterInDataUnits: marker diameter in axis units (linear)

    ax = ancestor(hPlot, 'axes');

    % Save and temporarily set axes units to pixels
    originalUnits = get(ax, 'Units');
    set(ax, 'Units', 'pixels');
    axPos = get(ax, 'Position');  % [left bottom width height]
    set(ax, 'Units', originalUnits);  % restore

    % Get axis limits and ranges
    xLimits = xlim(ax)
    yLimits = ylim(ax)
    xRange = diff(xLimits);
    yRange = diff(yLimits);

    % Pixels per data unit
    xScale = axPos(3) / xRange;  % pixels per x data unit
    yScale = axPos(4) / yRange;  % pixels per y data unit

    % Use average scale for isotropic marker size
    avgScale = mean([xScale, yScale]);

    % Convert diameter in data units → diameter in pixels → points
    markerDiameterPixels = desiredDiameterInDataUnits * avgScale;
    markerDiameterPoints = markerDiameterPixels * 72 / get(gcf, 'ScreenPixelsPerInch');  % px → pt

    % Set the plot marker size (in points)
    set(hPlot, 'MarkerSize', markerDiameterPoints);
end
