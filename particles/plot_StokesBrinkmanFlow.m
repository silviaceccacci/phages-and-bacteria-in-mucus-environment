function plot_StokesBrinkmanFlow(U_interp, x_min, x_max, y_min, y_max, micron, Omega_X, Omega_Y)
    % plotStokesBrinkmanFlow - Plot Stokes-Brinkman flow field magnitude
    %
    % Inputs:
    %   U_interp : function handle, takes [x, y] and returns [u, v]
    %   x_min, x_max, y_min, y_max : domain bounds
    %   micron   : scaling factor (m -> μm)
    %   Omega_X, Omega_Y : domain dimensions in m

    nplot = 100;
    hplot = (x_max - x_min) / nplot;
    [XX, YY] = meshgrid(x_min:hplot:x_max, y_min:hplot:y_max);

    UU = zeros(size(XX,1), size(XX,2), 2);
    for i = 1:size(XX,1)
        for j = 1:size(XX,2)
            UU(i,j,:) = U_interp([XX(i,j), YY(i,j)]);
        end
    end
    UU = sqrt(sum(UU .* UU, 3));

    % Convert to microns
    XX = XX ./ micron;
    YY = YY ./ micron;

    % Plot
    figure; hold on;
    fig = surf(XX, YY, UU);
    set(fig, 'edgecolor', 'none'); view(2)
    cb = colorbar();
    ylabel(cb, '$u \ (\mu m / s)$', 'Interpreter', 'LaTeX', 'FontSize', 16)

    xlabel('$\Omega_x \ (\mu m)$', 'Interpreter', 'LaTeX', 'FontSize', 16);
    ylabel('$\Omega_y \ (\mu m)$', 'Interpreter', 'LaTeX', 'FontSize', 16);
    title('Stokes-Brinkman flow', 'Interpreter', 'LaTeX', 'FontSize', 16);
    xlim([0 Omega_X / micron]);
    ylim([0 Omega_Y / micron]);
    set(gca, 'FontSize', 16);
end
