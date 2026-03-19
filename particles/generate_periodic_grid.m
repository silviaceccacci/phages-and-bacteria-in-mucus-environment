function positions = generate_periodic_grid(nx, ny, Lx, Ly)

    dx = Lx / nx;
    dy = Ly / ny;

    x = dx/2 : dx : Lx - dx/2;
    y = dy/2 : dy : Ly - dy/2;

    [X, Y] = meshgrid(x, y);
    positions = [X(:), Y(:)];
end
