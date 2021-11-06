% Grid parameters
L  = 1;     % Constant: Length in m
N  = 100;   % Constant: Number of cells
Lx = L;     % Grid length along the x-axis
Ly = L;     % Grid length along the y-axis
Nx = N;     % Number of cells along the x-axis
Ny = N;     % Number of cells along the x-axis
Dx = Lx/Nx; % Cell length along the x-axis
Dy = Ly/Ny; % Cell length along the y-axis

% Material parameters
K = 401;    % Thermal conductivity of copper in W/(mK)

% Boundary conditions
T_top       = 100;  % Temperature at the top side in degrees Celcius
T_bottom    = 0;    % Temperature at the bottom side in degrees Celcius
T_left      = 75;   % Temperature at the left side in degrees Celcius
T_right     = 50;   % Temperature at the right side in degrees Celcius

% Initialise the grid
grid = zeros(Ny, Nx);

% Initialise initial conditions
for i = 1:Nx
    grid(1, i)  = T_top;
    grid(Ny, i) = T_bottom;
end

for i = 1:Ny
    grid(i, 1)  = T_left;
    grid(i, Nx) = T_right;
end

% Stabilise
Err     = 1;        % Simulation error (initialize with any large enough value)
MaxErr  = 1e-6;     % Maximum simulation error
NumIter = 0;        % The number of iterations
while Err > MaxErr   % Run until the error is small enough (smaller than MaxErr)
    NumIter = NumIter + 1;
    Err = 0;
    grid_next = grid;
    for i = 2:Ny-1
        for j = 2:Nx-1
            grid_next(i, j) = (grid(i-1,j)    ...
                        + grid(i+1, j)   ...
                        + grid(i, j-1)   ...
                        + grid(i, j+1)) / 4;
            Err = Err + grid_next(i, j) - grid(i, j);
        end
    end
    grid = grid_next;
end
% Display the number of iterations needed
%fprintf("Simulation ends after %d iterations.\n", NumIter)

% Plot the temperature distribution
temp = heatmap(grid, 'Colormap', jet, 'GridVisible', 'off');
temp.XDisplayLabels = nan(size(temp.XDisplayData));
temp.YDisplayLabels = nan(size(temp.YDisplayData));

%disp(temp)
saveas(gca, 'Temperature.jpg')

% Calculate the heat flux at each point
[Tx, Ty] = gradient(grid);  % Temperature gradient in x and y direction, respectively
Fx = -K*Tx;                 % Heat flux in x direction in W/m^2
Fy = -K*Ty;                 % Heat flux in y direction in W/m^2

% Plot heat flux distribution
flux = quiver(Fx, Fy, 20);
set(gca, 'YDir','reverse')
axis equal

%disp(flux)
saveas(gca, 'Heat_Flux.jpg')