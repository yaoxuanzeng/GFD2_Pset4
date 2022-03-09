function DiffusionQ5 = DiffusionQ5(zeta, kappa, delta)
    % This function calculates the diffusion term. zeta: vorticity field.
    % kappa: diffusivity. delta: grid spacing
    nx = size(zeta, 1);
    ny = size(zeta, 2);
    DiffusionQ5 = zeros(nx,ny);
    for j=2:ny-1
        for i=2:nx-1
            DiffusionQ5(i,j) = zeta(i+1,j) + zeta(i-1,j) + zeta(i,j+1) + zeta(i,j-1) - 4*zeta(i,j);
        end
        % Applying periodic boundary conditions in x direction. Since
        % zeta=0 at y boundaries, there's no need to calculate the
        % diffusion terms there.
        DiffusionQ5(1,j) = (zeta(2,j) + zeta(nx,j) + zeta(1,j+1) + zeta(1,j-1) - 4*zeta(1,j));
        DiffusionQ5(nx,j) = (zeta(1,j) + zeta(nx-1,j) + zeta(nx,j+1) + zeta(nx,j-1) - 4*zeta(nx,j))/4;
    end
    DiffusionQ5 = kappa*DiffusionQ5/(delta^2);
end