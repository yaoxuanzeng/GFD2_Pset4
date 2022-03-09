function Diffusion = Diffusion(zeta, kappa, delta)
    % This function calculates the diffusion term. zeta: vorticity field.
    % kappa: diffusivity. delta: grid spacing.
    nx = size(zeta, 1);
    ny = size(zeta, 2);
    Diffusionx = zeros(nx,ny);
    Diffusiony = zeros(nx,ny);
    Diffusionx(2:nx-1,:) = zeta(3:nx,:) + zeta(1:nx-2,:) - 2*zeta(2:nx-1,:);
    Diffusiony(:,2:ny-1) = zeta(:,3:ny) + zeta(:,1:ny-2) - 2*zeta(:,2:ny-1);
    % Applying periodic boundary conditions.
    Diffusionx(1,:) = zeta(2,:) + zeta(nx,:) - 2*zeta(1,:);
    Diffusionx(nx,:) = zeta(1,:) + zeta(nx-1,:) - 2*zeta(nx,:);
    Diffusiony(:,1) = zeta(:,2) + zeta(:,ny) - 2*zeta(:,1);
    Diffusiony(:,ny) = zeta(:,1) + zeta(:,ny-1) - 2*zeta(:,ny);
    Diffusion = kappa*(Diffusionx + Diffusiony)/(delta^2);
end