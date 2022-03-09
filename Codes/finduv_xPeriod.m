function [u, v] = finduv_xPeriod(Psi, delta)
    % This function calculates velocity field from the stream function.
    % Psi: stream function. delta: grid spacing
    nx = size(Psi, 1);
    ny = size(Psi, 2);
    u = zeros(nx,ny);
    v = zeros(nx,ny);
    v(2:nx-1,:) = (Psi(3:nx,:)-Psi(1:nx-2,:))/(2*delta);
    u(:,2:ny-1) = -(Psi(:,3:ny)-Psi(:,1:ny-2))/(2*delta);
    % Applying periodic boundary condition in x direction. For y boundary
    % condition, the gradient is calculated by the first and second grid
    % points instead of center difference method.
    v(1,:) = (Psi(2,:)-Psi(nx,:))/(2*delta);
    v(nx,:) = (Psi(1,:)-Psi(nx-1,:))/(2*delta);
    u(:,1) = -(Psi(:,2)-Psi(:,1))/(delta);
    u(:,ny) = -(Psi(:,ny)-Psi(:,ny-1))/(delta);
end