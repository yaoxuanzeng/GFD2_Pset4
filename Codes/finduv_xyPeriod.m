function [u, v] = finduv_xyPeriod(Psi, delta)
    % This function calculates velocity field from the stream function.
    % Psi: stream function. delta: grid spacing
    nx = size(Psi, 1);
    ny = size(Psi, 2);
    u = zeros(nx,ny);
    v = zeros(nx,ny);
    u(:,2:ny-1) = -(Psi(:,3:ny)-Psi(:,1:ny-2))/(2*delta);
    v(2:nx-1,:) = (Psi(3:nx,:)-Psi(1:nx-2,:))/(2*delta);
    % Applying periodic boundary condition
    u(:,1) = -(Psi(:,2)-Psi(:,ny))/(2*delta);
    u(:,ny) = -(Psi(:,1)-Psi(:,ny-1))/(2*delta);
    v(1,:) = (Psi(2,:)-Psi(nx,:))/(2*delta);
    v(nx,:) = (Psi(1,:)-Psi(nx-1,:))/(2*delta);
end