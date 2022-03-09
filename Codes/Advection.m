function Advection = Advection(u, v, zeta, delta)
    % This function calculates the advection terms. u and v: zonal and
    % meridional velocity field. zeta: vorticity field. delta: grid
    % spacing.
    nx = size(zeta, 1);
    ny = size(zeta, 2);
    Fu = zeros(nx,ny);
    Fv = zeros(nx,ny);
    uzeta = u.*zeta;
    vzeta = v.*zeta;
    Fu(2:nx-1,:) = uzeta(3:nx,:) - uzeta(1:nx-2,:);
    Fv(:,2:ny-1) = vzeta(:,3:ny) - vzeta(:,1:ny-2);
    % Applying periodic boundary condition
    Fu(1,:) = uzeta(2,:) - uzeta(nx,:);
    Fu(nx,:) = uzeta(1,:) - uzeta(nx-1,:);
    Fv(:,1) = vzeta(:,2) - vzeta(:,ny);
    Fv(:,ny) = vzeta(:,1) - vzeta(:,ny-1);
    Advection = -(Fu+Fv)/(2*delta);    
 end