function zeta_Psi_xyPeriod = zeta_Psi_xyPeriod(zeta, Psi_old, delta, Thred)
    % This function solves the Poisson equation to get Psi from zeta. zeta:
    % vorticity field; Psi_old: streamfunction field of the last time step,
    % since "guessing" from the last time step should be closer to the
    % solution; delta: grid spacing; Thred: threshold for the accuracy.
    
    nx = size(zeta,1);
    ny = size(zeta,2);

    % Calculate the first step, and the difference
    Psix = zeros(nx,ny);
    Psiy = zeros(nx,ny);
    Psix(2:nx-1,:) = Psi_old(3:nx,:) + Psi_old(1:nx-2,:);
    Psiy(:,2:ny-1) = Psi_old(:,3:ny) + Psi_old(:,1:ny-2);
    % Applying periodic boundary condition
    Psix(1,:) = Psi_old(2,:) + Psi_old(nx,:);
    Psix(nx,:) = Psi_old(1,:) + Psi_old(nx-1,:);
    Psiy(:,1) = Psi_old(:,2) + Psi_old(:,ny);
    Psiy(:,ny) = Psi_old(:,1) + Psi_old(:,ny-1);
    Psi_new = (Psix + Psiy - zeta*delta^2)/4;
    Psi_difno = abs(Psi_new-Psi_old);
    Diff_newold = sum(sum(Psi_difno))/(nx*ny);

    % Calculate steps until difference is smaller than threshold

    while Diff_newold > Thred
        Psi_old = Psi_new;
        Psix = zeros(nx,ny);
        Psiy = zeros(nx,ny);
        Psix(2:nx-1,:) = Psi_old(3:nx,:) + Psi_old(1:nx-2,:);
        Psix(1,:) = Psi_old(2,:) + Psi_old(nx,:);
        Psix(nx,:) = Psi_old(1,:) + Psi_old(nx-1,:);
        Psiy(:,2:ny-1) = Psi_old(:,3:ny) + Psi_old(:,1:ny-2);
        Psiy(:,1) = Psi_old(:,2) + Psi_old(:,ny);
        Psiy(:,ny) = Psi_old(:,1) + Psi_old(:,ny-1);
        Psi_new = (Psix + Psiy - zeta*delta^2)/4;
        Psi_difno = abs(Psi_new-Psi_old);
        Diff_newold = sum(sum(Psi_difno))/(nx*ny);
    end
    zeta_Psi_xyPeriod = Psi_new;
end



