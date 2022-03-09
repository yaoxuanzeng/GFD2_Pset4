function zeta_Psi_xPeriod = zeta_Psi_xPeriod(zeta, Psi_old, delta, Thred)
    % This function solves the Poisson equation to get Psi from zeta. zeta:
    % vorticity field; Psi_old: streamfunction field of the last time step,
    % since "guessing" from the last time step should be closer to the
    % solution; delta: grid spacing; Thred: threshold for the accuracy.
    
    nx = size(zeta,1);
    ny = size(zeta,2);
    
    % Set the boundary condition
    Psi_new = Psi_old;
    
    % Calculate the first step, and the difference
    for j=2:ny-1
        for i=2:nx-1
            Psi_new(i,j) = (Psi_old(i+1,j)+Psi_old(i-1,j)+Psi_old(i,j+1)+Psi_old(i,j-1)-delta^2*zeta(i,j))/4;
        end
        % Applying periodic boundary condition in x. Since Psi is fixed at
        % the y boundary, there's no need to update Psi there
        Psi_new(1,j) = (Psi_old(2,j)+Psi_old(nx,j)+Psi_old(1,j+1)+Psi_old(1,j-1)-delta^2*zeta(1,j))/4;
        Psi_new(nx,j) = (Psi_old(1,j)+Psi_old(nx-1,j)+Psi_old(nx,j+1)+Psi_old(nx,j-1)-delta^2*zeta(nx,j))/4;
    end
    
    Psi_difno = abs(Psi_new-Psi_old);
    Diff_newold = sum(sum(Psi_difno))/(nx*ny);
    
    % Calculate steps until difference is smaller than threshold
    while Diff_newold > Thred
        Psi_old = Psi_new;
        for j=2:ny-1
           for i=2:nx-1
               Psi_new(i,j) = (Psi_old(i+1,j)+Psi_old(i-1,j)+Psi_old(i,j+1)+Psi_old(i,j-1)-delta^2*zeta(i,j))/4;
           end
           Psi_new(1,j) = (Psi_old(2,j)+Psi_old(nx,j)+Psi_old(1,j+1)+Psi_old(1,j-1)-delta^2*zeta(1,j))/4;
           Psi_new(nx,j) = (Psi_old(1,j)+Psi_old(nx-1,j)+Psi_old(nx,j+1)+Psi_old(nx,j-1)-delta^2*zeta(nx,j))/4;
        end
        Psi_difno = abs(Psi_new-Psi_old);
        Diff_newold = sum(sum(Psi_difno))/(nx*ny);
    end
    zeta_Psi_xPeriod = Psi_new;
end