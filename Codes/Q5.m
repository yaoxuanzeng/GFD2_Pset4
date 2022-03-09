clear;

% Set experiment parameters
nx = 34;
ny = 17;
delta = 312500;
x = 0:delta:(nx-1)*delta;
y = -(ny-1)/2*delta:delta:(ny-1)/2*delta;
dt = 900;
nt = 96*5+1;
Thred = 1e-6;

% Set initial conditions
kappa = 1e2;
zeta = zeros(nx,ny,nt);
zeta0x = zeros(nx,ny);
zeta0y = zeros(nx,ny);
Psi = zeros(nx,ny,nt);
u = zeros(nx,ny,nt);
v = zeros(nx,ny,nt);

Psi(:,:,1) = ReadPsi();
[u(:,:,1), v(:,:,1)] = finduv_Q5(Psi(:,:,1), delta);
zeta0x(2:nx-1,:) = Psi(3:nx,:,1) + Psi(1:nx-2,:,1) - 2*Psi(2:nx-1,:,1);
zeta0x(1,:) = Psi(2,:,1) + Psi(nx,:,1) - 2*Psi(1,:,1);
zeta0x(nx,:) = Psi(1,:,1) + Psi(nx-1,:,1) - 2*Psi(nx,:,1);
zeta0y(:,2:ny-1) = Psi(:,3:ny,1) + Psi(:,1:ny-2,1) - 2*Psi(:,2:ny-1,1);
zeta(:,:,1) = (zeta0x + zeta0y)/delta^2;

% Calculating the first time step
Fa = AdvectionQ5(u(:,:,1), v(:,:,1), zeta(:,:,1), delta);
Fd = DiffusionQ5(zeta(:,:,1), kappa, delta);
zeta_trend = Fa + Fd;
zeta(:,:,2) = zeta(:,:,1) + zeta_trend * dt;
Psi(:,:,2) = zeta_Psi_Q5(zeta(:,:,2), Psi(:,:,1), delta, Thred);
[u(:,:,2), v(:,:,2)] = finduv_xPeriod(Psi(:,:,2), delta);
zeta(:,1,2) = 0;
zeta(:,ny,2) = 0;

% Using leap frog for the rest of the time steps
for k=3:nt
    Fa = AdvectionB(u(:,:,k-1), v(:,:,k-1), zeta(:,:,k-1), delta);
    Fd = DiffusionB(zeta(:,:,k-1), kappa, delta);
    zeta_trend = Fa + Fd;
    zeta(:,:,k) = LeapFrog(zeta(:,:,k-2), zeta_trend, dt);
    Psi(:,:,k) = zeta_Psi_xPeriod(zeta(:,:,k), Psi(:,:,k-1), delta, Thred);
    [u(:,:,k), v(:,:,k)] = finduv_xPeriod(Psi(:,:,k), delta);
    zeta(:,1,k) = 0;
    zeta(:,ny,k) = 0;
end


% Plotting the results

Delaytime = 0.5;
FontSizeF = 14;
FontLabelF = 14;

x = x/1e6;
y = y/1e6;
zeta = zeta*1e4;
Psi = Psi/1e8;

for i=1:4*2:nt
    formatSpec1 = 'T = %.3d hrs';
    str1 = sprintf(formatSpec1, (i-1)/4);

    colormap(redblue);
    % Plot each timestep
    subplot(2,2,1);
    contourf(x,y,zeta(:,:,i)');
    xlabel('x (10^3 km)');
    ylabel('y (10^3 km)');
    title(str1);
    colo = colorbar;
    set(get(colo,'title'),'string','\zeta (10^-^4 s^-^1)')
    lim=caxis;
    ax=gca;
    ax.FontSize =FontSizeF;
    caxis([-2 2]);
    
    subplot(2,2,2);
    contourf(x,y,Psi(:,:,i)');
    xlabel('x (10^3 km)');
    ylabel('y (10^3 km)');
    title(str1);
    colo = colorbar;
    set(get(colo,'title'),'string','\Psi (10^8 m^2/s)')
    lim=caxis;
    ax=gca;
    ax.FontSize =FontSizeF;
    caxis([4 6]);
    
    subplot(2,2,4);
    contourf(x,y,v(:,:,i)');
    xlabel('x (10^3 km)');
    ylabel('y (10^3 km)');
    title(str1);
    colo = colorbar;
    set(get(colo,'title'),'string','v (m/s)')
    lim=caxis;
    ax=gca;
    ax.FontSize =FontSizeF;
    caxis([-20 20]);
    
    subplot(2,2,3);
    contourf(x,y,u(:,:,i)');
    xlabel('x (10^3 km)');
    ylabel('y (10^3 km)');
    title(str1);
    colo = colorbar;
    set(get(colo,'title'),'string','u (m/s)')
    lim=caxis;
    ax=gca;
    ax.FontSize =FontSizeF;
    caxis([-50 50]);
    
    % Making animation
    drawnow;
    F = getframe(gcf);
    I = frame2im(F);
    [I,map] = rgb2ind(I,256);
    
    formatSpec = 'Q5kappa=%.3d.gif';
    str = sprintf(formatSpec, kappa);
    if i == 1
        imwrite(I,map,str,'gif','Loopcount',inf,'DelayTime',Delaytime);
    else
        imwrite(I,map,str,'gif','WriteMode','append','DelayTime',Delaytime);
    end
end




