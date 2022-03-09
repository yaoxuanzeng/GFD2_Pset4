clear;

% Set experiment parameters
Lx = 1;
nx = 32;
Ly = Lx;
ny = nx+1;
delta = Lx/nx;
x = 0:delta:(Lx-delta);
y = 0:delta:Ly;
dt = 1e-3;
nt = 1000;
Thred = 1e-6;

% Set initial conditions
U0 = 10;
alpha = 0.5;
n = 4;
kappa = 1e-5;
zeta = zeros(nx,ny,nt);
Psi = zeros(nx,ny,nt);
u = zeros(nx,ny,nt);
v = zeros(nx,ny,nt);

for i=1:nx
    for j=1:ny
        Psi(i,j,1) = -U0*y(j) + alpha*sin(n*pi*x(i)/Lx)*sin(pi*y(j)/Ly);
        u(i,j,1) = U0 + alpha*pi/Ly*sin(n*pi*x(i)/Lx)*cos(pi*y(j)/Ly);
        v(i,j,1) = alpha*n*pi/Lx*cos(n*pi*x(i)/Lx)*sin(pi*y(j)/Ly);
        zeta(i,j,1) = -alpha*pi^2*(n^2/Lx^2+1/Ly^2)*sin(n*pi*x(i)/Lx)*sin(pi*y(j)/Ly);
    end
end
% Set the y boundary conditions
Psi(:,1,1) = 0;
Psi(:,ny,1) = -U0*Ly;


% Calculating the first time step
Fa = AdvectionB(u(:,:,1), v(:,:,1), zeta(:,:,1), delta);
Fd = DiffusionB(zeta(:,:,1), kappa, delta);
zeta_trend = Fa + Fd;
zeta(:,:,2) = zeta(:,:,1) + zeta_trend * dt;
Psi(:,:,2) = zeta_Psi_xPeriod(zeta(:,:,2), Psi(:,:,1), delta, Thred);
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

% Check the integrated zeta over the domain
zeta_sum = squeeze(sum(squeeze(sum(zeta,1)),1))/(nx*ny);

% Plotting the results

Delaytime = 0.5;
LineWid = 2;
FontSizeF = 18;
FontLabelF = 18;
tim = (0:nt-1)*dt;

for i=1:nt/100:nt
    colormap(redblue);
    % Plot each timestep
    subplot(3,2,1);
    contourf(x,y,zeta(:,:,i)');
    xlabel('x');
    ylabel('y');
    colo = colorbar;
    set(get(colo,'title'),'string','\zeta')
    lim=caxis;
    ax=gca;
    ax.FontSize =FontSizeF;
    caxis([-100 100]);
    
    subplot(3,2,2);
    contourf(x,y,Psi(:,:,i)');
    xlabel('x');
    ylabel('y');
    colo = colorbar;
    set(get(colo,'title'),'string','\Psi')
    lim=caxis;
    ax=gca;
    ax.FontSize =FontSizeF;
    caxis([-10 0]);
    
    subplot(3,2,3);
    contourf(x,y,u(:,:,i)');
    xlabel('x');
    ylabel('y');
    colo = colorbar;
    set(get(colo,'title'),'string','u')
    lim=caxis;
    ax=gca;
    ax.FontSize =FontSizeF;
    caxis([5 15]);
    
    subplot(3,2,4);
    contourf(x,y,v(:,:,i)');
    xlabel('x');
    ylabel('y');
    colo = colorbar;
    set(get(colo,'title'),'string','v')
    lim=caxis;
    ax=gca;
    ax.FontSize =FontSizeF;
    caxis([-5 5]);
    
    subplot(3,2,4);
    contourf(x,y,v(:,:,i)');
    xlabel('x');
    ylabel('y');
    colo = colorbar;
    set(get(colo,'title'),'string','v')
    lim=caxis;
    ax=gca;
    ax.FontSize =FontSizeF;
    caxis([-5 5]);
    
    subplot(3,1,3);
    plot(tim(1:i),zeta_sum(1:i),'-b','LineWidth',LineWid);
    xlabel('Time');
    ylabel('[\zeta]');
    axis([tim(1) tim(nt) min(zeta_sum) max(zeta_sum)]);
    set(get(gca,'XLabel'),'FontSize',FontLabelF);
    set(get(gca,'YLabel'),'FontSize',FontLabelF);
    set(get(gca,'Title'),'FontSize',FontLabelF);
    ax=gca;
    ax.FontSize =FontSizeF;
    
    % Making animation
    drawnow;
    F = getframe(gcf);
    I = frame2im(F);
    [I,map] = rgb2ind(I,256);
    
    formatSpec = 'Q4Bkappa=%.3d.gif';
    str = sprintf(formatSpec, kappa);
    if i == 1
        imwrite(I,map,str,'gif','Loopcount',inf,'DelayTime',Delaytime);
    else
        imwrite(I,map,str,'gif','WriteMode','append','DelayTime',Delaytime);
    end
end





