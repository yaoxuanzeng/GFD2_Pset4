clear;

% Set experiment parameters
Lx = 1;
nx = 33;
Ly = Lx;
ny = nx;
delta = Lx/(nx-1);
x = 0:delta:Lx;
y = 0:delta:Ly;
dt = 0.5;
nt = 250;
Thred = 1e-6;

% Set initial conditions
x1 = 0.45;
y1 = 0.5;
x2 = 0.55;
y2 = 0.5;
eps = 1e-2;
kappa = 1e-5;
zeta = zeros(nx,ny,nt);
Psi = zeros(nx,ny,nt);
u = zeros(nx,ny,nt);
v = zeros(nx,ny,nt);


for i=1:nx
    for j=1:ny
        zeta(i,j,1) = exp(-((x(i)-x1)^2+(y(j)-y1)^2)/eps) - exp(-((x(i)-x2)^2+(y(j)-y2)^2)/eps);
    end
end
Psi(:,:,1) = zeta_Psi_xyPeriod(zeta(:,:,1), Psi(:,:,1)*0, delta, Thred);
[u(:,:,1), v(:,:,1)] = finduv_xyPeriod(Psi(:,:,1), delta);

% Calculating the first time step
Fa = Advection(u(:,:,1), v(:,:,1), zeta(:,:,1), delta);
Fd = Diffusion(zeta(:,:,1), kappa, delta);
zeta_trend = Fa + Fd;
zeta(:,:,2) = zeta(:,:,1) + zeta_trend * dt;
Psi(:,:,2) = zeta_Psi_xyPeriod(zeta(:,:,2), Psi(:,:,1), delta, Thred);
[u(:,:,2), v(:,:,2)] = finduv_xyPeriod(Psi(:,:,2), delta);

% Using leap frog for the rest of the time steps
for k=3:nt
    Fa = Advection(u(:,:,k-1), v(:,:,k-1), zeta(:,:,k-1), delta);
    Fd = Diffusion(zeta(:,:,k-1), kappa, delta);
    zeta_trend = Fa + Fd;
    zeta(:,:,k) = LeapFrog(zeta(:,:,k-2), zeta_trend, dt);
    Psi(:,:,k) = zeta_Psi_xyPeriod(zeta(:,:,k), Psi(:,:,k-1), delta, Thred);
    [u(:,:,k), v(:,:,k)] = finduv_xyPeriod(Psi(:,:,k), delta);
end

% Check the integrated zeta over the domain
zeta_sum = squeeze(sum(squeeze(sum(zeta,1)),1))/(nx*ny);

% Plotting the results

Delaytime = 0.5;
LineWid = 2;
FontSizeF = 18;
FontLabelF = 18;
tim = (1:nt)*dt;

for i=1:2:nt
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
    caxis([-0.5 0.5]);
    
    subplot(3,2,2);
    contourf(x,y,Psi(:,:,i)');
    xlabel('x');
    ylabel('y');
    colo = colorbar;
    set(get(colo,'title'),'string','\Psi')
    lim=caxis;
    ax=gca;
    ax.FontSize =FontSizeF;
    caxis([-2e-3 2e-3]);
    
    subplot(3,2,4);
    contourf(x,y,v(:,:,i)');
    xlabel('x');
    ylabel('y');
    colo = colorbar;
    set(get(colo,'title'),'string','v')
    lim=caxis;
    ax=gca;
    ax.FontSize =FontSizeF;
    caxis([-0.02 0.02]);
    
    subplot(3,2,3);
    contourf(x,y,u(:,:,i)');
    xlabel('x');
    ylabel('y');
    colo = colorbar;
    set(get(colo,'title'),'string','u')
    lim=caxis;
    ax=gca;
    ax.FontSize =FontSizeF;
    caxis([-0.01 0.01]);
    
    subplot(3,1,3);
    plot(tim(1:i),zeta_sum(1:i),'-b','LineWidth',LineWid);
    xlabel('Time');
    ylabel('[\zeta]');
    axis([tim(1) tim(end) min(zeta_sum) max(zeta_sum)]);
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
    
    formatSpec = 'Q4Akappa=%.3d.gif';
    str = sprintf(formatSpec, kappa);
    if i == 1
        imwrite(I,map,str,'gif','Loopcount',inf,'DelayTime',Delaytime);
    else
        imwrite(I,map,str,'gif','WriteMode','append','DelayTime',Delaytime);
    end
end





