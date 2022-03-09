clear;

% Set experiment parameters
Lx = 1;
nx = 32+1;
Ly = 1.5*Lx;
ny = 1.5*(nx-1)+1;
delta = Lx/(nx-1);
x = 0:delta:Lx;
y = 0:delta:Ly;
Thred = 1e-6;

% Initialize the field
zeta = zeros(nx,ny);
zeta = zeta + 4; % Initialize the vorticity field
% Store the difference between time steps and with the accurate solution
Diff_oldnew = [];
Diff_newacc = [];
Psi_acc = zeros(nx,ny);
% Set the "accurate" solution
for i=1:nx
    for j=1:ny
        Psi_acc(i,j) = x(i)^2 + y(j)^2;
    end
end
% Set the boundary condition
Psi_old = zeros(nx,ny);
Psi_old(1,:) = x(1)^2 + y.^2;
Psi_old(nx,:) = x(nx)^2 + y.^2;
Psi_old(:,1) = x.^2 + y(1)^2;
Psi_old(:,ny) = x.^2 + y(ny)^2;
Psi_new = Psi_old;

% Calculate the first step, and the difference between previous and current
% steps, and between numerical solution and the accurate solution
Psi_old = Psi_new;
for i=2:nx-1
    for j=2:ny-1
        Psi_new(i,j) = (Psi_old(i+1,j)+Psi_old(i-1,j)+Psi_old(i,j+1)+Psi_old(i,j-1)-delta^2*zeta(i,j))/4;
    end
end
Psi_difno = abs(Psi_new-Psi_old);
Psi_difna = abs(Psi_new-Psi_acc);
Diff_newold(1) = sum(sum(Psi_difno))/(nx*ny);
Diff_newacc(1) = sum(sum(Psi_difna))/(nx*ny);

% Calculate steps until difference is smaller than threshold
m = 1;
while Diff_newold(m) > Thred
    Psi_old = Psi_new;
    for i=2:nx-1
        for j=2:ny-1
            Psi_new(i,j) = (Psi_old(i+1,j)+Psi_old(i-1,j)+Psi_old(i,j+1)+Psi_old(i,j-1)-delta^2*zeta(i,j))/4;
        end
    end
    m = m+1;
    Psi_difno = abs(Psi_new-Psi_old);
    Psi_difna = abs(Psi_new-Psi_acc);
    Diff_newold(m) = sum(sum(Psi_difno))/(nx*ny);
    Diff_newacc(m) = sum(sum(Psi_difna))/(nx*ny);
end

% Plot the results
LineWid = 2;
FontSizeF = 18;
FontLabelF = 18;
step = 1:m;

figure(1);
colormap(flipud(hot));
subplot(2,2,1);
plot(step,Diff_newold,'-b','LineWidth',LineWid);
xlabel('Steps');
ylabel('Difference');
title('Difference between steps');
set(get(gca,'XLabel'),'FontSize',FontLabelF);
set(get(gca,'YLabel'),'FontSize',FontLabelF);
set(get(gca,'Title'),'FontSize',FontLabelF);
ax=gca;
ax.FontSize =FontSizeF;

subplot(2,2,2);
plot(step,Diff_newacc,'-b','LineWidth',LineWid);
xlabel('Steps');
ylabel('Difference');
title('Numerical - Accurate');
set(get(gca,'XLabel'),'FontSize',FontLabelF);
set(get(gca,'YLabel'),'FontSize',FontLabelF);
set(get(gca,'Title'),'FontSize',FontLabelF);
ax=gca;
ax.FontSize =FontSizeF;

levelPsi = 0:0.2:4;
subplot(2,3,4);
contourf(x,y,Psi_new',levelPsi);
colo = colorbar;
set(get(colo,'title'),'string','\Psi')
xlabel('x');
ylabel('y');
title('Numerical solution')
lim=caxis
ax=gca;
ax.FontSize =FontSizeF;
caxis([0 3]);

subplot(2,3,5);
contourf(x,y,Psi_acc',levelPsi);
colo = colorbar;
set(get(colo,'title'),'string','\Psi')
xlabel('x');
ylabel('y');
title('Accurate solution')
lim=caxis
ax=gca;
ax.FontSize =FontSizeF;
caxis([0 3]);

subplot(2,3,6);
contourf(x,y,Psi_difna');
colo = colorbar;
set(get(colo,'title'),'string','\Psi_n-\Psi_a')
xlabel('x');
ylabel('y');
title('Numerical-Accurate')
lim=caxis
ax=gca;
ax.FontSize =FontSizeF;
%caxis([-2e-6 2e-6]);


