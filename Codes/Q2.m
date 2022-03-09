clear;

% Set experiment parameters
L = 5;
c = 1;
dx = 0.1;
dt = 0.01;
nx = 2*L/dx;
nt = 4000;
sigma = c*dt/dx;

% Initialize the field
t = 0:dt:(nt-1)*dt;
x = -L:dx:L-dx;
q = zeros(nt,nx);
q(1,:) = exp(-x.^2);


% Integrating the equation

% For the first time step, use the forward stepping
q_trend = -c*Center_Dif(q(1,:), dx);
q(2,:) = q(1,:) + dt * q_trend;

% Using leap frog for the rest of the time steps
for i=3:nt
    q_trend = -c*Center_Dif(q(i-1,:), dx);
    q(i,:) = LeapFrog(q(i-2,:), q_trend, dt);
end


% Calculating the analytic solution q=q0(x-ct)
qa = zeros(nt,nx);
for i=1:nt
    for j=1:nx
        tnow = (i-1) * dt;
        while x(j)-c*tnow < -L
            tnow = tnow - 2*L/c;
        end
        qa(i,j) = exp(-(x(j)-c*tnow)^2);
    end
end


% Plotting the results

Delaytime = 0.5;
LineWid = 2;
FontSizeF = 18;
FontLabelF = 18;

for i=1:nt/50:nt
    % Plot each timestep
    plot(x,q(i,:),'-k','LineWidth',LineWid);
    hold on;
    plot(x,qa(i,:),'--r','LineWidth',LineWid);
    hold off;
    axis([-L L 0 1]);
    set(gca,'XTick',-L:L/5:L);
    set(gca,'YTick',0:0.2:1);
    xlabel('x');
    ylabel('q');
    formatSpec = 'sigma = %.3f';
    str = sprintf(formatSpec, sigma);
    title(str);
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
    
    formatSpec = 'Q2Sigma=%.3f.gif';
    str = sprintf(formatSpec, sigma);
    if i == 1
        imwrite(I,map,str,'gif','Loopcount',inf,'DelayTime',Delaytime);
    else
        imwrite(I,map,str,'gif','WriteMode','append','DelayTime',Delaytime);
    end
end










