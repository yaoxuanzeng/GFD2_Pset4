function Center_Dif = Center_Dif(x, dx)
    % Center differences calculation
    % x: target array; nx: array length; dx: interval
    nx = size(x,2);
    Center_Dif = zeros(1,nx);
    Center_Dif(2:nx-1) = (x(3:nx) - x(1:nx-2))/(2*dx);
    Center_Dif(1) = (x(2)-x(nx))/(2*dx);
    Center_Dif(nx) = (x(1)-x(nx-1))/(2*dx);
end