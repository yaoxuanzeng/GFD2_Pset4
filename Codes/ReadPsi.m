function ReadPsi = ReadPsi()
    % This function reads the initial geopotential data and converts it to
    % the stream function as Psi = g*Phi/f_0.
    Phi = importdata('./msc11.txt');
    f0 = 2*(2*pi/86400)*sin(pi/4);
    g = 9.8;
    ReadPsi = g*Phi'/f0;
end