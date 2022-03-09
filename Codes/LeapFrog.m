function LeapFrog = LeapFrog(x, x_trend, dt)
    % Leap frog time stepping
    % x: value at t-2*dt; x_trend: tendency at t-dt; dt: timestep
    LeapFrog = x + 2*dt*x_trend;
end