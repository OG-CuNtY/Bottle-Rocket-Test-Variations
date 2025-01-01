%{
Author: Mark Anifowose 
Assignment: Water Bottle Rocket Flight.
Creation Date: 11/15/2024
Inputs: project2verification.mat
Outputs: A model of the flight of a bottle rocket in comparison to a verification test case. 
              A flight trajectory (z vs. x position) and thrust profile (F vs. t).
Purpose: The goal of this project is to model the trajectory of the bottle rocket launch, 
              using the numerical integration of a system of ordinary differential equation.
%}

%% Initialize test/verification data
clc; clear; close all; load project2verification.mat;

%% Get Constants
const = getConst();

%% Set Height, Distance, Thrust, and Time data for verification
heigthVerify = verification.height();
distanceVerify = verification.distance();
thrustVerify = verification.thrust();
timeVerify = verification.time();

%% Initial State Vector
statevector_0 = [const.x0, 0, const.z0, 0, const.mr0, const.Vair, const.Mair];

%% Time Span
tspan = [0 5];

%% Solve Combined Phases with ode45
[t, state] = ode45(@(t, statevector) rocket_phases(t, statevector, const), tspan, statevector_0);


%% Calculate Thrust Profile
Thrust = zeros(length(t), 1);
for i = 1:length(t)
    [~, Thrust(i)] = rocket_phases(t(i), state(i, :), const);
end

%% Plot Results

% Save Values
maxThrust = max(Thrust);
maxHeight = max(state(:, 3));
maxDistance = max(state(:, 1));

% Display Results
fprintf("Max Height (current data): %2.1f ", maxHeight);
fprintf("\nMax Distance (current data): %2.1f ", maxDistance);
fprintf("\nPeak Thrust (current data): %2.1f ", maxThrust);


% Trajectory Plot
figure;
plot(distanceVerify, heigthVerify, 'b-', 'LineWidth', 4)
hold on;
plot(state(:, 1), state(:, 3), 'r-', 'LineWidth', 1.5);
xlabel('Distance (m)');
ylabel('Height (m)');
title('Rocket Trajectory');
xlim([0 max(state(:, 1))+5]);
ylim([0 20]);
legend("Verification Data", "Current Data", Location="best");
grid on;
grid minor;

% Thrust vs Time Plot
figure;
plot(timeVerify, thrustVerify, 'b-', 'LineWidth', 4);
hold on;
plot(t, Thrust, 'r-', 'LineWidth', 1.5); 
xlabel('Time (s)');
ylabel('Thrust (N)');
title('Thrust vs Time');
xticks(0:0.05:0.3);
xlim([0 0.3]);
ylim([0 max(Thrust)]);
legend("Verification Data", "Current Data", Location="best");
grid on;
grid minor;

%% Rocket Phases Function
function [d_statevector_dt, F_thrust] = rocket_phases(t, statevector, const)
    % Extract state variables
    x = statevector(1); % x-position
    vx = statevector(2); % x-velocity
    z = statevector(3); % z-position
    vz = statevector(4); % z-velocity
    mr = statevector(5); % Mass of rocket
    Vair = statevector(6); % Volume of air
    mAir = statevector(7); % Mass of air
    
    vel = [vx; vz];
    velMag = norm(vel);

    % Determine the heading vector
    if norm([x, z - const.z0]) > const.ls
        h_hat = vel/velMag;
    else
        h_hat = [cosd(const.theta); sind(const.theta)];
    end

    % Compute Drag Force
    D = 0.5 * const.rho_air * (velMag^2) * const.CD * const.A_body;
    
    % Set Pressure for Phase change
    if ~Vair < const.Vb
        pEnd = const.p0 * ((const.Vair / const.Vb)^const.gamma);
        p = pEnd * ((mAir / const.Mair)^const.gamma);
    end


    % Phase determination
    if Vair < const.Vb % Phase 1: Water exhaustion

        p = const.p0 * ((const.Vair / Vair)^const.gamma);

        F_thrust = 2*const.cdis*const.A_throat*(p - const.pa);

        mDot_r = -1*const.cdis*const.A_throat*sqrt(2*const.rho_w*(p - const.pa)); % change in rocket mass

        vDot_Air = const.cdis * const.A_throat *sqrt((2/const.rho_w)*(const.p0*((const.Vair/Vair)^const.gamma)-const.pa)); % change in air volume

        % Net forces
        Fx = (F_thrust * h_hat(1)) - (D * h_hat(1)); 
        Fz = (F_thrust * h_hat(2)) - (D * h_hat(2)) - (mr * const.g);

        % State derivatives
        ax = Fx / mr;
        az = Fz / mr;
        

        % Return derivatives
        d_statevector_dt = [vx; ax; vz; az; mDot_r; vDot_Air; 0];

    elseif p > const.pa % Phase 2: Air exhaustion

        rhoODE = mAir / const.Vb;

        Temp = p / (rhoODE * const.Rair);

        pCritical = p*(2/(const.gamma+1))^(const.gamma/(const.gamma-1)); % Critical Pressure 

        if pCritical > const.pa

            pExit = pCritical;

            TempExit = (2 / (const.gamma + 1)) * Temp;

            vExit = sqrt(const.gamma * const.Rair * TempExit);

            rhoExit = pExit / (const.Rair * TempExit);

        elseif pCritical < const.pa
            pExit = const.pa;

            MachExit = sqrt((2/(const.gamma-1))*(((p/const.pa)^((const.gamma-1)/const.gamma))-1));

            TempExit = Temp/(1+((const.gamma-1)/2)*MachExit^2);
            
            rhoExit = pExit / (const.Rair * TempExit);

            vExit = MachExit * sqrt(const.gamma * const.Rair * TempExit);
        end

        mDot_air = -1*const.cdis * rhoExit * const.A_throat * vExit; % change in air mass

        F_thrust = -1*mDot_air * vExit + (pExit - const.pa) * const.A_throat;

        % Net forces
        Fx = (F_thrust * h_hat(1)) - (D * h_hat(1));
        Fz = (F_thrust * h_hat(2)) - (D * h_hat(2)) - (mr * const.g);

        % Accelerations
        ax = Fx / mr;
        az = Fz / mr;

        % Return derivatives
        d_statevector_dt = [vx; ax; vz; az; mDot_air; 0; mDot_air];

    else % Phase 3: Ballistic flight
        % Thrust is zero
        F_thrust = 0;

        % Net forces
        Fx = -D * h_hat(1);
        Fz = -D * h_hat(2) - (mr * const.g);

        % Accelerations
        ax = Fx / mr;
        az = Fz / mr;

        if z < 0 && x > 0 
            vx = 0;
            vz = 0;
        end

        % Return derivatives
        d_statevector_dt = [vx; ax; vz; az; 0; 0; 0];
    end
end

%% Constants Function
function const = getConst()
    const.g = 9.81; % Gravity (m/s^2)
    const.cdis = 0.78; % Discharge coefficient
    const.rho_air = 0.961; % Air density (kg/m^3)
    const.Vb = 0.002; % Bottle volume (m^3)
    const.pa = 83426.563088; % Atmospheric pressure (Pa)
    const.gamma = 1.4; % Specific heat ratio
    const.rho_w = 1000; % Water density (kg/m^3)
    const.de = 0.021; % Nozzle diameter (m)
    const.dB = 0.105; % Bottle diameter (m)
    const.Rair = 287; % Specific gas constant for air (J/(kg·K))
    const.mB = 0.15; % Mass of empty bottle (kg)
    const.CD = 0.311; % Drag coefficient                          %THIS VALUE CHANGED FROM 0.425
    const.p0 = 494390 + const.pa; % Initial air pressure (Pa)
    const.Vi_water = 0.0009; % Initial water volume (m^3)
    const.T0 = 310; % Initial air temperature (K)
    const.v0 = 0.0; % Initial velocity of rocket (m/s)
    const.theta = 25; % Launch angle (degrees)
    const.x0 = 0.0; % Initial x-position (m)
    const.z0 = 0.25; % Initial z-position (m)
    const.ls = 0.5; % Length of test stand (m)
    const.A_throat = pi * (const.de / 2)^2; % Nozzle throat area (m^2)
    const.A_body = pi * (const.dB / 2)^2; % Bottle cross-sectional area (m^2)
    const.Vair = const.Vb - const.Vi_water; % Initial volume of air
    const.Mair = (const.p0 * const.Vair) / (const.Rair * const.T0); % Initial mass of air
    const.mWater = const.rho_w * (const.Vb - const.Vair); % Initial mass of water
    const.mr0 = const.mB + (const.rho_w*(const.Vb - const.Vair)) + const.Vair*(const.p0/(const.Rair * const.T0)); % Initial mass of rocket
end
