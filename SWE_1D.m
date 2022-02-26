%--------------------------------------------------------------------------
% NUMERICAL INTEGRATION OF 1D SWE WITH MACCORMACK SCHEME
%--------------------------------------------------------------------------

% initialize code
clear; close all; clc; %tic

%% INPUT DATA

% CONSTANTS
g  = 9.806;     % gravity acceleration (m/s2)
Cc = 0.611;     % contraction coefficient

% GEOMETRY DATA
Ks = 20;        % Strickler roughness coefficient (m1/3 s-1)
S0 = 0.001;     % longitudinal bed slope (-)
n  = 1000;      % number of nodes (-)
dx = 1;         % node spacing (m)

% BOUNDARY CONDITIONS
gate_opening   = [2, 2, 0.1, 0.1, 2];      % sluice gate inizial opening (m)
gate_times     = [0, 1, 2, 3, 4] * 3600;  % sluice gate times (s)
Hgate          = 8;                        % upstream gate head

% INITIAL CONDITIONS
qus0 = 2;       % unit discharge (m3/s/m) Q = B * qus0

% SIMULATION SETTING
TT  = 12 * 3600;  % total time (s)
dt0 = 60;         % time step
itmax = 10^8;     % maximun number of allowed iteration

%% GEOMETRY BUILDING

% local roughness
Ks = repmat(Ks,n,1);              
% Ks(n/2:end) = Ks(n/2:end)/1.5;

% local bed slope
S0 = repmat(S0,n,1);             
% S0(n/2:end) = S0(n/2:end)*2;

% longitudinal coordinate (m)
x  = (0 : dx : (n-1)*dx)'; 

% bed elevation (m s.m.m.)
z = zeros(n,1);                 
for i = n-1:-1:1
    z(i) = z(i+1) + (x(i+1)-x(i)) * ( S0(i+1) + S0(i))/2;
end

%% INITIAL CONDITIONS: UNIFORM FLOW

q = repmat(qus0,n,1);
D = uniform_flow_depth(q, Ks, S0, 'q');
H = D + z;
U = qus0 ./ D;
Fr = U ./ sqrt(g * D);
Se = friction_slope(Ks, U, D);

%% TIME INTEGRATION

% initialize time and iteration counter
t  = 0;
it = 0;

% initialize gate pointer
jgate = 1;

% loop over time
while t < TT && it < itmax;
   
% update time step
    dt = min(dt0, cfl(g, dx, U, D));
    dt = min(dt, TT-t);
    
% upstream boundary conditions - compute gate opening
    [gate_a, jgate] = gate_operations(t, gate_times, gate_opening, jgate);
    
% upstream boundary conditions - compute unit discharge and velocity
    if gate_a > uniform_flow_depth(qus0, Ks(1), S0(1), 'q');
        Dus = uniform_flow_depth(qus0, Ks(1), S0(1), 'q');
        qus = qus0;
    else
        Dus = Cc * gate_a;
        qus = sluice_gate(g, Cc, gate_a, Hgate);        
    end
    
    Uus = qus / Dus;

% downstream boundary conditions - normal depth
    qds = q(end-1);
%     Dds = uniform_flow_depth(qds, Ks(end), S0(1), 'q');
    Dds = 1.1 * critical_flow_depth(g, q(end));
    Uds = qds / Dds;

% PREDICTION STEP

% friction slope
    Se = friction_slope(Ks, U, D);    
    
% prediction terms
    [Dp, Up] = maccormack_prediction(g, dt, dx, S0, n, ...
                                    U, D, U, D, Se, Dus, Uus, Dds, Uds);
    
% ARTIFICIAL VISCOSITY
    [Dp, Up] = maccormack_viscosity(n, Dp, Up);
                                
% CORRECTION STEP

% friction slope
    Se = friction_slope(Ks, Up, Dp);
    
% correction step
	[Dc, Uc] = maccormack_correction(g, dt, dx, S0, n, ...
                                    U, D, Up, Dp, Se, Dus, Uus, Dds, Uds);
    
% SOLUTION STEP
    [D, U] = maccormack_solution(Dp, Up, Dc, Uc);
    
% VARIABLES
    H = z + D;
    q = D .* U;
    Fr = U ./ sqrt(g * D);
    
% update time and iteration counter
    t  = t + dt;
    it = it + 1;
 
% PLOT
    if mod(it, 1000) == 1 || t == TT
        
% print summary
        fprintf('dt = %5.2f s | t = %5.2f h | a = %5.2f m | Du = %5.2f m | q = %5.2f m2/s\n', ...
            dt, t/3600, gate_a, Dus, qus);
    
        subplot(4,1,1)
        plot(x,H);
        hold on
        plot(x,z,'k')
        hold off
%         ylim([0 z(1)+max(gate_opening)])
        title('WSE'); 
        
        subplot(4,1,2)
        plot(x,U);
%         ylim([0 2])
        title('Velocity')
        
        subplot(4,1,3)
        plot(x,q);
%         ylim([0 qus0*1.5])
        title('Discharge')
        
        subplot(4,1,4)
        plot(x,Fr);
        title('Froude number')
        ylim([0 2])
%         pause(0.5);    
        drawnow
    end
    
    if isreal(D) == 0 || isreal(U) == 0
        break
    end
    
% end of the loop over time
end

