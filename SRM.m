%% Solid Rocket Motor Simulator - by Ftps

% Propellant
prop = "Propellant/sorbitol_fine.br";	% For now, only available propellant
Cc = 0.95;			% Combustion efficiency
n_c = 0.95;         % Real/Ideal density

% Grain Geometry (cylindrical)
Lg = 121e-3;		% Grain length
Dg = 130e-3;			% Grain diamter
Dcore = 60e-3;		% Core diameter
Seg = 3;			% Number of grain segments

% Burn Type
core = true;		% core burning
ends = true;		% ends burning
outer = false;		% outer surface burning

%Nozzle
Dt = 15e-3;			% Nozzle throat diameter
De = 45e-3;			% Nozzle eixt diameter
Cn = 0.5*(1 + cosd(12));	% Nozzle losses (Divergence losses: 0.5*(1+cosd(a)), a = divergence half-angle

% Chamber Geometry
Lc = 400e-3;		% Chamber length
Dc = 130e-3;			% Chamber diameter

% Simulation Time-Step
dt = 1e-4;

% Erosive Burning
K = 0;			% Erosive burning coefficient (K = 0 for no erosive burning, value usualy in 0.05 < K < 0.3)
M_erosive = 0.9;	% Mach number in the chamber at which erosive burning starts happening, usually in 0.8 < M < 1

% Sea level (1) or vacuum (0)
Sea_level = 1;

% Estimated Burn Time
t_est = 4;

%% NO CHANGE FROM HERE

% Simulation Set Up and Run
b = core*(2^0) + ends*(2^1) + outer*(2^2);
m = Motor(prop, Lg, Dg, Dcore, Seg, b, Dt, De, Lc, Dc, Cc, Cn, Sea_level, K, M_erosive, n_c);
m.simulation(dt, t_est);


% Plot Display
figure('Name', 'SRM Simulator', 'NumberTitle', 'off', 'Position', [50, 200 1600, 350])
subplot(1, 3, 1);
plot(m.t, m.Th, 'r');
hold on;
title("Thrust");
xlabel("Time - s");
ylabel("T - N");

subplot(1, 3, 2);
plot(m.t, m.pc, 'b');
title("Chamber Pressure");
xlabel("Time - s");
ylabel("p_c - bar");

subplot(1, 3, 3);
plot(m.t, m.m_dot, 'm');
title("Mass Flow Rate");
xlabel("Time - s");
ylabel("m_d - kg/s");
hold off;

g0 = 9.81;
I = 0;
F_av = 0;
F_max = 0;
[~,l] = size(m.Th);
for i = 1:l
	if m.Th(i) > F_max
		F_max = m.Th(i);
	end
	F_av = F_av + m.Th(i);
	I = I + m.Th(i)*dt;
end
F_av = F_av/l;
Isp = I/(m.p.ro*(Seg*Lg*pi*(Dg^2 - Dcore^2)/4)*g0);

disp(" ");
disp("%%%%% MOTOR DATA %%%%%");
disp(" ");
disp("Total Impulse = " + I + " Ns");
disp("Isp = " + Isp + " s");
disp("Max Thrust = " + F_max + " N");
disp("Average Thrust = " + F_av + " N");
disp("Burn Time = " + m.t_burn + " s");
disp("Thrust Time = " + m.t_t + " s");




