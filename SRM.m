%% Solid Rocket Motor Simulator - by Ftps

clear

% Propellant
prop = "Propellant/sorbitol_coarse.br";	% For now, only available propellant
Cc = 0.95;			% Combustion efficiency
n_c = 0.95;         % Real/Ideal density

% Grain Geometry (cylindrical)
Lg = 107e-3;		% Grain length
Dg = 76.6e-3;			% Grain diamter
Dcore = 25e-3;		% Core diameter
Seg = 1;			% Number of grain segments

% Burn Type
core = true;		% core burning
ends = false;		% ends burning
outer = false;		% outer surface burning

%Nozzle
Dt = 7e-3;			% Nozzle throat diameter
De = 14e-3;			% Nozzle eixt diameter
Cn = 0.5*(1 + cosd(12));	% Nozzle losses (Divergence losses: 0.5*(1+cosd(a)), a = divergence half-angle

% Chamber Geometry
Lc = 118.9e-3;		% Chamber length
Dc = 80e-3;			% Chamber diameter

% Simulation Time-Step
dt = 1e-3;

% Sea level (1) or vacuum (0)
Sea_level = 1;      % DON'T CHANGE, NO BURN MODEL FOR VACCUM PRESSURES AND NULL-MASS ERRORS

% Estimated Burn Time
t_est = 4;

%% NO CHANGE FROM HERE

% Simulation Set Up and Run
b = core*(2^0) + ends*(2^1) + outer*(2^2);
m = Motor(prop, Lg, Dg, Dcore, Seg, b, Dt, De, Lc, Dc, Cc, Cn, Sea_level, n_c);
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
