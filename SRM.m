g0 = 9.81;

% Propellant
prop = "Propellant/sorbitol_fine.br";
Cc = 0.95;

% Grain Geometry
Lg = 100e-3;
Dg = 76e-3;
Dcore = 15e-3;
Seg = 3;

% Burn Type
core = true;	% core burning
ends = true;	% ends burning
outer = false;	% outer surface burning
b = core*(2^0) + ends*(2^1) + outer*(2^2);	% NO CHANGE THIS

%Nozzle
Dt = 15e-3;
De = 40e-3;
Cn = 0.5*(1 + cosd(12));	% Divergence losses

% Chamber Geometry
Lc = 400e-3;
Dc = 76e-3;

% Simulation Time-Step
dt = 1e-3;

% Erosive Burning
K = 0;
M_erosive = 0.8;

%%%%%%% NO CHANGE FROM HERE %%%%%%%

% Simulation Set Up and Run
m = Motor(prop, Lg, Dg, Dcore, Seg, b, Dt, De, Lc, Dc, Cc, Cn, 1, K, M_erosive);
m.simulation(dt);


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




