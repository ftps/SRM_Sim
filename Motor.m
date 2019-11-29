classdef Motor<handle	
	properties
		amb, p, At, Ae, e, Ceff, Vol;
		pc, Tc, Th, m_dot, m, t;
	end
	
	methods
		function m = Motor(fileName, Lg, Dg, Dcore, Seg, b, Dt, De, Lc, Dc, Cc, Cn, a)
			% ambient conditions
			if a == 1
				% sea level
				m.amb.p = 1;
				m.amb.M = 28.9647e-3;
				k = 1.4;
				R = 8.3144/m.amb.M;
				m.amb.cp = k*R/(k-1);
				m.amb.ro = 1.225;
				m.amb.T = 273.15;
			else
				% vacuum
				m.amb.p = 0;
				m.amb.M = 0;
				m.amb.cp = 0;
				m.amb.ro = 0;
				m.amb.T = 0;
			end
			
			% chemestry and propellant
			m.p = Propellant(fileName, Lg, Dg, Dcore, Seg, b);
			
			% motor geometry
			m.Vol = Lc*pi*(Dc/2)^2;
			m.At = pi*(Dt/2)^2;
			m.Ae = pi*(De/2)^2;
			m.e = m.Ae/m.At;
			m.Ceff = Cc*Cn;
			
			% initial conditions
			m.pc = m.amb.p;
			m.m_dot = 0;
			m.Tc = m.amb.T;
			m.Th = 0;
			m.t = 0;
			Vol_i = Seg*Lg*pi*(Dg^2 - Dcore^2)/4;
			m.m = m.amb.ro*(m.Vol - Vol_i);
		end
		
		function [] = simulation(m, dt)
			phi = 0;
			bar = 1e5;
			cp = m.amb.cp;
			[dm, Vc] = m.p.burn(m.pc(end), m.Vol, dt);
			
			while true
				m.m = [m.m, m.m(end) + dm - m.m_dot(end)*dt];
				phi = (phi*(m.m(end)-dm) + dm)/m.m(end);		% mass percentage of propelant
				M = phi*m.p.M + (1-phi)*m.amb.M;
				R = 8.3144/M;
				cp_old = cp;
				cp = phi*m.p.cp + (1-phi)*m.amb.cp;
				k = cp/(cp-R);
				G = sqrt(k)*(2/(k+1))^((k+1)/(2*(k-1)));
				
				m.Tc = [m.Tc, (dm*m.p.cp*m.p.Tf + (m.m(end-1) - m.m_dot(end)*dt)*cp_old*m.Tc(end))/(m.m(end)*cp)];
				m.pc = [m.pc, m.m(end)*R*m.Tc(end)/(Vc*bar)];
				m.t = [m.t, m.t(end) + dt];
				
				if Vc ~= m.Vol
					[dm, Vc] = m.p.burn(m.pc(end), m.Vol, dt);
				else
					dm = 0;
				end
				
				if m.is_subsonic(k)
					Me = sqrt(2*((m.pc(end)/m.amb.p)^((k-1)/k) - 1)/(k-1));
					m.m_dot = [m.m_dot, m.Ae*Me*m.amb.p*bar*sqrt(k*(1+0.5*(k-1)*Me^2)/(R*m.Tc(end)))];
					m.Th = [m.Th, m.Ceff*m.Ae*m.amb.p*bar*k*Me^2];
				else
					Me = fzero(@(x) (1/x)*sqrt((1+((k-1)/2)*x^2)/(1+((k-1)/2)))^((k+1)/(k-1)) - m.e, [1, 10]);
					tt = 1+0.5*(k-1)*Me^2;
					p_exit = m.pc(end)/(tt^(k/(k-1)));
					m.m_dot = [m.m_dot, m.pc(end)*bar*m.At*G/sqrt(R*m.Tc(end))];
					m.Th = [m.Th, m.Ceff*m.m_dot(end)*Me*sqrt(k*R*m.Tc(end)/tt) + m.Ae*(p_exit - m.amb.p)*bar];
				end
				
				if m.Th(end) < 0
					m.Th(end) = 0;
				end
				
				if abs(m.pc(end) - m.pc(end-1)) < 0.1 && abs(m.Th(end)) < 5
					break;
				end
			end
		end
		
		function [boo] = is_subsonic(m, k)
			Me = fzero(@(x) (1/x)*sqrt((1+((k-1)/2)*x^2)/(1+((k-1)/2)))^((k+1)/(k-1)) - m.e, [1e-4, 1]);
			pp = m.amb.p*((1+0.5*(k-1)*Me^2)^(k/(k-1)));
			if m.pc < pp
				boo = true;
			else
				boo = false;
			end
		end
	end
	
end

