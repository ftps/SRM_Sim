classdef Motor<handle
	properties
		amb, p, At, Ae, e, Ceff, Vol;
		pc, Tc, Th, m_dot, m, t, t_burn, t_size, t_t;
	end

	methods
		function m = Motor(fileName, Lg, Dg, Dcore, Seg, b, Dt, De, Lc, Dc, Cc, Cn, a, n_c)
			% chemestry and propellant
			m.p = Propellant(fileName, Lg, Dg, Dcore, Seg, b, n_c);

			% motor geometry
			m.Vol = Lc*pi*(Dc/2)^2;
			m.At = pi*(Dt/2)^2;
			m.Ae = pi*(De/2)^2;
			m.e = m.Ae/m.At;
			m.Ceff = Cc*Cn;

            % ambient conditions
			if a == 1
				% sea level
				m.amb.p = 1.01325;
				m.amb.M = 28.9647e-3;
				m.amb.k = 1.4;
				m.amb.R = 8.3144/m.amb.M;
				m.amb.cp = m.amb.k*m.amb.R/(m.amb.k-1);
				m.amb.ro = 1.225;
				m.amb.T = 273.15;
			else
				% vacuum
				m.amb.p = 1e-3;
				m.amb.M = 1;
				m.amb.k = 1.4;
				m.amb.R = 8.3144/m.p.M;
				m.amb.cp = m.p.cp;
				m.amb.ro = 0;
				m.amb.T = 0;
			end

			% initial conditions
			m.pc = m.amb.p;
			m.m_dot = 0;
			m.Tc = m.amb.T;
			m.Th = 0;
			m.t = 0;
			Vol_i = Seg*Lg*pi*(Dg^2 - Dcore^2)/4;
			m.m = m.amb.ro*(m.Vol - Vol_i);
			m.t_burn = 0;
            m.t_t = 0;
		end

		function [] = simulation(m, dt, t_est)
			phi = 0;
			bar = 1e5;
            cv = m.amb.cp - m.amb.R;
			%C = 0;
			[mp, Sr] = m.p.burn(m.pc, dt);
			Vc = m.Vol - m.p.Vol;

            extra = floor(t_est/dt);
            max = extra;

            m.pc = [m.pc, zeros(1, max-1)];
            m.Tc = [m.Tc, zeros(1, max-1)];
            m.m_dot = [m.m_dot, zeros(1, max-1)];
            m.Th = [m.Th, zeros(1, max-1)];
            m.t = [m.t, zeros(1, max-1)];
            m.m = [m.m, zeros(1, max-1)];
            i = 1;

			while (m.pc(i) > m.amb.p && m.m(i) > 0) || Sr ~= 0
                i = i+1;
                if i > max

                    max = max + extra;

                    m.pc = [m.pc, zeros(1, extra)];
                    m.Tc = [m.Tc, zeros(1, extra)];
                    m.m_dot = [m.m_dot, zeros(1, extra)];
                    m.Th = [m.Th, zeros(1, extra)];
                    m.t = [m.t, zeros(1, extra)];
                    m.m = [m.m, zeros(1, extra)];
                end

				m.t(i) = m.t(i-1) + dt;
				m.m(i) = m.m(i-1) + (mp - m.m_dot(i-1))*dt;

				if Sr ~= 0
                    %phi = (phi*(m.m(i-1)-dm) + dm)/m.m(i-1);		% mass percentage of propelant
					phi = phi + (1-phi)*mp*dt/m.m(i-1);
                    M = 1/(phi/m.p.M + (1-phi)/m.amb.M);
					R = 8.3144/M;
                    cp = phi*m.p.cp + (1-phi)*m.amb.cp;
                    cv = cp - R;
					k = cp/cv;
                    G = sqrt(k)*(2/(k+1))^((k+1)/(2*(k-1)));
                    dcv = (1-phi)*(mp/m.m(i))*(cp-m.amb.cp + 8.3144*(1/m.p.M - 1/m.amb.M));
                    Vc = Vc + Sr*dt;
					m.Tc(i) = m.Tc(i-1) + ((mp*m.p.cp*m.p.Tf - m.pc(i)*bar*Sr)/(m.m(i)*cv) - m.Tc(i-1)*(k*(m.m_dot(i-1)/m.m(i)) + (mp-m.m_dot(i-1))/m.m(i) + dcv/cv))*dt;
                    %m.pc(i) = m.pc(i-1) + ((R*m.Tc(i-1)/Vc)*(dm1 - m.m_dot(i-1)*dt) + (1/m.Tc(i-1))*(m.Tc(i) - m.Tc(i-1)))/bar;
                    m.pc(i) = m.m(i)*R*m.Tc(i)/(Vc*bar);
                    [mp, Sr] = m.p.burn(m.pc(i), dt);
                else
                    if m.t_burn == 0
                       m.t_burn = m.t(i);
                       %C = m.pc(i-1)/(m.Tc(i-1)^(k/(k-1)));
                    end
					m.Tc(i) = m.Tc(i-1) - (k-1)*m.Tc(i-1)*m.m_dot(i-1)*dt/m.m(i);
                    m.pc(i) = m.m(i)*R*m.Tc(i)/(Vc*bar);
                    %m.pc(i) = C*m.Tc(i)^(k/(k-1));

                end

				Me = fzero(@(x) (1/x)*sqrt((1+((k-1)/2)*x^2)/(1+((k-1)/2)))^((k+1)/(k-1)) - m.e, [1, 10]);
				tt = 1+0.5*(k-1)*Me^2;
				p_exit = m.pc(i)/(tt^(k/(k-1)));
				m.m_dot(i) = bar*m.pc(i)*m.At*G/sqrt(R*m.Tc(i));
                m.Th(i) = m.Ceff*m.m_dot(i)*Me*sqrt(k*R*m.Tc(i)/tt) + m.Ae*(p_exit - m.amb.p)*bar;

				if m.Th(i) < 0
					m.Th(i) = 0;
                    if m.t_t == 0 && Sr == 0
                       m.t_t = m.t(i);
                    end
                end
            end

            m.t = m.t(1:i);
            m.Tc = m.Tc(1:i);
            m.pc = m.pc(1:i);
            m.m = m.m(1:i);
            m.m_dot = m.m_dot(1:i);
            m.Th = m.Th(1:i);
            disp(i);
		end
	end

end
