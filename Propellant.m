classdef Propellant<handle
	properties
		cp, ro, M, Tf, r, M_er, K;
		Lg, Dg, Dcore, Seg, b, S;
	end
	
	methods
		function obj = Propellant(fileName, Lg, Dg, Dcore, Seg, b, K, M_er, n_c)
			r = @(p) 0;

			fp = fopen(fileName, 'r');
			
			pp = fscanf(fp, "%f %f %f %f\n", [4, 1]);

			obj.ro = n_c*pp(1);
			obj.Tf = pp(2);
			k = pp(3);
			obj.M = pp(4);
			
			R = 8.3144/obj.M;
			obj.cp = k*R/(k - 1);
			
			pp = fscanf(fp, "%f %f %f\n", [3, Inf]);
			pp = pp';
			[l, ~] = size(pp);

			p_old = 0;

			for i = 1:l
				if i ~= l
					r = @(p) r(p) + (1e-3)*pp(i,2)*((p/10)^pp(i,3))*(p > p_old && p <= pp(i,1));
					p_old = pp(i,1);
				else
					r = @(p) r(p) + (1e-3)*pp(i,2)*((p/10)^pp(i,3))*(p > p_old);
				end
			end
			obj.r = r;
			obj.Lg = Lg;
			obj.Dg = Dg;
			obj.Dcore = Dcore;
			obj.Seg = Seg;
			obj.b = b;
			obj.K = K;
			obj.M_er = M_er;
            obj.S = 0;
            if bitand(obj.b, 1)
				obj.S = pi*obj.Dcore*obj.Lg;
			end
			if bitand(obj.b, 2)
				obj.S = obj.S + (pi/2)*(obj.Dg^2 - obj.Dcore^2);
			end
			if bitand(obj.b, 4)
				obj.S = obj.S + pi*obj.Dg*obj.Lg;
            end
            obj.S = obj.S*obj.Seg;
			
		end
		
		function [dm, Vc] = burn(obj, pc, Vol, dt, At, k, R, Tc)
			
			if obj.K ~= 0
                A_core = pi*(obj.Dcore/2)^2;
                if A_core < At
                    Me = 1;
                else
                    Me = fzero(@(x) (1/x)*sqrt((1+((k-1)/2)*x^2)/(1+((k-1)/2)))^((k+1)/(k-1)) - A_core/At, [1e-4, 1]);
                end

                if Me > obj.M_er
                    dd = 2*obj.r(pc)*(1 + obj.K*(Me - obj.M_er)*sqrt(k*R*Tc/(1+0.5*(k-1))))/1000;
                else
                    dd = 2*obj.r(pc)*dt;
                end
            else
                dd = 2*obj.r(pc)*dt;
            end
			
			Gvol = obj.Seg*obj.Lg*pi*(obj.Dg^2 - obj.Dcore^2)/4;
            
			if bitand(obj.b, 1)
				obj.Dcore = obj.Dcore + dd;
			end
			if bitand(obj.b, 2)
				obj.Lg = obj.Lg - dd;
			end
			if bitand(obj.b, 4)
				obj.Dg = obj.Dg - dd;
            end

            obj.S = 0;
            if bitand(obj.b, 1)
				obj.S = pi*obj.Dcore*obj.Lg;
			end
			if bitand(obj.b, 2)
				obj.S = obj.S + (pi/2)*(obj.Dg^2 - obj.Dcore^2);
			end
			if bitand(obj.b, 4)
				obj.S = obj.S + pi*obj.Dg*obj.Lg;
            end
            obj.S = obj.S*obj.Seg;
            
			Gvol_new = obj.Seg*obj.Lg*pi*(obj.Dg^2 - obj.Dcore^2)/4;

			if Gvol_new < 0
				Gvol_new = 0;
			end

			dm = obj.ro*(Gvol - Gvol_new);
			%dm1 = (obj.ro - pc*(1e5)/(R*Tc))*(Gvol - Gvol_new);
            
            if dm < 0
				dm = 0;
                obj.S = 0;
			end
			Vc = Vol - Gvol_new;
		end
	end
end