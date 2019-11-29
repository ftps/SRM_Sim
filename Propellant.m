classdef Propellant<handle
	properties
		cp, ro, M, Tf, r;
		Lg, Dg, Dcore, Seg, b;
	end
	
	methods
		function obj = Propellant(fileName, Lg, Dg, Dcore, Seg, b)
			r = @(p) 0;

			fp = fopen(fileName, 'r');
			
			pp = fscanf(fp, "%f %f %f %f\n", [4, 1]);

			obj.ro = pp(1);
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
			
		end
		
		function [dm, Vc] = burn(obj, pc, Vol, dt)
			dd = 2*obj.r(pc)*dt;
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

			Gvol_new = obj.Seg*obj.Lg*pi*(obj.Dg^2 - obj.Dcore^2)/4;

			if Gvol_new < 0
				Gvol_new = 0;
			end

			dm = obj.ro*(Gvol - Gvol_new);
			if dm < 0
				dm = 0;
			end
			Vc = Vol - Gvol_new;
		end
	end
end

