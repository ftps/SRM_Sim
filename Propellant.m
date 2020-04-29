classdef Propellant<handle
	properties
		cp, ro, M, Tf, r;
		Lg, Dg, Dc, Seg, Vol, b;
	end

	methods
		function obj = Propellant(fileName, Lg, Dg, Dcore, Seg, b, n_c)
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
			obj.Dc = Dcore;
			obj.Seg = Seg;
			obj.b = b;
            obj.Vol = 0;
		end

		function [mp, Sr] = burn(obj, pc, dt)
            obj.Vol = (pi/4)*(obj.Dg^2 - obj.Dc^2)*obj.Lg;
            Sr = 0;
            if obj.Vol < 0
                mp = 0;
            else
                if bitand(obj.b, 1)
                    Sr = pi*obj.Dc*obj.Lg;
                end
                if bitand(obj.b, 2)
                    Sr = Sr + (pi/2)*(obj.Dg^2 - obj.Dc^2);
                end
                if bitand(obj.b, 4)
                   Sr = Sr + pi*obj.Dg*obj.Lg;
                end
                r = obj.r(pc);
                Sr = obj.Seg*Sr*r;
                mp = obj.ro*Sr;

                if bitand(obj.b, 1)
                    obj.Dc = obj.Dc + 2*r*dt;
                end
                if bitand(obj.b, 2)
                    obj.Lg = obj.Lg - 2*r*dt;
                end
                if bitand(obj.b, 4)
                    obj.Dg = obj.Dg - 2*r*dt;
                end
            end
		end
	end
end
