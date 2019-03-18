classdef NJM < handle
	%NJM Negyedjármûmodell
	%   
	
	properties (Constant)
		
		g = 9.81;			% [m/s^2]
		
		% Busz
		m_R = 4500;			% [kg]
		c_R = 300000;		% [N/m]
		k_R = 20000;		% [N/(m/s)]
		
		m_0 = 500;
		c_0 = 1600000;		% [N/m]
		k_0 = 150;			% [N/(m/s)]
		
		% Állapot
		% x = [v_R; v_0; z_R, z_0]
		
		% Rendszermátrix
		A = [ ...
			[-NJM.k_R,	NJM.k_R,				-NJM.c_R,	NJM.c_R]/NJM.m_R; ...
			[ NJM.k_R,	-(NJM.k_0 + NJM.k_R),	 NJM.c_R,	-(NJM.c_0 +NJM.c_R)]/NJM.m_0; ...
			1, 0, 0, 0; ...
			0, 1, 0, 0 ...
			];
		
		% Bemenet
		% [F; g; v; z]
		
		% Bemeneti mátrix
		B = [ ...
			 1/NJM.m_R,	-1, 0, 0; ...
			-1/NJM.m_0,	-1, NJM.k_0/NJM.m_0, NJM.c_0/NJM.m_0; ...
			0, 0, 0, 0; ...
			0, 0, 0, 0 ...
			];
		
		% Kimenet
		% y = [a_R; v_R; v_0; z_R, z_0]
		
		% Kimeneti mátrix
		C = [ ...
			[-NJM.k_R,	NJM.k_R,	-NJM.c_R,	NJM.c_R]/NJM.m_R ...
			];
		
		% Elõrecsatolási mátrix
		D = [ ...
			1/NJM.m_R,	-1, 0, 0 ...
			];
		
	end
	
	properties
		u = @(t, x) [0; NJM.g; 0; 0];
		
		% A szimuláció hossza
		T_0 = 1;		% [s]
		
		% A mintavételi idõ
		T_S = 1e-3;		% [s]
		
		% Kezdeti érték
		x_0 = [0; 0; 0; 0];
		
		% A megoldások
		T;
		X;
		Y;
	end
	
	methods
		
		function this = NJM(u, x_0)
			this.u = u;
			this.x_0 = x_0;
			
		end
		
		function Simulate(this)
			
			% Az állapot deriváltja
			dxdt = @(t, x) NJM.A * x + NJM.B * this.u(t, x);
			
			% Megoldás
			[this.T, this.X] = ode45(dxdt, 0:this.T_S:this.T_0, this.x_0);
			
			% Kimenet, y = C x + D u
			for i = 1:length(this.T)
				% y' = x' C' + u' D'
				this.Y(i, :) = this.X(i, :) * NJM.C' + this.u(this.T(i), this.X(i, :)')' * NJM.D';
			end
		end
		
		function Plot(this)
			figure(round(1000*rand(1)));
			
			% Bemenetek
			subplot(2, 1, 1); hold on;
			title('Bemenetek');
			
			U = zeros(length(this.T), 4);
			for i = 1:length(this.T)
				U(i, :) = this.u(this.T(i), this.X(i, :)')';
			end
			
			pf = plot(this.T, U(:, 1), 'LineWidth', 2);
			pg = plot(this.T, U(:, 2), 'LineWidth', 2);
			pv = plot(this.T, U(:, 3), 'LineWidth', 2);
			pz = plot(this.T, U(:, 4), 'LineWidth', 2);
			
			legend([pf, pg, pv, pz], ...
				{'F(t)', 'g', 'v(t)', 'z(t)'} ...
				);
			
			% Állapotok
			subplot(2, 1, 2); hold on;
			title('Állapotok');
			
			pvr = plot(this.T, this.X(:, 1), 'LineWidth', 2);
			pv0 = plot(this.T, this.X(:, 2), 'LineWidth', 2);
			pzr = plot(this.T, this.X(:, 3), 'LineWidth', 2);
			pz0 = plot(this.T, this.X(:, 4), 'LineWidth', 2);
			
			legend([pvr, pv0, pzr, pz0], ...
				{'v_R(t)', 'v_0(t)', 'z_R(t)', 'z_0(t)'} ...
				);
		end
		
	end
	
	methods (Static)
		
		function njm = Run()
			
			% Bemenet
			u = @(t, x) [0; NJM.g; 0; 0];
			
			njm = NJM(u, [0; 0; 0; 0]);
			
			njm.Simulate();
			njm.Plot();
			
		end
		
		function njm = RunWithExcitation()
			
			% Fekvõrendõr
			A = 0.055;	% [m]
			
			w = 2*pi/0.2;
			z = @(t) (t >= 1) * (t <= 1.2) * ( -A*cos(w*t) + A );
			v = @(t) (t >= 1) * (t <= 1.2) * ( A*w*sin(w*t) );
			
			% Bemenet
			u = @(t, x) [0; NJM.g; v(t); z(t)];
			
			njm = NJM(u, NJM.SteadyState([0; NJM.g; 0; 0]));
			njm.T_0 = 5;
			
			njm.Simulate();
			njm.Plot();
			
		end
		
		function [x, y] = SteadyState(u)
			
			% 0 = A x + B u, x = ?
			% A x = -B u
			
			% Állandósult állapot
			x = NJM.A \ (-NJM.B * u);
			
			% Állandósult kimenet
			y = NJM.C * x + NJM.D * u;
			
		end
		
		function Legendre()
			
			L = @(x) 2*x^2 + 29;
			
			p = 1;
			i = 0;
			while p
				N = L(i);
				p = isprime(N);
				
				if p
					disp([num2str(N), ' prím']);
				else
					disp([num2str(N), ' nem prím']);
				end
				
				i = i + 1;
			end
			
		end
		
	end
	
end

