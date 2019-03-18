classdef LT < handle
	%LT Lengõ tömeg
	
	% Tulajdonságok
	properties
		
		% Rendszer paraméterei
		J = 1;	% [kg m^2 = Nm / rad/s^2]
		k = 1;	% [Nm / rad/s]
		c = 1;	% [Nm / rad]
		
		% Kezdeti értékek és bemenet
		phi_0 = 0;		% [rad]
		w_0 = 0;		% [rad/s]
		M_K = @(t)0;	% [Nm]
		
		% Eredmények
		T_S = 1e-3;		% [s]
		T_0 = 1;		% [s]
		T = [];			% Idõvektor, oszlop
		X = [];			% Megoldás, [w oszlop, phi oszlop]
		
		% Megjelenítés
		a = 1.5;		% Oldalhossz
		F;				% Az animációs ablak
		L;				% A négyzet
		
		% Forgatási mátrix
		R = @(th)[ ...
			cos(th), -sin(th); ...
			sin(th),  cos(th) ...
			];
		
	end
	
	% Dinamikus tagfüggvények
	methods
		
		% Konstruktor
		function this = LT(st, x_0, varargin)
			% Paraméterek
			this.J = st.J;
			this.k = st.k;
			this.c = st.c;
			
			% Kezdeti értékek
			this.w_0 = x_0(1);
			this.phi_0 = x_0(2);
			
			% Bemenet
			if nargin >= 3
				this.M_K = varargin{1};
			end
		end
		
		function dxdt = Modell(this, t, x)
			% x = [w; phi]
			dxdt = [ ...
				1/this.J * (-this.k*x(1) - this.c*x(2) + this.M_K(t)); ...
				x(1) ...
				];
		end
		
		function Simulation(this)
			[this.T, this.X] = ode45( ...
				@this.Modell, ...
				0:this.T_S:this.T_0, ...
				[this.w_0; this.phi_0] ...
				);
		end
		
		function Plot(this)
			figure(round(1000*rand())); 
			
			subplot(2, 1, 1); hold on;
			title('Szögsebesség');
			plot(this.T, this.X(:, 1), 'b-', 'LineWidth', 3);
			
			subplot(2, 1, 2); hold on;
			title('Elfodulás');
			plot(this.T, this.X(:, 2), 'r-', 'LineWidth', 3);
		end
		
		function Animate(this)
			this.F = figure(round(1000*rand())); 
			hold on;
			t = title('Forgó test, t = 0 s');
			
			axis equal;
			set(gca, 'XLim', this.a*[-2, 2], 'YLim', this.a*[-2, 2]);
			
			% A négyzet egy vonal objektum
			this.L = line('XData', [], 'YData', []);
			
			% A kezdeti helyzetbe forgatás
			this.Rotate(this.phi_0);
			
			% Képkockák elõállítása
			fps = 25;
			
			for i = 1:10:size(this.T, 1)
				% Az idõ kiírása a címbe
				t.String = ['Forgó test, t = ', num2str(this.T(i)) ,' s'];
				
				% Forgatás
				this.Rotate(this.X(i, 2));
				
				% Az ábra újrarajzolásának kényszerítése
				drawnow;
				
				% Várakozás
				pause(1/fps);
			end
			
		end
		
		function Rotate(this, theta)
			x = this.a*[-1, -1, 1, 1, -1]';
			y = this.a*[1, -1, -1, 1, 1]';
			
			v = zeros(length(x), 2);
			for i = 1:length(x)
				v(i, :) = (this.R(theta) * [x(i); y(i)])';
			end
			
			this.L.XData = v(:, 1);
			this.L.YData = v(:, 2);
		end
		
	end
	
	% Statikus tagfüggvények
	methods (Static)
		
		function lt = Run()
			st = struct( ...
				'J', 2, ...
				'k', 1, ...
				'c', 5 ...
				);
			
			phi_0 = 0;
			w_0 = 0;
			x_0 = [w_0; phi_0];
			
			% Bemenet
			M_K = @(t)5*cos(2*pi*0.5*t);
			
			lt = LT(st, x_0, M_K);
			lt.Simulation();
			lt.Plot();
		end
		
		function lt = Animation()
			st = struct( ...
				'J', 2, ...
				'k', 1, ...
				'c', 5 ...
				);
			
			phi_0 = 0;
			w_0 = 0;
			x_0 = [w_0; phi_0];
			
			% Bemenet
			M_K = @(t)25*cos(2*pi*0.5*t);
			
			lt = LT(st, x_0, M_K);
			lt.T_0 = 5;
			lt.Simulation();
			lt.Animate();
		end
		
	end
	
end

