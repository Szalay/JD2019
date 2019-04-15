classdef EJM < handle
	%EJM Egynyomvonal� j�rm�modell
	%
	
	properties
		v = 10;					% [m/s] = 36 km/h
		T_S = 1e-3;				% [s]
		
		T;
		X;						% [beta, dpszi/dt, pszi, x, y]
		Y;
		
		T_0 = 10;				% [s]
		x_0 = [0; 0; 0; 0; 0];	% [beta; dpszi/dt; pszi; x; y]
	end
	
	properties (Constant)
		% Kanyarmerevs�gek a tengelyekre
		c_1 = 3300;			% [N/rad]
		c_2 = 2800;			% [N/rad]
		
		m = 1290;			% [kg]
		J = 1750;			% [kg m^2]
		
		l_1 = 1.2;			% [m]
		l_2 = 1.4;			% [m]
	end
	
	methods (Static)
		function ejm = Run()
			
			ejm = EJM();
		end
	end
	
	methods
		function this = EJM()
			
			
		
		end
		
		function delta = Driver(this, t, x)
			%DRIVER A vezet�
			delta = 5*pi/180;
		end
		
		function dxdt = Modell(this, t, x)
			%MODELL Az j�rm�modell a t�megk�z�ppont mozg�s�val kieg�sz�tve
			% x = [
			%      b�ta;		�sz�si sz�g
			%      dpszi/dt;	legyez�si sz�gsebess�g
			%	   pszi;		legyez�s sz�g
			%      x;			a t�megk�z�ppont helye
			%      y
			%	  ]
			
			% Rendszerm�trix
			A = [ ...
				-(EJM.c_1 + EJM.c_2)/(EJM.m * this.v), ...
				(EJM.c_2*EJM.l_2 - EJM.c_1*EJM.l_1)/(EJM.m * this.v^2) - 1; ...
				(EJM.c_2*EJM.l_2 - EJM.c_1*EJM.l_1)/EJM.J, ...
				-(EJM.c_1*EJM.l_1^2 + EJM.c_2*EJM.l_2^2)/(EJM.J*this.v) ...
				];
			
			% Bemeneti m�trix
			B = [ ...
				EJM.c_1/(EJM.m * this.v); ...
				EJM.c_1*EJM.l_1/EJM.J ...
				];
			
			% Az �llapot deriv�ltja, 5�1
			dxdt = [ ...
				A * x([1, 2], 1); ... % 2�2
				x(2); ...
				this.v * cos(x(3) + x(1)); ...
				this.v * sin(x(3) + x(1)) ...
				] ...
				+ ...
				[B; 0; 0; 0] * this.Driver(t, x);
		end
		
		function Simulate(this)
			% Szimul�ci�
			[this.T, this.X] = ode45(@this.Modell, 0:this.T_S:this.T_0, this.x_0);
		end
		
	end
end

