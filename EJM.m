classdef EJM < handle
	%EJM Egynyomvonal� j�rm�modell
	%
	
	properties
		v = 10;					% [m/s] = 36 km/h
		T_S = 1e-3;				% [s]
		
		T;
		X;						% [beta, dpszi/dt, pszi, x, y]
		Y;						% [a_y]
		
		T_0 = 60;				% [s]
		x_0 = [0; 0; 0; 0; 0];	% [beta; dpszi/dt; pszi; x; y]
		
		% Megjelen�t�s
		Window;
		Menu;
		Verbose = true;
		TrajectoryTitle;
		
		Car;
		VA;						% A j�rm� hossztengelye
		FA;						% Els� tengely
		RA;						% H�ts� tengely
		VV;						% Sebess�gvektor
		AY;						% Oldalir�ny� gyorsul�s
	end
	
	properties (Constant)
		% BMW 3 (F30) adatai r�szben
		
		% Kanyarmerevs�gek a tengelyekre
		c_1 = 3300;			% [N/rad]
		c_2 = 2800;			% [N/rad]
		
		m = 1290;			% [kg]
		J = 1750;			% [kg m^2]
		
		l_1 = 1.2;			% [m]
		l_2 = 1.61;			% [m]
		
		% Tengelyt�v, nyomt�v
		l = EJM.l_1 + EJM.l_2;
		b = 1.543;			% [m]
		
		% Hossz �s a t�ll�g�sok figyelembe v�tele
		L = 4.624;
		L_1 = EJM.l_1 + 0.776;
		L_2 = EJM.l_2 + 1.038;
		
		% Forgat�si m�trix (s�kban)
		R = @(phi) [cos(phi), -sin(phi); sin(phi), cos(phi)];
	end
	
	methods (Static)
		function ejm = Run()
			% Az el�z� futtat�s sor�n l�trehozott v�ltoz� t�rl�se
			evalin('base', 'clear ejm');
			
			% �j p�ld�ny
			ejm = EJM();
			ejm.CreateWindow();
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
			this.Disp('Szimul�ci�');
			
			% Szimul�ci�
			[this.T, this.X] = ode45(@this.Modell, 0:this.T_S:this.T_0, this.x_0);
			
			% Kimenetek
			for i = 1:length(this.T)
				% a_y
				this.Y(i, 1) = ...
					-(EJM.c_1 + EJM.c_2)/(EJM.m * this.v) * this.X(i, 1) + ...
					(EJM.c_2*EJM.l_2 - EJM.c_1*EJM.l_1)/(EJM.m * this.v) * this.X(i, 2) + ...
					EJM.c_1/(EJM.m * this.v) * this.Driver(this.T(i), this.X(i, :)');
			end
		end
		
		function CreateWindow(this)
			if isempty(this.Window) || ~ishandle(this.Window)
				f = 123;
				while ishandle(f)
					f = f + 1;
				end
				this.Window = figure(f);
				this.Window.Name = 'Egynyomvonal� j�rm�modell';
				
				this.Menu = uimenu(this.Window, 'Text', 'Egynyomvonal� j�rm�modell');
				uimenu(this.Menu, 'Text', 'Szimul�ci�', ...
					'MenuSelectedFcn', @this.OnSimulation);
				
				uimenu(this.Menu, 'Text', 'Anim�ci�', ...
					'MenuSelectedFcn', @this.OnAnimation);
				
			end
		end
		
		function OnSimulation(this, ~, ~)
			this.Simulate();
		end
		
		function OnAnimation(this, ~, ~)
			if isempty(this.T)
				this.Simulate();
			end
			
			tr = subplot(1, 2, 1);
			hold on; axis equal;
			this.TrajectoryTitle = title('Halad�s a p�lyag�rb�n');
			
			plot(this.X(:, 4), this.X(:, 5), 'k-', 'LineWidth', 3);
			
			this.Car = patch( ...
				[-EJM.L_2; EJM.L_1; EJM.L_1; -EJM.L_2], ...
				[-EJM.b/2; -EJM.b/2; EJM.b/2; EJM.b/2], ...
				[0, 0, 1] ...
				);
			
			drawnow;
			tr.XLim = [round(tr.XLim(1)/10 - 1)*10, round(tr.XLim(2)/10 + 1)*10];
			
			vm = subplot(1, 2, 2);
			hold on; axis equal;
			title('J�rm�mozg�s');
			vm.XLim = [-3, 3];
			vm.YLim = [-3, 3];
			
			this.VA = plot(0, 0, 'k-', 'LineWidth', 3);
			this.FA = plot(0, 0, 'k-', 'LineWidth', 3);
			this.RA = plot(0, 0, 'k-', 'LineWidth', 3);
			this.VV = plot(0, 0, 'b-', 'LineWidth', 3);
			this.AY = plot(0, 0, 'r-', 'LineWidth', 3);
			this.Vehicle(this.T(1), this.X(1, :), this.Y(1, :));
			
			this.Animate();			
		end
		
		function Animate(this)
			fps = 25;
			
			for i = 1:40:length(this.T)
				if ~ishandle(this.Window)
					break;
				end
				
				this.TrajectoryTitle.String = sprintf( ...
					'Halad�s a p�lyag�rb�n, %3.2f s', this.T(i) ...
					);
				
				xc = [-EJM.L_2; EJM.L_1; EJM.L_1; -EJM.L_2];
				yc = [-EJM.b/2; -EJM.b/2; EJM.b/2; EJM.b/2];
				
				this.Car.Vertices = this.X(i, 4:5) + (EJM.R(this.X(i, 3)) * [xc'; yc'])';
				
				this.Vehicle(this.T(i), this.X(i, :), this.Y(i, :));
				
				drawnow;
				pause(1/fps);
			end
		end
		
		function Vehicle(this, t, x, y)
			
			% Hossztengely
			xy = EJM.R(x(3)) * [-EJM.L_2, EJM.L_1; 0, 0];
			this.VA.XData = xy(1, :);
			this.VA.YData = xy(2, :);
			
			% Els� tengely
			xy = EJM.R(x(3)) * [EJM.l_1, EJM.l_1; EJM.b/2, -EJM.b/2];
			this.FA.XData = xy(1, :);
			this.FA.YData = xy(2, :);
			
			% H�ts� tengely
			xy = EJM.R(x(3)) * [-EJM.l_2, -EJM.l_2; EJM.b/2, -EJM.b/2];
			this.RA.XData = xy(1, :);
			this.RA.YData = xy(2, :);
			
			% Sebess�g vektor, R(b�ta + pszi)
			xy = EJM.R(x(1) + x(3)) * [0, this.v/10; 0, 0];
			this.VV.XData = xy(1, :);
			this.VV.YData = xy(2, :);
			
			% Oldalir�ny� gyorsul�s
			xy = EJM.R(x(3)) * [0, 0; 0, y(1)];
			this.AY.XData = xy(1, :);
			this.AY.YData = xy(2, :);
		end
		
		function Disp(this, str)
			if this.Verbose
				disp(str);
			end
		end
	end
end

