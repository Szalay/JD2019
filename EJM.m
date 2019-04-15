classdef EJM < handle
	%EJM Egynyomvonalú jármûmodell
	%
	
	properties
		v = 10;					% [m/s] = 36 km/h
		T_S = 1e-3;				% [s]
		
		T;
		X;						% [beta, dpszi/dt, pszi, x, y]
		Y;						% [a_y]
		
		T_0 = 60;				% [s]
		x_0 = [0; 0; 0; 0; 0];	% [beta; dpszi/dt; pszi; x; y]
		
		% Megjelenítés
		Window;
		Menu;
		Verbose = true;
		TrajectoryTitle;
		
		Car;
		VA;						% A jármû hossztengelye
		FA;						% Elsõ tengely
		RA;						% Hátsó tengely
		VV;						% Sebességvektor
		AY;						% Oldalirányú gyorsulás
	end
	
	properties (Constant)
		% BMW 3 (F30) adatai részben
		
		% Kanyarmerevségek a tengelyekre
		c_1 = 3300;			% [N/rad]
		c_2 = 2800;			% [N/rad]
		
		m = 1290;			% [kg]
		J = 1750;			% [kg m^2]
		
		l_1 = 1.2;			% [m]
		l_2 = 1.61;			% [m]
		
		% Tengelytáv, nyomtáv
		l = EJM.l_1 + EJM.l_2;
		b = 1.543;			% [m]
		
		% Hossz és a túllógások figyelembe vétele
		L = 4.624;
		L_1 = EJM.l_1 + 0.776;
		L_2 = EJM.l_2 + 1.038;
		
		% Forgatási mátrix (síkban)
		R = @(phi) [cos(phi), -sin(phi); sin(phi), cos(phi)];
	end
	
	methods (Static)
		function ejm = Run()
			% Az elõzõ futtatás során létrehozott változó törlése
			evalin('base', 'clear ejm');
			
			% Új példány
			ejm = EJM();
			ejm.CreateWindow();
		end
	end
	
	methods
		function this = EJM()
			
		end
		
		function delta = Driver(this, t, x)
			%DRIVER A vezetõ
			delta = 5*pi/180;
		end
		
		function dxdt = Modell(this, t, x)
			%MODELL Az jármûmodell a tömegközéppont mozgásával kiegészítve
			% x = [
			%      béta;		úszási szög
			%      dpszi/dt;	legyezési szögsebesség
			%	   pszi;		legyezés szög
			%      x;			a tömegközéppont helye
			%      y
			%	  ]
			
			% Rendszermátrix
			A = [ ...
				-(EJM.c_1 + EJM.c_2)/(EJM.m * this.v), ...
				(EJM.c_2*EJM.l_2 - EJM.c_1*EJM.l_1)/(EJM.m * this.v^2) - 1; ...
				(EJM.c_2*EJM.l_2 - EJM.c_1*EJM.l_1)/EJM.J, ...
				-(EJM.c_1*EJM.l_1^2 + EJM.c_2*EJM.l_2^2)/(EJM.J*this.v) ...
				];
			
			% Bemeneti mátrix
			B = [ ...
				EJM.c_1/(EJM.m * this.v); ...
				EJM.c_1*EJM.l_1/EJM.J ...
				];
			
			% Az állapot deriváltja, 5×1
			dxdt = [ ...
				A * x([1, 2], 1); ... % 2×2
				x(2); ...
				this.v * cos(x(3) + x(1)); ...
				this.v * sin(x(3) + x(1)) ...
				] ...
				+ ...
				[B; 0; 0; 0] * this.Driver(t, x);
		end
		
		function Simulate(this)
			this.Disp('Szimuláció');
			
			% Szimuláció
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
				this.Window.Name = 'Egynyomvonalú jármûmodell';
				
				this.Menu = uimenu(this.Window, 'Text', 'Egynyomvonalú jármûmodell');
				uimenu(this.Menu, 'Text', 'Szimuláció', ...
					'MenuSelectedFcn', @this.OnSimulation);
				
				uimenu(this.Menu, 'Text', 'Animáció', ...
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
			this.TrajectoryTitle = title('Haladás a pályagörbén');
			
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
			title('Jármûmozgás');
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
					'Haladás a pályagörbén, %3.2f s', this.T(i) ...
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
			
			% Elsõ tengely
			xy = EJM.R(x(3)) * [EJM.l_1, EJM.l_1; EJM.b/2, -EJM.b/2];
			this.FA.XData = xy(1, :);
			this.FA.YData = xy(2, :);
			
			% Hátsó tengely
			xy = EJM.R(x(3)) * [-EJM.l_2, -EJM.l_2; EJM.b/2, -EJM.b/2];
			this.RA.XData = xy(1, :);
			this.RA.YData = xy(2, :);
			
			% Sebesség vektor, R(béta + pszi)
			xy = EJM.R(x(1) + x(3)) * [0, this.v/10; 0, 0];
			this.VV.XData = xy(1, :);
			this.VV.YData = xy(2, :);
			
			% Oldalirányú gyorsulás
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

