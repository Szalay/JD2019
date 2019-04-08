classdef EKM
	%EKM Egyker�kmodell
	
	properties (Constant)
		
		m = 1000;			% [kg]
		g = 9.81;			% [m/s^2]
		
		% Ker�k
		R = 0.35;			% [m]
		J = 1.5;			% [kg m^2]
		B = 1;				% [Nm/(rad/s)]
		
		% L�gellen�ll�s
		c_W = 0.3;
		rho = 1.2;			% [kg/m^3]
		A = 2;				% [m^2]
		
		% F�krendszer cm^3 -> kPa
		C_q = 1.4;			% [cm^3/(s*sqrt(kPa))]
		pV0 = 762;			% [kPa]
		pV1 = 2;
		pV2 = 1.41;
		pV3 = 0.097;
		
		V_0 = 0.59;			% [cm^3], haszontalan t�rfogat
		T_H = 10e-3;		% [s], f�kholtid�
		
		A_F = 4e-4;			% [m^2], f�kmunkahenger keresztmetszete
		R_F = 0.2;			% [m], a f�kt�rcsa �tlagos sugara
		mu_F = 1;			% a s�rl�d�si t�nyez� a t�rcsa �s a pofk k�z�tt
		
		p_0 = 20000;		% [kPa]
		
		% Kezdeti �rt�kek
		v_0 = 50/3.6;			% [m/s]
		w_0 = EKM.v_0/EKM.R;	% [rad/s]
	end
	
	properties
		T;
		
		c;
		s_x;
		mu_x;
		
		w;
		v_x;
	end
	
	methods
		
		function this = EKM()
			
			% Szimul�ci�
			[this.T, ~, Y] = sim('ekmsl');
			
			this.c = Y(:, 1);
			this.s_x = Y(:, 2);
			this.mu_x = Y(:, 3);
			
			this.w = Y(:, 4);
			this.v_x = Y(:, 5);
			
		end
		
		function Plot(this)
			figure(489);
			
			% Ker�ksebess�g
			yyaxis left;
			ylim([-5, 40]);
			hold on;
			title('V�szf�kez�s ABS-szel');
			xlabel('Id�, {\itt}, [s]');
			ylabel('Sebess�g �s ker�ksebess�g');
			
			pw = plot(this.T, this.w, 'b-', 'LineWidth', 3);
			pvx = plot(this.T, this.v_x, 'r-', 'LineWidth', 3);
			
			plot(this.T, this.v_x/EKM.R, 'r--', 'LineWidth', 1);
			plot(this.T(this.T >= 0.2), this.v_x(this.T >= 0.2)*(-0.1 + 1)/EKM.R, 'b--');
			
			k = [10, 40];
			plot([0.2, 0.2], k, 'k--');
			text(0.21, mean(k), {'A f�kez�s ', 'kezdete'});
			
			% Cs�sz�s
			yyaxis right;
			ylim([-1, 8]);
			ylabel('Vez�rl�jel, cs�sz�s, s�rl�d�si t�nyez�');
			
			pc = plot(this.T, this.c, 'k-', 'LineWidth', 3);
			psx = plot(this.T, this.s_x, 'y-', 'LineWidth', 3);
			pmu = plot(this.T, this.mu_x, 'm-', 'LineWidth', 3);
			
			legend([pw, pvx, pc, psx, pmu], { ...
				'\omega', 'v_x', 'c', 's_x', '\mu_x' ...
				});
		end
		
	end
	
	methods (Static)
		
		function mu_x = Pacejka(s_x)
			b = 3.76;
			c = 2.7;
			%d = 1;			% mu_max
			d = 0.7;
			e = 1;
			
			% Magic Formula
			mu_x = d * sin( ...
				c*atan( ...
					b * ( (1-e)*s_x + e/b*atan(b*s_x) ) ...
					) ...
				);
		end
		
		function p_2 = F(V)
			if V > EKM.V_0
				p_2 = EKM.pV0 * ( V - EKM.pV1 * (1 - exp( -(V - EKM.pV3)/EKM.pV2 )) );
			else
				% M�g nem �rnek a f�kpof�k a f�kt�rcs�hoz
				p_2 = 0;
			end
		end
		
		function m_f = M_F(w, M, M_F0)
			if w ~= 0
				% A ker�k forog
				m_f = -sign(w) * M_F0;
			else
				% A ker�k nem forog
				if abs(M) > M_F0
					m_f = -sign(M) * M_F0;
				else
					m_f = -sign(M) * M;
				end
			end
		end
		
		function PlotPacejka()
			s_x = -1:0.01:1;
			mu_x = EKM.Pacejka(s_x);
			
			figure(576);
			hold on;
			grid on;
			
			title('A Pacejka-modell szerinti {\it\mu_x}({\its_x}) g�rbe');
			
			xlabel('Cs�sz�s, {\its_x}, [1]');
			ylabel('Hosszir�ny� er��tad�si t�nyez�, {\it\mu_x}, [1]');
			
			plot(s_x, mu_x, 'k-', 'LineWidth', 3);
		end
		
		function ekm = Run()
			
			ekm = EKM();
			ekm.Plot();
			
		end
		
	end
	
end

