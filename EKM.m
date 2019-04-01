classdef EKM
	%EKM Egykerékmodell
	
	properties (Constant)
		
		m = 1000;			% [kg]
		g = 9.81;			% [m/s^2]
		
		% Kerék
		R = 0.35;			% [m]
		J = 1.5;			% [kg m^2]
		B = 1;				% [Nm/(rad/s)]
		
		% Légellenállás
		c_W = 0.3;
		rho = 1.2;			% [kg/m^3]
		A = 2;				% [m^2]
		
		% Fékrendszer cm^3 -> kPa
		C_q = 1.4;			% [cm^3/(s*sqrt(kPa))]
		pV0 = 762;			% [kPa]
		pV1 = 2;
		pV2 = 1.41;
		pV3 = 0.097;
		
		V_0 = 0.59;			% [cm^3]
		T_H = 5e-3;			% [s]
		
		A_F = 4e-4;			% [m^2], fékmunkahenger keresztmetszete
		R_F = 0.2;			% [m], a féktárcsa átlagos sugara
		mu_F = 1;			% a súrlódási tényezõ a tárcsa és a pofk között
		
		p_0 = 20000;		% [kPa]
		
		% Kezdeti értékek
		v_0 = 20;				% [m/s]
		w_0 = EKM.v_0/EKM.R;	%  [rad/s]
	end
	
	methods (Static)
		
		function mu_x = Pacejka(s_x)
			b = 3.76;
			c = 2.7;
			d = 1;			% mu_max
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
				p_2 = 0;
			end
		end
		
		function m_f = M_F(w, M, M_F0)
			if w ~= 0
				% A kerék forog
				m_f = -sign(w) * M_F0;
			else
				% A kerék nem forog
				if abs(M) > M_F0
					m_f = -sign(M) * M_F0;
				else
					m_f = -sign(M) * M;
				end
			end
		end
		
	end
end

