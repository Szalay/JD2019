classdef KEM < handle
	%KEM Kerékerõ meghatározás
	%   - leírás
	
	properties
		
	end
	
	methods (Static)
		
		function LPA1()
			%LPA1 Lejtõn parkoló autó 1 tengelyen fékezve
			
			syms m g phi l_1 l_2 real;
			
			% Egyenletrendszer, A x = b
			% 2 F_x							     = m g sin(fi)
			%        2 F_z1 + 2 F_z2		     = m g cos(fi)
			%                 2 F_z2 (l_1 + l_2) = m g cos(fi) l_1 + m g sin(fi) l_2
			
			A = [ ...
				2, 0, 0; ...
				0, 2, 2; ...
				0, 0, 2*(l_1 + l_2) ...
				];
			
			b = [ ...
				m*g*sin(phi); ...
				m*g*cos(phi); ...
				m*g*( l_1*sin(phi) + l_2*cos(phi) ) ...
				];
			
			x = A \ b;
			
			x
			
		end
		
	end
end

