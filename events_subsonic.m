function [value,isterminal,direction] = events_subsonic(x,M_sq)
value = (M_sq>=1);  
isterminal = 1;        
direction = 0;  
end
