function [x,y,z,new_elev,new_azim] = fwrotation(x,y,z,elev,azim)
%[x,y,z,elev,azim] = rotation(x,y,z,elev,azim)
% rotation rotates the position of a point by two angles
% Vertical rotation is relative to the XY plane. Horizontal rotation
% relative to the X axis.


% Rotacion por elevacion
new_elev= atan2(z,x) + (elev-pi/2);
proj_r= sqrt(x.^2 + z.^2);
x= proj_r.* cos(new_elev);
z= proj_r.* sin(new_elev);

% Rotacion por azimut
new_azim= atan2(y,x) + azim;
proj_r= sqrt(x.^2 + y.^2);
x= proj_r .* cos(new_azim);
y= proj_r .* sin(new_azim);

end

   