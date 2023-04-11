function [X,Y,Z] = fwcylinder(base,tip,diam,N,color,transp)
%[X,Y,Z] = fwcylinder(base,tip,diam,N,color,transp)
%FWCYLINDER plots a cilinder of diameter 'diam' between two points in 3D
% space defined by 3-element vectors 'base' and 'tip'. The smoothnes of
% the surface is defined by the number of points 'N'. The color is
% specified in 3-element vector 'color' and is rendered with transparency
% 'transp'.

% Check inputs
if numel(base) < 3, error('The base of the cylinder must have three values'); end
if numel(tip)  < 3, error('The tip of the cylinder must have three values'); end
if numel(color)< 3, error('Color vector must have three values'); end
if any(color > 1),  error('Color values must be < 1'); end
if any(color < 0),  error('Color values must be > 0'); end
transp= max(0,transp);
transp= min(1,transp);


% Generate base cylinder
[X,Y,Z] = cylinder(diam,N);


% Strech it to length
length= sqrt(abs(sum((tip-base).^2)));
Z= Z*length;


% Rotate to correct angle assuming base is at 0,0,0
tip0= tip - base;
azim= atan2(tip0(2),tip0(1));
elev= atan2(tip0(3),sqrt(tip0(1)^2 + tip0(2)^2));
[X,Y,Z]= fwrotation(X,Y,Z,elev,azim);


% Shift to position
X= X + base(1);
Y= Y + base(2);
Z= Z + base(3);


% Plot the cylinder
surf(X,Y,Z,'EdgeColor','none','FaceColor',color,'FaceAlpha',transp)

end

