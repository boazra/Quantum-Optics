function [x,y,z] = Torus(LargeRadius,SmallRadius,n)

[theta, phi]  = ndgrid(linspace(0,2*pi,n),linspace(0,2*pi,n));
x = (LargeRadius + SmallRadius*cos(theta)) .* cos(phi);
y = (LargeRadius + SmallRadius*cos(theta)) .* sin(phi);
z = SmallRadius * sin(theta);

end