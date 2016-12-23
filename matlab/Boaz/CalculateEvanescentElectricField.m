function [ElectricField,r,phi] = CalculateEvanescentElectricField(lambda,A_in,phi_0,FiberRadius)

n1 = 1.5;            % Refractive Index Inside Waveguide
n2 = 1.0;            % RefractiveIndexOutsideWaveguide
beta = 5.0;            % Mode Propagation Constant
k0 = 2*pi/lambda;    % WaveNumber
h11 = sqrt(k0^2*n1^2-beta^2);        % Characteristic Decay Length Inside Fiber
q11 = sqrt(beta^2-k0^2*n2^2);             % Characteristic Decay Length Outside Fiber
a =FiberRadius;             % Fiber Radius nanometer

D_besselk = @(nu,z) 0.5*(besselk(nu-1,z)-besselk(nu+1,z));
D_besselj = @(nu,z) 0.5*(besselj(nu-1,z)-besselj(nu+1,z));

s11 = (1/(h11*a)^2+1/(q11*a)^2)*(D_besselj(1,h11*a)/(h11*a*besselj(1,h11*a))+...
    D_besselk(1,q11*a)/(q11*a*besselk(1,q11*a)));

A = A_in*beta*besselj(1,h11*a)/(2*q11*besselk(1,q11*a));
B = 1i*A_in*besselj(1,h11*a)/besselk(1,q11*a);

[r, phi]  = ndgrid(linspace(a*1.1,3*a,1000),linspace(0,2*pi,100));

Ex = A*((1-s11)*besselk(0,q11*r).*cos(phi_0)...
    +(1+s11)*besselk(2,q11*r).*cos(2*phi-phi_0));

Ey = A*((1-s11)*besselk(0,q11*r)*sin(phi_0)...
    +(1+s11)*besselk(2,q11*r).*sin(2*phi-phi_0));

Ez = B*besselk(1,q11*r).*cos(phi-phi_0);

ElectricField = abs(Ex.^2+Ey.^2+Ez.^2);



end