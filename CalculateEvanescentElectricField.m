function [ElectricFieldEnergy,r,phi] = CalculateEvanescentElectricField(lambda,A_in,phi_0,FiberRadius,mode_area,RadialVector,AzimuthalVector)
% Returs the energy of the electric field (µk) and the polar coordinates used to
% plot it correctly.
% The electric field calculations are taken from kimble article - http://iopscience.iop.org/article/10.1088/1367-2630/14/2/023056/pdf
% lambda          - Input field wavelength
% A_in            - Input field intensity
% phi_0           - polarization angle relative to the fiber
% FiberRadius     - Radius of waveguide
% beta            - Mode propagation constant
% mode_area       - Effective area (cm^2) occupied by the mode
% RadialVector    - A vector that holds the points in the radial direction
% AzimuthalVector - A vecotr that holds the angles for the field calculation

n1 = 1.447;                           % Refractive Index Inside Waveguide
n2 = 1.0;                           % Refractive Index Outside Waveguide
k0 = 2*pi/lambda;                   % Wave number
beta = (n1+n2)/2*k0;
h11 = sqrt(k0^2*n1^2-beta^2);       % Characteristic Decay Length Inside Fiber
q11 = sqrt(beta^2-k0^2*n2^2);       % Characteristic Decay Length Outside Fiber
a =FiberRadius;                     % Fiber Radius nanometer (just changing letters for bervity)
Kb = 1.38e-23;                      % Boltzmann constant Joule/Kelvin
Gamma = 5.746e6;                        % [Hz] in the 1*gamma formulations (books, not ofer)
c = 3e8;                            % Speed of light  [m/sec]
finesse = 2.5e4;                    % finnes of the cavity
lambda_res1 = 795e-9;               % rubidium87 resonance D1 [meter]
lambda_res2 = 780e-9;               % rubidium87 resonance D2 [meter]
omega0 = 2*pi*c/lambda_res1;        % Optical transition frequency of D1
omega1 = 2*pi*c/lambda_res2;        % Optical transition frequency of D2
omegaIn = 2*pi*c/lambda;

D_besselk = @(nu,z) 0.5*(besselk(nu-1,z)-besselk(nu+1,z));  % Auxillary function for the derivative of besselk function
D_besselj = @(nu,z) 0.5*(besselj(nu-1,z)-besselj(nu+1,z));  % Auxillary function for the derivative of besselj function

%% Calculate filed intensity I(r)
s11 = (1/(h11*a)^2+1/(q11*a)^2)*(D_besselj(1,h11*a)/(h11*a*besselj(1,h11*a))+...
    D_besselk(1,q11*a)/(q11*a*besselk(1,q11*a)));

A = A_in*beta*besselj(1,h11*a)/(2*q11*besselk(1,q11*a));
B = 1i*A_in*besselj(1,h11*a)/besselk(1,q11*a);

if length(AzimuthalVector) > 1
    [r, phi]  = ndgrid(RadialVector,AzimuthalVector);
else
    phi = AzimuthalVector;
    r = RadialVector;    
end

Ex = A*((1-s11)*besselk(0,q11*r).*cos(phi_0)+(1+s11)*besselk(2,q11*r).*cos(2*phi-phi_0));
Ey = A*((1-s11)*besselk(0,q11*r)*sin(phi_0)+(1+s11)*besselk(2,q11*r).*sin(2*phi-phi_0));
Ez = B*besselk(1,q11*r).*cos(phi-phi_0);

ElectricFieldIntensity = sqrt(abs(Ex).^2+abs(Ey).^2+abs(Ez).^2);

%% Calculate energy U(r)

I = ElectricFieldIntensity/mode_area*finesse;
delta0 = omegaIn - omega0;
delta1 = omegaIn - omega1;
U = (pi*c^2*Gamma/(2*omega0^3))*(2/delta1+1/delta0)*I;
ElectricFieldEnergy = 1e6*U/Kb; % [µK]

end