phi_0 = 0.0;         % Input Polariztion Angle
lambdaRed = 1085e-9;        % WaveLength
lambdaBlue = 729e-9;        % WaveLength
Kb = 1.38e-23;          % Boltzmann constant
a = 250e-9;             % Fiber Radius nanometer
A_in_Red = 1000;
A_in_Blue = 330;
h_bar = 1.054e-34;
c3 = 1.16*h_bar*1e3*(1e6)^3;

field1 = CalculateEvanescentElectricField(lambdaRed,A_in_Red,phi_0,a);
field2 = CalculateEvanescentElectricField(lambdaRed,A_in_Red,phi_0,a);
[field3,r,phi] = CalculateEvanescentElectricField(lambdaBlue,A_in_Blue,phi_0,a);

VanDerWaals = -c3./((r-a).^3);

field = -5*field1-field2+field3;
[x,y] = pol2cart(phi,r);
mesh(x,y,field3-field1-VanDerWaals/50000.0);


