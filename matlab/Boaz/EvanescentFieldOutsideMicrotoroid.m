phi_0 = 0.0;         % Input Polariztion Angle
lambdaRed = 1085e-9;        % WaveLength
lambdaBlue = 729e-9;        % WaveLength
a = 250e-9;             % Fiber Radius nanometer
A_in_Red = 1000;
A_in_Blue = 330;
h_bar = 1.054e-34;
c3 = 1.16*h_bar*1e3*(1e6)^3/50000.0; % TODO: find c3
ToroidRadius = 3971.7e-9;
k = round(2*pi*ToroidRadius/lambdaRed);

field1 = CalculateEvanescentElectricField(lambdaRed,A_in_Red,phi_0,a,1.0);
[field3,r,phi] = CalculateEvanescentElectricField(lambdaBlue,A_in_Blue,phi_0,a,1.0);

VanDerWaals = -c3./((r-a).^3);
field = field3-field1-VanDerWaals;
ToroidField = repmat(field(:,1),1,size(phi,2)) .* cos(k*phi);
[x,y] = pol2cart(phi,r);
mesh(x,y,ToroidField);
colormap jet
shading interp;