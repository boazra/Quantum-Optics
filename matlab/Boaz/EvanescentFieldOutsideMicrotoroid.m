close all
clear all

Phi_0 = 0.0;                                 % Input Polariztion Angle
LambdaRed = 1085e-9;                         % Red WaveLength [meter]
LambdaBlue = 729e-9;                         % Blue WaveLength [meter]
ToroidSmallRadius = 250e-9;                  % Fiber Radius [meter]
Kb = 1.38e-23;                               % Joule/Kelvin
A_in_Red = 1000e-3;                          % Joule/sec
A_in_Blue = 330e-3;                         % Joule/sec
h_bar = 1.054e-34;                           % meter 
ToroidLargeRadius = 3971.7e-9;               % meter
c3 = 1e6*(1.6e-19*4.9*(1e-10)^3)*(1e9)^3/Kb; % [microKelvin*nm^3]
lambda_res1 = 780e-9;                        % Atom resonance 1 [meter]
lambda_res2 = 795e-9;                        % Atom resonance 2 [meter]
mode_area = 300e-8/(40*pi);                  % Area occupied by the mode [cm^2] 
ROILength = 2000e-9;                         % length of Region Of Interest away from the fiber
Margin = 33e-9;                              % Too close to the surface the plot explode because of vdw
RadialVector = linspace(ToroidSmallRadius+Margin,ToroidSmallRadius+Margin+ROILength,1000);
NumOscillationsInToroid = round(2*pi*ToroidLargeRadius/LambdaRed); 

RedField1 = CalculateEvanescentElectricField(LambdaRed,lambda_res1,A_in_Red,Phi_0,ToroidSmallRadius,1.0,mode_area,RadialVector);
RedField2 = CalculateEvanescentElectricField(LambdaRed,lambda_res2,A_in_Red,Phi_0,ToroidSmallRadius,1.0,mode_area,RadialVector);
BlueField1 = CalculateEvanescentElectricField(LambdaBlue,lambda_res1,A_in_Blue,Phi_0,ToroidSmallRadius,1.0,mode_area*0.67,RadialVector);
[BlueField2,r,phi] = CalculateEvanescentElectricField(LambdaBlue,lambda_res2,A_in_Blue,Phi_0,ToroidSmallRadius,1.0,mode_area*0.67,RadialVector);
VanDerWaals = -c3./(((r-ToroidSmallRadius)*1e9).^3);
TotalField = BlueField1+RedField1+BlueField2+RedField2+VanDerWaals;

ToroidField = repmat(TotalField(:,1),1,size(phi,2)) .* cos(NumOscillationsInToroid*phi);
[x,y] = pol2cart(phi,r+ToroidLargeRadius);

figure(1)
surf(x,y,ToroidField,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
colormap jet

figure(2)
hold on
plot(r(:,1),RedField1(:,1),'r')
plot(r(:,1),BlueField1(:,1),'b')
plot(r(:,1),RedField2(:,1),'r--')
plot(r(:,1),BlueField2(:,1),'b--')
plot(r(:,1),VanDerWaals(:,1),'g')
plot(r(:,1), ToroidField(:,1),'k');

