close all
clear all
%% constants

Phi_0 = 0.0;                                 % Input Polariztion Angle
LambdaRed = 1085e-9;                         % Red WaveLength [meter]
LambdaBlue = 729e-9;                         % Blue WaveLength [meter]
ToroidSmallRadius = 750e-9;                  % Fiber Radius [meter]
Kb = 1.38e-23;                               % Joule/Kelvin
A_in_Red = 556e-3;                          % Joule/sec
A_in_Blue = 75e-3;                          % Joule/sec
h_bar = 1.054e-34;                           % meter 
ToroidLargeRadius = 3971.7e-9;               % meter
c3 = 1e6*(1.6e-19*4.9*(1e-10)^3)*(1e9)^3/Kb; % [microKelvin*nm^3]
mode_area = 300e-8/(40*pi);                  % Area occupied by the mode [cm^2] 
ROILength = 2500e-9;                         % length of Region Of Interest away from the fiber
Margin = 15e-9;                              % Too close to the surface the plot explode because of vdw

RadialVector = linspace(ToroidSmallRadius+Margin,ToroidSmallRadius+Margin+ROILength,1000);
Optimization = false;
%% Find optimal values for the intensity of the blue and red fields
if Optimization
    OptimizationN = 501;
    TrapDepth = zeros(OptimizationN,OptimizationN);
    TrapSpatialWidth = zeros(OptimizationN,OptimizationN);
    Power1 = linspace(10,500,OptimizationN) *1e-3 / A_in_Blue;
    Power2 = linspace(100,1500,OptimizationN) *1e-3 / A_in_Red;
    VanDerWaals = -c3./(((RadialVector-ToroidSmallRadius)*1e9).^3);
    AzimuthalVector = 0;
    RedFieldBase = CalculateEvanescentElectricField(LambdaRed,A_in_Red,Phi_0,ToroidSmallRadius,mode_area,RadialVector,0);                       
    [BlueFieldBase,r,phi] = CalculateEvanescentElectricField(LambdaBlue,A_in_Blue,Phi_0,ToroidSmallRadius,mode_area*0.67,RadialVector,0);                    
    for i = 1:OptimizationN         
        for j = 1:OptimizationN            
            RedField = RedFieldBase * Power2(j);
            BlueField = BlueFieldBase * Power1(i);
            TotalField = RedField+BlueField+VanDerWaals;            
            TrapDepth(i,j) = GetTrapDepth(TotalField);
        end
        if mod(i,10)==0
            i
        end
    end
    
    %% plot optimization process results
    Power1 = Power1 * A_in_Blue;
    Power2 = Power2 * A_in_Red;
    imagesc(Power1*1000,Power2*1000, TrapDepth')
    title('Trap Depth(탃) VS Field Intensities');
    xlabel('Blue Field Intensity(mW)');
    ylabel('Red Field Intensity(mW)');
    set(gca,'Ydir','Normal');
    set(gca,'Xdir','Normal');
    set(gca,'FontSize',20)
    colormap jet

    [C,I] = max(TrapDepth(:));
    [blueInd,redInd] = ind2sub(size(TrapDepth),I);    
    
else
    %% calculate for a specific intensity
    AzimuthalVector = linspace(0,2*pi,1000);    
    RedField = CalculateEvanescentElectricField(LambdaRed,A_in_Red,Phi_0,ToroidSmallRadius,mode_area,RadialVector,AzimuthalVector);       
    [BlueField,r,phi] = CalculateEvanescentElectricField(LambdaBlue,A_in_Blue,Phi_0,ToroidSmallRadius,mode_area*0.67,RadialVector,AzimuthalVector);        
    VanDerWaals = -c3./(((r-ToroidSmallRadius)*1e9).^3);
    TotalField = BlueField+RedField+VanDerWaals;
       
    [x,y] = pol2cart(phi,r);
    %% plot fiber cross section
    figure(1)    
    surf(x*1e9,y*1e9,BlueField+RedField,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    view(0,90)
    colormap(jet(1024))
    xlim([-ToroidSmallRadius*2e9 ToroidSmallRadius*2e9])
    ylim([-ToroidSmallRadius*2e9 ToroidSmallRadius*2e9])
    daspect([max(daspect)*[1 1] 1]);
    colorbar();    
    title('Field Cross-Section of a fiber');
    xlabel('X(nm)');
    ylabel('Y(nm)');
    set(gca,'Ydir','Normal');
    set(gca,'FontSize',20)  
    
    %% plot radial dependence at maximum
    figure(2)
    hold on    
    if length(AzimuthalVector) > 1
        TotalFieldSlice = TotalField(:,1);
        plot(r(:,1)',RedField(:,1)','r')
        plot(r(:,1)',BlueField(:,1)','b')
        plot(r(:,1)',VanDerWaals(:,1)','g')
        plot(r(:,1)', TotalField(:,1)','k');
        legend('red','blue','Van Der Waals','Total');
        y1 = ylim;
        ylim([y1(1)/2 max(TotalFieldSlice)*1.1]);        
        title(sprintf('Trap Potential VS Radial Distance. Blue: %d (mW). Red: %d (mW)',A_in_Blue*1000,A_in_Red*1000));
        xlabel('Radial Distance(nm)');
        ylabel('Potentail(탃)');
        set(gca,'Ydir','Normal');
        set(gca,'FontSize',20)  
        axes('position',[.55 .175 .35 .35])
        box on     
        hold on
        title('Trap Inset')
        plot(r(:,1)',RedField(:,1)','r')
        plot(r(:,1)',BlueField(:,1)','b')           
        plot(r(:,1)',VanDerWaals(:,1)','g')
        plot(r(:,1)', TotalField(:,1)','k');
        peak = findpeaks(-TotalFieldSlice(40:end));
        ylim([-peak*1.1 peak/3]);
    else    
        plot(r,RedField,'r')
        plot(r,BlueField,'b')      
        plot(r,VanDerWaals,'g')
        plot(r, TotalField,'k');
        legend('red','blue','Van Der Waals','Total');
        y1 = ylim;
        ylim([y1(1)/2 max(TotalField)*1.1])  
        title(sprintf('Trap Potential VS Radial Distance. Blue: %d (mW). Red: %d (mW)',A_in_Blue*1000,A_in_Red*1000));
        xlabel('Radial Distance(nm)');
        ylabel('Potentail(탃)');
        set(gca,'Ydir','Normal');
        set(gca,'FontSize',20)  
        axes('position',[.55 .175 .35 .35])
        box on     
        hold on
        title('Trap Inset')
        plot(r,RedField,'r')
        plot(r,BlueField,'b')        
        plot(r,VanDerWaals,'g')
        plot(r, TotalField,'k');
        peak = findpeaks(-TotalField(40:end));
        ylim([-peak*1.1 peak/3])  
    end
    
    
    %% plot 3d potential
    figure(3)
    NumOscillationsInToroid = round(2*pi*ToroidLargeRadius/LambdaRed);
    ToroidField = repmat(TotalField(:,1),1,size(phi,2)) .* cos(NumOscillationsInToroid*phi);    
    [x,y] = pol2cart(phi,r+ToroidLargeRadius);
    surf(x*1e9,y*1e9,ToroidField,'EdgeColor','none','LineStyle','none','FaceLighting','phong');    
    colormap(jet(1024))
    title('Total Toroid Potential');
    xlabel('X(nm)');
    ylabel('Y(nm)');  
    zlabel('Potential (탃)');
    set(gca,'FontSize',20)         
    daspect([max(daspect)*[0.01 0.01] 0.1])
    caxis([-peak*1.1 peak/3])
    zlim([-peak*6 peak*7])    
    
    %% plot toroid along with potentail      
    hold on
    ToroidColormap = gray(256);
    ToroidColormap = ToroidColormap(1:200,:);
    ToroidColormap = cat(1,ToroidColormap,flip(ToroidColormap,1));
    [ToroidBodyX,ToroidBodyY,ToroidBodyZ] = Torus(ToroidLargeRadius*1e9,ToroidSmallRadius*1e9,360);
    ToroidBodyZ = ToroidBodyZ/max(max(ToroidBodyZ))*max(max(ToroidBodyX))*ToroidSmallRadius/ToroidLargeRadius*2;
    surf(ToroidBodyX,ToroidBodyY,ToroidBodyZ,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
        
end




