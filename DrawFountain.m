[x,y,z] = sphere;

figure
surf(x,y,z,'EdgeColor','none','LineStyle','none','FaceLighting','phong')
hold on

for i = 1:30    
    surf(x+5*rand(1)-2.5,y+5*rand(1)-2.5,z+5*rand(1)-2.5,'EdgeColor','none','LineStyle','none','FaceLighting','phong') 
end

for i = 1:100    
    surf(x+2*rand(1)-1,y+2*rand(1)-1,z+i+rand(1)-0.5,'EdgeColor','none','LineStyle','none','FaceLighting','phong') 
end

for i = 1:100
    theta = 2*pi*rand(1);    
    r = 20*rand(1);
    x0 = (2+r)*cos(theta);
    y0 = (2+r)*sin(theta);
    z0 = 0.03*r^2;
    surf(x+x0,y+y0,z+100-z0,'EdgeColor','none','LineStyle','none','FaceLighting','phong') 
end

camlight right
axis equal
colormap([1 0 1])