clearvars
clc
tic

%% INPUTS

num_segments = 200;                  % number of segments
S = 25;                              % Area of wing m^2
AR = 12;                             % aspect ratio
lambda = 0.4;                        % taper ratio
i_w = 2;                             % wind setting angle (deg)
alpha_twist = -1;                    % twist angle
height=50000;                        % altitude of the flight in feet
Velo_start=50;                       % initital velocity of the airfoil in m/s
Velo_end=150;                        % final velocity of the airfoil in m/s
Num_runs=1;                          % how many runs should be completed if 1 will only measure Velo_Start
sweep=0;                             % sweep angle at midchord in degrees
Re_Crit=500000;                      % Critical Reynolds Number
Airfoil = '0012';                    % Type of NACA 4 series airfoil
Flap_L= 0.2;                         % Length of flap in L/c (Breaks down for cambered airfoils >.3 and <.2) We tried
Flap_a= 3;                           % Flap deflection angle must be greater than 0
flapstart = .5;                      % Where on the span do the flaps start? (must be greater than or equal to 0.2) x/b
flapend = .75;                       % Where on the span do the flaps end? (must be less than or equal to 0.8) x/b 
choice = 'n';                        % Use Sadraey flap deflection relationship (Y/N)
LoverD_name = 'L/D of Airfoil';      % Name of L/D plot
LiftDist_Name = 'Lift Distribution'; % Name of lift distribution plot
Planform_name = 'Wing Planform';     % Name of wing planform plot
Table_name = 'Important Values';     % Name of Table

%%
non_flap_1 = flapstart*num_segments; 
non_flap_2 = flapend*num_segments;
a_2d = 2*pi;                        % lift curve slope (1/rad)
font = 13; %Font size for all plots                                            
N=num_segments-1; %number of segments -1
b = sqrt(AR*S); %wing span (m)
MAC = S/b; %mean aerodynamic chord(m)
Croot = (1.5*(1+lambda)*MAC)/(1+lambda+lambda^2); %root chord (m)
theta = pi/(2*N):pi/(2*N):pi/2;
alpha = i_w+alpha_twist:-alpha_twist/(N-1):i_w;
    %segment's angle of attack
step=.000025;   
%%
thickness= Airfoil;
thickness(1:2)=[];
Tc= str2num(thickness)*.01;
camber=str2num(Airfoil(1))*.01;
camberLoc=str2num(Airfoil(2))*.1;
if Airfoil(1) == '0'
    camberLoc=.3;
end
if Airfoil(1) == '0'
alpha_0 = 0;
Wing_Length=cosd(asind((sind(Flap_a)*Flap_L)/(1-Flap_L)))*(1-Flap_L);
Height=sind(Flap_a)*Flap_L;
FlapLength=cosd(Flap_a)*Flap_L;
Slope_1=Height/Wing_Length ;
Slope_2=-Height/FlapLength;
FlapStart_c=acos(1-2*(Wing_Length/(Wing_Length+FlapLength)));
alpha_flap_0=(-1/pi)*((Slope_1*(sin(FlapStart_c)-sin(0)-FlapStart_c+0))+(Slope_2*(sin(pi)-sin(FlapStart_c)-pi+FlapStart_c)))*(180/pi);
else
start_angle=atand(camber/(1-camberLoc));
Flap_true_L=Flap_L/(cosd(start_angle));
height_change_flap=Flap_true_L*(sind(Flap_a));
second_half_L=(1-camberLoc-Flap_L)/(cosd(start_angle));
angle_change=atand(height_change_flap/second_half_L);
new_length_pt1=camberLoc*cosd(angle_change);
new_length_pt2=(1-camberLoc-Flap_L)*cosd(angle_change);
new_length_flap=Flap_L*cosd(Flap_a);
new_length_total=new_length_pt1+new_length_pt2+new_length_flap;
counting=1;
for i=0:step:camberLoc
    y_camber_1(counting)=(camber/(camberLoc^2))*(2*camberLoc*i-i^2);
    x_camber_1(counting)=acos(1-(2*i));
    counting=counting+1;
end
counting=1;
for i=camberLoc:step:(1-Flap_L)
    y_camber_2(counting)=(camber/((1-camberLoc)^2))*((1-2*camberLoc)+2*camberLoc*i-i^2);
    x_camber_2(counting)=acos(1-(2*i));
    counting=counting+1;
end
counting=1;
for i=(1-Flap_L):step:1
    y_camber_3(counting)=(camber/((1-camberLoc)^2))*((1-2*camberLoc)+2*camberLoc*i-i^2);
    x_camber_3(counting)=acos(1-(2*i));
    counting=counting+1;
end
x_camber_1_new=((x_camber_1*cosd(angle_change))-(y_camber_1*sind(angle_change)))/new_length_total;
y_camber_1_new=((x_camber_1*sind(angle_change))+(y_camber_1*cosd(angle_change)))/new_length_total;
x_camber_2_new=((x_camber_2*cosd(angle_change))-(y_camber_2*sind(angle_change)))/new_length_total;
y_camber_2_new=((x_camber_2*sind(angle_change))+(y_camber_2*cosd(angle_change)))/new_length_total;
x_camber_3_new=((x_camber_3*cosd(-1*Flap_a))-(y_camber_3*sind(-1*Flap_a)))/new_length_total;
y_camber_3_n=((((x_camber_3*sind(-1*Flap_a*((camberLoc*10)^.5)))+(y_camber_3*cosd(-1*Flap_a*((camberLoc*10)^.5))))/new_length_total)); %best attempt to get near it
y_camber_3_new=y_camber_3_n-y_camber_3_n(length(y_camber_3_n));
slope1=y_camber_1_new(length(y_camber_1_new))/x_camber_1_new(length(x_camber_1_new));
slope2=(y_camber_2_new(length(y_camber_2_new))-y_camber_1_new(length(y_camber_1_new)))/(x_camber_2_new(length(x_camber_2_new))-x_camber_1_new(length(x_camber_1_new)));
slope3=y_camber_2_new(length(y_camber_2_new))/(1-x_camber_2_new(length(x_camber_2_new)));
alpha_flap_0=(-1/pi)*(slope1*(sin(y_camber_1_new(length(y_camber_1_new)))-sin(0)-(y_camber_1_new(length(y_camber_1_new)))+0))+(slope2*(sin(y_camber_2_new(length(y_camber_2_new)))-sin((y_camber_1_new(length(y_camber_1_new))))-(y_camber_2_new(length(y_camber_2_new)))+(y_camber_1_new(length(y_camber_1_new)))))+(slope3*(sin(pi)-sin((y_camber_2_new(length(y_camber_2_new))))-pi+(y_camber_2_new(length(y_camber_2_new)))))*(180/pi)*-1;
syms C; %uses gradient of camber + thin airfoil theroy to determine unflapped alfa 0 lift
M=camber;
p=camberLoc;
xTRANS=(1-cos(C))/2;
dydx1=((2*M)/(p^2))*(p-xTRANS);
dydx2=((2*M)/((1-p)^2))*(p-xTRANS);
important_location=acos(1-(2*p));
Upper = int(dydx1*(cos(C)-1)*(-1/pi),[0 important_location]);
Downer = int(dydx2*(cos(C)-1)*(-1/pi),[important_location pi]);
symbolicGARBAGE=(Upper+Downer)*(180/pi);
alpha_0 = double(symbolicGARBAGE);
end
if Flap_a == 0
    alpha_flap_0=alpha_0
end
% alpha_0       =   0 lift angle of attack for locations without flaps
% alpha_flap_0  =   0 lift angle of attack for locations with flaps

if choice == 'y' || choice == 'Y'
    alpha_flap_0=(-1.15*(Flap_L)*Flap_a)+alpha_0;
end

%%

for q=0:1:Num_runs-1

Velocity=Velo_start+(q*((Velo_end-Velo_start)/(Num_runs-1)));

if Num_runs == 1
    Velocity=Velo_start;
end

z = (b/2)*cos(theta);
c = Croot * (1-(1-lambda)*cos(theta)); %mean aerodynamics chord at each segment (m)
mu = c * a_2d / (4*b);
%%

for i = 1:non_flap_1
LHS(i) = mu(i) .* (alpha(i)-alpha_0)/57.3;
end
for i = (non_flap_1 + 1):non_flap_2
LHS(i) = mu(i) .* (alpha(i)-alpha_flap_0)/57.3;
end
for i = (non_flap_2 + 1):num_segments-1
LHS(i) = mu(i)  .* (alpha(i)-alpha_0)/57.3;
end
%%
%solving N equations to find coefficients A(i);
for i = 1:N
    for j = 1:N
        B(i,j) = sin((2*j-1)*theta(i))*(1+(mu(i)*(2*j-1))/sin(theta(i)));
    end
end
A=B\transpose(LHS);
for i = 1:N
    sum1(i) = 0;
    sum2(i) = 0;
    for j = 1:N
        sum1(i) = sum1(i) + (2*j-1) * A(j)*sin((2*j-1)*theta(i));
        sum2(i) = sum2(i) + A(j)*sin((2*j-1)*theta(i));
    end
end
CL = 4*b*sum2 ./ c;
CL1(1)=0;
y_s(1)=b/2;
for i = 1:1:N
    CL1(i+1)=CL(i);
    y_s(i+1)=z(i);
end
CL_wing = pi *AR * A(1);
%%
if height < 36152 %all from https://www.grc.nasa.gov/WWW/K-12/rocket/atmos.html
    temp=(59-(.00356*height));
    pressure = 2116*((temp+459.7)/518.6)^5.256;
    density = (pressure/(1718*(temp+459.7)))*515.379;
else
    temp=-70;
    pressure=473.1*2.718281828459045^(1.73-(.000048*height));
    density = (pressure/(1718*(temp+459.7)))*515.379;
end
p1 =  -6.947e-42;  %Curve fit ran with dynamic viscocity to atmospheric height
       p2 =    1.76e-36;  
       p3 =  -1.527e-31 ; 
       p4 =   3.413e-27  ;
       p5 =    2.57e-22 ;
       p6 =  -1.874e-17  ;
       p7 =   4.628e-13  ;
       p8 =  -3.637e-09  ;
       p9 =  -2.559e-05  ;
       p10 =       1.793  ;
visc_d = (p1*height^9 + p2*height^8 + p3*height^7 + p4*height^6 + p5*height^5 + p6*height^4 + p7*height^3 + p8*height^2 + p9*height + p10)*10^-5;
thickness= Airfoil;
thickness(1:2)=[];
Tc= str2num(thickness)*.01;
camber=str2num(Airfoil(1))*.01;
camberLoc=str2num(Airfoil(2))*.1;
if Airfoil(2) == '0'
    camberLoc=.3;
end
if height < 36000 %determines speed of sound based on altitude
    SOS=0.514444*(29.06*((518.7-(3.57*(height/1000)))^.5));
else
    SOS=0.514444*(29.06*((518.7-(3.57*(36000/1000)))^.5));
end
    
id_coff=0; %Finds induced drag coeffient
for i=3:2:N
    id_coff=id_coff+(A(i)/A(1))^2;
end
F=(1+((.6/camberLoc)*(Tc))+(100*((Tc)^4)))*(1.34*((Velocity/SOS)^.18)*((cosd(sweep))^.28)); %Wing form factor
Re=(Velocity*MAC*density)/visc_d; %wing reynolds number
Cf_Turb=0.074/((Re)^.2); %turbulent skin friction
Cf_Lam=1.328/((Re)^.5); %laminar skin friction
xcr=((visc_d*Re_Crit)/(density*Velocity))/MAC; %critical reynolds number transition location
skin_Drag_Coeff=(xcr*Cf_Lam)+Cf_Turb-(xcr*Cf_Turb); %skin drag per example 4.10 in anderson text
pressure_drag_coeff=skin_Drag_Coeff*(F-1); %pressure drag from form factor
Drag_Induced_coeff=((CL_wing^2)/(pi*AR))*(1+id_coff); %wing induced drag
drag_coeff=skin_Drag_Coeff+pressure_drag_coeff+Drag_Induced_coeff; %total drag coeff

if Flap_a ~= 0

k1 = 0.0007142857 + 4.160714*Flap_L + 1.785714*(Flap_L^2) + 12.5*(Flap_L^3); % one single outcome generalized for plain flaps and all t/c
k2 = -0.001757651 + 0.001567825*Flap_a + 0.00001560353*(Flap_a^2); % one single outcome generalized for plain flaps and all t/c

drag_coeff=drag_coeff+(k1*k2*((MAC*((flapend-flapstart)*(b/2)))/S)); %additional flap drag
% Need the length over the wing that the flaps occupy

end

L_D = CL_wing/drag_coeff; % Wing lift/drag value
LD(q+1)=L_D;
LDx(q+1)=Velocity;

if Num_runs > 1

fprintf('\nRun %.0f\n',q+1)
fprintf('Velocity    = %.1f m/s \n',Velocity)
fprintf('Drag Coef   = %f \n',drag_coeff)
fprintf('Lift Coef   = %f \n',CL_wing)
fprintf('Lift/Drag   = %f \n',L_D)
fprintf('\n')
end

end
%%

y_Planform=[0,c(length(c)),(c(length(c))*.5)+(.5*c(1))-(sind(sweep)*b/2),(c(length(c))*.5)-(.5*c(1))-(sind(sweep)*b/2),0];
x_Planform=[0,0,b/2,b/2,0];
figure(1)
tiledlayout(2,1)
nexttile %plot of lift curve
plot(y_s/b,CL1,'LineWidth',2.25,'Color',[0.4940 0.1840 0.5560])
hold on
plot(-y_s/b,CL1,'LineWidth',2.25,'Color',[0.4940 0.1840 0.5560])
xlim([-b/2*1.05/b b/2*1.05/b])
hold off
 grid
title(LiftDist_Name,'FontSize',font)
xlabel('Wingspan location','FontSize',font)
ylabel('Lift coefficient','FontSize',font)
nexttile %plot of wing planform
plot(x_Planform/b,y_Planform/b,'LineWidth',2.25,'Color',[0 0.4470 0.7410])
hold on
plot(-1*x_Planform/b,y_Planform/b,'LineWidth',2.25,'Color',[0 0.4470 0.7410])
hold off
daspect([1 1 1])
xlim([-b/2*1.05/b b/2*1.05/b])
if (c(length(c))*.5)-(.5*c(1))-(sind(sweep)*b/2)<0
    ylim([((c(length(c))*.5)-(.5*c(1))-(sind(sweep)*b/2)-(c(length(c)))*.05)/b (c(length(c))*1.05)/b])
else
    ylim([0-((c(length(c)))*.05)/b (c(length(c))*1.05)/b])
end
title(Planform_name,'FontSize',font)
ylabel('Chord Length','FontSize',font)
xlabel('Wing Span','FontSize',font)
set(gca,'FontSize',font)

if Num_runs > 1
    figure(2)
plot(LDx,LD,'LineWidth',2.25,'Color',[0.9290 0.6940 0.1250])
xlabel('Velocity in m/s','FontSize',font)
ylabel('Lift/Drag','FontSize',font)
title(LoverD_name,'FontSize',font)
end

%%
fprintf(Table_name)
fprintf('\n')
fprintf('Airfoil is NACA %s \n',Airfoil)
fprintf('Wingspan    = %f m \n',b)
fprintf('Mean Chord  = %f m \n',MAC)
fprintf('Taper Ratio = %f \n',lambda)
if (sweep~= 0)
fprintf('Mid Sweep   = %f degrees \n',sweep)
end
fprintf('Altitude    = %f m \n',height*.3048)
fprintf('Twist Angle = %f degrees \n',alpha_twist)
fprintf('Wind Angle  = %f degrees \n',i_w)
if (Flap_a ~= 0)
fprintf('Flap Length = %f m \n',MAC*Flap_L)
fprintf('Flap Angle  = %f degrees \n',Flap_a)
end
fprintf('α L=0 Wing  = %f degrees \n',alpha_0)
if (Flap_a ~= 0)
fprintf('α L=0 Flap  = %f degrees \n',alpha_flap_0)
fprintf('Flap Start  = %f m\n',flapstart*(b/2))
fprintf('Flap End    = %f m\n',flapend*(b/2))
end

if Num_runs == 1
fprintf('Velocity    = %.1f m/s \n',Velocity)
fprintf('Drag Coef   = %f \n',drag_coeff)
fprintf('Lift Coef   = %f \n',CL_wing)
fprintf('Lift/Drag   = %f \n',L_D)
end
fprintf('\n')
toc
