function [fitresult, gof] = MeasuredFS()
%CREATEFIT(XM,YM)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : xm
%      Y Output: ym
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 29-Sep-2019 23:01:34


%% Fit: 'untitled fit 1'.
clc
openfig('Group_5.fig')
shading interp
hold on
xm = [0.1648
    0.1558
    0.1513
    0.1441
    0.1298
    0.1280
    0.1217
    0.1136
    0.1073
    0.1010
    0.0948
    0.0858
    0.0822
    0.0777
    0.0723
    0.0669
    0.0588
    0.0561
    0.0516
    0.0480
    0.0436
    0.0391
    0.0364
    0.0301
    0.0274
    0.0193
    0.0175
    0.0130
    0.0076
   -0.0013
   -0.0076
   -0.0148
   -0.0193
   -0.0265
   -0.0346
   -0.0463
   -0.0588
   -0.0678
   -0.0804
   -0.0876
   -0.0948
   -0.1019
   -0.1127
   -0.1190
   -0.1307
   -0.1379
   -0.1432
   -0.1477
   -0.1558
   -0.1612
   -0.1657
   -0.1720
   -0.1828
   -0.1846
   -0.1855
   -0.1855];

ym = [0.1450
    0.1468
    0.1504
    0.1495
    0.1406
    0.1379
    0.1334
    0.1253
    0.1199
    0.1118
    0.1046
    0.0930
    0.0885
    0.0849
    0.0768
    0.0687
    0.0606
    0.0552
    0.0454
    0.0373
    0.0265
    0.0184
    0.0121
    0.0022
   -0.0031
   -0.0103
   -0.0184
   -0.0274
   -0.0382
   -0.0480
   -0.0588
   -0.0714
   -0.0768
   -0.0849
   -0.0948
   -0.1064
   -0.1154
   -0.1199
   -0.1280
   -0.1280
   -0.1298
   -0.1298
   -0.1271
   -0.1271
   -0.1271
   -0.1235
   -0.1226
   -0.1199
   -0.1199
   -0.1172
   -0.1172
   -0.1172
   -0.1199
   -0.1244
   -0.1298
   -0.1343];
syms m_u K x y R d dy                              % SYMBOLIC DECLARATION OF VARIABLES
rho = 2500;                                        % density of particles
g   = 9.81;                                        % gravitational acceleration
R   = 0.238;                                       % Radius of mill
P_crit = 0.6;                                      % 60 % of critical speed
crit_speed = 42.3/(sqrt(2*R));                     % Mill critical speed
speed = P_crit * crit_speed;                       % Mill actual speed in RPM
omega = (speed * 2*pi)/60;                         % omega od drum in rad/sec
eta  = 6e-3;                                       % Granular viscosity
d_p = 3e-3;                                        % Uniform particle diameter
h0 = 1.3*d_p;                                      % Depth of flowing layer
phi = 0.58;                                        % Average solids fraction
theta = deg2rad(40);                               % Angle of repose measured in the anticlockwise direction (Hence 180-40 notation)
p0 = rho*phi*g*h0*cos(theta);   
                                                   % Pressure at depth h0
alpha = 0.4;                                       % fill fraction  
m_u = 0.35;                                        % Mohr columb friction
K = (2*p0^3)/(3*eta * omega * rho^2*g^2);          % Bulk reighnological parameter
%Tangent_vector = linspace(-R,R,100);               % Pre-defining domain for tangent line
robust = 0;                                        % Intersection function usage 
%% =========================================================================================================================
ode = K*((dy)^3 - m_u*dy^2 +dy - m_u) == R^2 - y^2 - x^2;
derivatives = solve(ode,dy,'MaxDegree',3);
root1 = derivatives(1);
%% =========================================================================================================================
%% SOLVING NUMERICALLY FOR THE INFLECTION POINT.....
    eqn1 = K*((-x/y)^3 - m_u*(-x/y)^2 + (-x/y) - m_u) - R^2 + x^2 + y^2;
    d_p = sqrt(x^2 + y^2);
    zeta = acos(d_p/R);
    eqn2 = zeta/pi - (d_p*sin(zeta))/(pi*R) - alpha;
    [X,Y] = vpasolve([eqn1,eqn2],[x y],[0 R;-R R]);
    R1 = matlabFunction(root1,'File','R1_optimized','Optimize',true);                         % optimising the real root
    XX = double(X);
    YY = double(Y);
    XX_check = length(XX);
    if  XX_check == 0
        disp('                                        FAILED TO FIND SOLUTION TO INFLECTION POINT')
    end
    domain = linspace(-R,R,1000);
    domain1 = domain;
    range = sqrt(R^2 - domain.^2);
    range1 = -sqrt(R^2 - domain1.^2);
    [xx,yy] = ode45(R1,[XX,R],YY);
    [xx1,yy1] = ode45(R1,[XX,-R],YY);
    X_intercept_upper = intersection_function(domain,range,xx,yy,robust);
    X_intercept_lower = intersection_function(domain,range1,xx1,yy1,robust);
    [xx,yy] = ode45(R1,[XX,X_intercept_upper(1)],YY);
    [xx1,yy1] = ode45(R1,[XX,X_intercept_lower(1)],YY);
    if length(X_intercept_upper) > 1
        disp('                                                THE DRUM IS CENTRIFUGING')
    end
    plot(xx,yy,'r')
    hold on
    plot(xx1,yy1,'r'); 
    hold on
   
[xData, yData] = prepareCurveData( xm, ym );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.999999731511752;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
 
   
plot( fitresult, xData, yData);
title(['ODE plot with PEPT data study for (\alpha) = ',num2str(alpha)],'FontSize',25); 
legend('off')
clc



