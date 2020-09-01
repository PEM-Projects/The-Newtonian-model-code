%% THIS THE MAIN SCRIPT IN WHICH THE ENTIRE CODE IS EMBBEDED. ALL SUBSEQUENT FUNCTIONS ARE INVOKED IN THIS SCRIPT; KIND REGARDS, GROUP 5
%% =======================================================================================================================================
clc;
clear figure()
clear;
close all;
MeasuredFS();
figure();
%%                                                  TO GOD BE THE GLORY
%% GROUP MEMBERS
%% 1.  Phelelani Eshmael Mamba (217058111)
%% 2.  Tashmira Ramjan         (217003070)
%% 3.  Taryn Gore              (217003065)
%% 4.  Prishalan Moodley       (217002787)
%% ========================================================================================================================================
%%                                                 THE DECLARATION OF VARIABLES
%% VARIABLE                                        VARIABLE DESCRIPTION
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
theta = deg2rad(40);                               % Angle of repose measured in the anticlockwise direction 
p0 = rho*phi*g*h0*cos(theta);                      % Pressure at depth h0
alpha = 0.25;                                       % fill fraction  
m_u = 0.35;                                        % Mohr columb friction
K = (2*p0^3)/(3*eta * omega * rho^2*g^2);          % Bulk reighnological parameter
Tangent_vector = linspace(-R,R,100);               % Pre-defining domain for tangent line
org_x = 0;                                         % X point of the origin
org_y = 0;                                         % Y point of the origin
robust = 0;                                        % Intersection function usage 
%% =========================================================================================================================
ode = K*((dy)^3 - m_u*dy^2 +dy - m_u) == R^2 - y^2 - x^2;
derivatives = solve(ode,dy,'MaxDegree',3);
root1 = derivatives(1);
root2 = derivatives(2);
root3 = derivatives(3);

%% =========================================================================================================================
%% SOLVING NUMERICALLY FOR THE INFLECTION POINT.....
if alpha < 0.5 && alpha >=0
    disp('                                     SOLVING INFLECTION POINT FOR ALPHA LESS THAN 0.5  ')
    eqn1 = K*((-x/y)^3 - m_u*(-x/y)^2 + (-x/y) - m_u) - R^2 + x^2 + y^2;
    d = sqrt(x^2 + y^2);
    zeta = acos(d/R);
    eqn2 = zeta/pi - (d*sin(zeta))/(pi*R) - alpha;
    [X,Y] = vpasolve([eqn1,eqn2],[x y],[0 R;-R R]);
    R1 = matlabFunction(root1,'File','R1_optimized','Optimize',true);                         % optimising the real root
    XX = double(X);
    YY = double(Y);
    XX_check = length(XX);
    if  XX_check == 0
        disp('                                        FAILED TO FIND SOLUTION TO INFLECTION POINT')
    end
    R_vector = linspace(0,XX,50);
    grad_diam = ((0-YY)/(0-XX));
    diam_line = grad_diam.*R_vector;
    Tangent_grad = -1/grad_diam;
    Y_int = YY - Tangent_grad*XX;
    Tangent_line = Tangent_grad.*Tangent_vector + Y_int;
    angle_of_repose = abs(rad2deg(atan(-1/grad_diam)))
    statictics = odeset('Stats','on');
    [xx,yy] = ode45(R1,[XX,R],YY);
    [xx1,yy1] = ode45(R1,[XX,-R],YY);
    domain = linspace(-R,R,1000);
    domain1 = domain;
    range = sqrt(R^2 - domain.^2);
    range1 = -sqrt(R^2 - domain1.^2);
    Tangent_line_intercept_upper = intersection_function(domain,range,Tangent_vector,Tangent_line,robust);
    Tangent_line_intercept_lower = intersection_function(domain,range1,Tangent_vector,Tangent_line,robust);
   [X_intercept_upper,yu] = intersection_function(domain,range,xx,yy,robust)
    if length(X_intercept_upper) > 1
        disp('                                                THE DRUM IS CENTRIFUGING')
    end
    [X_intercept_lower,ydn] = intersection_function(domain,range1,xx1,yy1,robust)
    
    Tangent_vector = linspace(Tangent_line_intercept_lower,Tangent_line_intercept_upper,200);
    Tangent_line = Tangent_grad.*Tangent_vector + Y_int;
    line1_domain = linspace(0,R,100);
    line2_domain = linspace(0,Tangent_line_intercept_upper,100);
        i = 1;
    while xx(i) < X_intercept_upper 
        XX_vector(i) = xx(i);
        i = i +1;
    end
    YY_vector = yy(1:length(XX_vector));
    plot(XX_vector,YY_vector,'color','g','LineWidth',2);
    hold on
    j = 1;
    while xx1(j) > X_intercept_lower
        XX1_vector(j) = xx1(j);
        j = j +1;
    end
    YY1_vector = yy1(1:length(XX1_vector));
    plot(XX1_vector,YY1_vector,'color','r','LineWidth',2);
    len = sqrt(XX^2 + YY^2);
    R = double(R);
    zetDa =  acos(len/R);
    zetha = deg2rad(angle_of_repose) + pi/2 - zetDa;
    zetha = double(zetha);
    x2=(-R*cos(zetha));
    y2=(-R*sin(zetha));
    Y_tang_line = Tangent_grad*Tangent_line_intercept_upper + Y_int;
    R_vector1 = [0 Tangent_line_intercept_upper];
    Tang_line = ((0 - Y_tang_line)/(0 - Tangent_line_intercept_upper)).*R_vector1; 
    plot(domain,range,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
    plot(domain,range1,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    plot([0 x2],[0 y2],'color','k','HandleVisibility','off')
    hold on
    plot(R_vector1,Tang_line,'color','k','HandleVisibility','off')
    hold on
    plot(R_vector,diam_line,'color','k','LineWidth',1);
    hold on
    plot(Tangent_vector,Tangent_line,'color','k','LineWidth',1);
    hold on
    plot(org_x,org_y,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k') 
    hold on
     text(XX, YY, sprintf('   INF P = %f, %f', XX , YY))
     text(X_intercept_upper, yu , sprintf('  \\leftarrow  %f, %f', X_intercept_upper, yu))
      text(X_intercept_lower, ydn +0.005 , sprintf('      \\leftarrow  %f, %f', X_intercept_lower, ydn))
       text(0,0, sprintf(' \\leftarrow \\phi = %f',rad2deg(zetDa) ))
     title(['Drum free surface plot for (\alpha) = ',num2str(alpha)],'FontSize',25); 
   % RETURN CONDITIONS ARE SET TO FALSE SO SOLVE CAN INVOKE VPASOLVE SOLVE AUTOMATICALLY HENCE SOLVING PROBLEM NUMERICALLY. 
%% ===========================================================================================================================
%% ===========================================================================================================================
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
elseif alpha > 0.5 && alpha <=1
disp('                                     SOLVING INFLECTION POINT FOR ALPHA GREATER THAN 0.5')
    eqn1 = K*((-x/y)^3 - m_u*(-x/y)^2 + (-x/y) - m_u) - R^2 + x^2 + y^2;
    d = sqrt(x^2 + y^2);
    zeta = acos(d/R);
    eqn2 = 1- zeta/pi + (d*sin(zeta))/(pi*R) - alpha;
    [X,Y] = vpasolve([eqn1,eqn2],[x y],[-R 0;-R R]);
    R1 = matlabFunction(root1,'File','R1_optimized','Optimize',true);                         % optimising the real root
    XX = double(X);
    YY = double(Y);
     XX_check = length(XX);
    if  XX_check == 0
        disp('FAILED TO FIND SOLUTION TO INFLECTION POINT')
    end
    R_vector = linspace(0,XX,50);
    grad_diam = ((0-YY)/(0-XX));
    diam_line = grad_diam.*R_vector;
    Tangent_grad = -1/grad_diam;
    Y_int = YY - Tangent_grad*XX;
    Tangent_line = Tangent_grad.*Tangent_vector + Y_int;
    angle_of_repose = abs(rad2deg(atan(-1/grad_diam)));
    statictics = odeset('Stats','on');
    [xx,yy] = ode45(R1,[XX,R],YY);
    [xx1,yy1] = ode45(R1,[XX,-R],YY);
    domain = linspace(-R,R,1000);
    domain1 = domain;
    range = sqrt(R^2 - domain.^2);
    range1 = -sqrt(R^2 - domain1.^2);
    plot(domain,range,'color','b','LineWidth',2)
    axis equal
    hold on
    plot(domain,range1,'color','b','LineWidth',2)
    axis equal
    Tangent_line_intercept_upper = intersection_function(domain,range,Tangent_vector,Tangent_line,robust);
    Tangent_line_intercept_lower = intersection_function(domain,range1,Tangent_vector,Tangent_line,robust);
    [X_intercept_upper,yu] = intersection_function(domain,range,xx,yy,robust);
    [X_intercept_lower,ydn] = intersection_function(domain,range1,xx1,yy1,robust);
    Tangent_vector = linspace(Tangent_line_intercept_lower,Tangent_line_intercept_upper,200);   
    Tangent_line = Tangent_grad.*Tangent_vector + Y_int;
    hold on
    i = 1;
    while xx(i) < X_intercept_upper
        XX_vector(i) = xx(i);
        i = i +1;
    end
    YY_vector = yy(1:length(XX_vector));
    plot(XX_vector,YY_vector,'color','g','LineWidth',2);
    hold on
    j = 1;
    while xx1(j) > X_intercept_lower
        XX1_vector(j) = xx1(j);
        j = j +1;
    end
    YY1_vector = yy1(1:length(XX1_vector));
    R_vector1 = [0 Tangent_line_intercept_lower];
    len = sqrt(XX^2 + YY^2);
    R = double(R);
    zetDa =  acos(len/R);
    zetha = deg2rad(angle_of_repose) + pi/2 - zetDa;
    zetha = double(zetha);
    x2=(R*cos(zetha));
    y2=(R*sin(zetha));
    Y_tang_line = Tangent_grad*Tangent_line_intercept_lower + Y_int;
    R_vector1 = [0 Tangent_line_intercept_lower];
    Tang_line = ((0 - Y_tang_line)/(0 - Tangent_line_intercept_lower)).*R_vector1;
       
    plot([0 x2],[0 y2],'color','k','HandleVisibility','off')
    hold on
    plot(R_vector1,Tang_line,'color','k','HandleVisibility','off')
    hold on
    plot(R_vector,diam_line,'color','k','LineWidth',1);
    hold on
    plot(XX1_vector,YY1_vector,'color','r','LineWidth',2);
    hold on
    plot(R_vector,diam_line,'--','color','k','LineWidth',1);
    hold on
    plot(Tangent_vector,Tangent_line,'color','k','LineWidth',1);
    hold on
    plot(org_x,org_y,'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor','k') 
  text(XX, YY, sprintf('   INF P = %f, %f', XX , YY))
     text(X_intercept_upper, yu , sprintf('  \\leftarrow  %f, %f', X_intercept_upper, yu))
      text(X_intercept_lower, ydn +0.005 , sprintf('  \\leftarrow  %f, %f', X_intercept_lower, ydn))
       text(0,0, sprintf(' \\leftarrow \\phi = %f',rad2deg(zetDa) ))
    title(['Drum free surface plot for (\alpha) = ',num2str(alpha)],'FontSize',25); 
% RETURN CONDITIONS ARE SET TO FALSE SO SOLVE CAN INVOKE VPASOLVE SOLVE AUTOMATICALLY HENCE SOLVING PROBLEM NUMERICALLY.
%% ===============================================================================================================================
elseif alpha == 0.5
disp('                                          Trivial solution of fill fraction (0,0) is used for 50% fill fraction')
    XX = 0;
    YY = 0;
    statictics = odeset('Stats','on');
    clear figure()
    y = YY;
    x = XX;
    domain = linspace(-R,R,1000);
    domain1 = domain;
    range = sqrt(R^2 - domain.^2);
    range1 = -sqrt(R^2 - domain1.^2);
    plot(domain,range,'color','b','LineWidth',2)
    axis equal
    hold on
    plot(domain,range1,'color','b','LineWidth',2)
    axis equal
    R1 = matlabFunction(root1,'File','R1_optimized','Optimize',true);
    [xx,yy] = ode45(R1,[XX,R],YY);
    [xx1,yy1] = ode45(R1,[XX,-R],YY);
    X_intercept_upper = intersection_function(domain,range,xx,yy,robust);
    X_intercept_lower = intersection_function(domain,range1,xx1,yy1,robust);
     i = 1;
    while xx(i) < X_intercept_upper
        XX_vector(i) = xx(i);
        i = i +1;
    end
    YY_vector = yy(1:length(XX_vector));
    plot(XX_vector,YY_vector,'color','b','LineWidth',2);
    hold on
    j = 1;
    while xx1(j) > X_intercept_lower
        XX1_vector(j) = xx1(j);
        j = j +1;
    end
    YY1_vector = yy1(1:length(XX1_vector));
    plot(XX1_vector,YY1_vector,'color','r','LineWidth',2);
   
    plot(XX_vector,YY_vector,'color','g','LineWidth',2);
    
    plot(XX1_vector,YY1_vector,'color','r','LineWidth',2);
   title(['Drum free surface plot for (\alpha) = ',num2str(alpha)],'FontSize',25); 
    % Insert block code after optimisation has been established.
else 
    disp('                                      INVALID VALUE OF ALPHA, FILL FRACTION RANGES FROM 0 T0 UNITY')
end
%% ================================================================================================================================

figure('Name','mohr columb effect')
%% THIS THE MAIN SCRIPT IN WHICH THE ENTIRE CODE IS EMBBEDED. ALL SUBSEQUENT FUNCTIONS ARE INVOKED IN THIS SCRIPT; KIND REGARDS, GROUP 5
%% =======================================================================================================================================
clc;
clear
legend('-DyanamicLegend');
clc
%%                                                  TO GOD BE THE GLORY
%% GROUP MEMBERS
%% 1.  Phelelani Eshmael Mamba (217058111)
%% 2.  Tashmira Ramjan         (217003070)
%% 3.  Taryn Gore              (217003065)
%% 4.  Prishalan Moodley       (217002787)
%% ========================================================================================================================================
%%                                                 THE DECLARATION OF VARIABLES
%% VARIABLE                                        VARIABLE DESCRIPTION

for k = 1:5
colours = 'rgbmc';
P = [0.2 0.4 0.6 0.8 1];
used_color = colours(k);
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
p0 = rho*phi*g*h0*cos(theta);                      % Pressure at depth h0
alpha = 0.4;                                       % fill fraction  
m_u = P(k);                                        % Mohr columb friction
K = (2*p0^3)/(3*eta * omega * rho^2*g^2);          % Bulk reighnological parameter
Tangent_vector = linspace(-R,R,100);               % Pre-defining domain for tangent line
robust = 0;                                        % Intersection function usage 
%% =========================================================================================================================
ode = K*((dy)^3 - m_u*dy^2 +dy - m_u) == R^2 - y^2 - x^2;
derivatives = solve(ode,dy,'MaxDegree',3);
root1 = derivatives(1);
root2 = derivatives(2);
root3 = derivatives(3);

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
    plot(xx,yy,used_color,'DisplayName',num2str(P(k)),'HandleVisibility','off'); 
    hold on
    plot(xx1,yy1,used_color,'DisplayName',num2str(P(k))); 
    plot(domain,range,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
    plot(domain,range1,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
    title(['Mohr Coulomb friction parametric study for (\alpha) = ',num2str(alpha)],'FontSize',25); 
    legend('show')  
%% ===========================================================================================================================
%% ===========================================================================================================================
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
end
 
figure('Name','Radius parametric study')
%% THIS THE MAIN SCRIPT IN WHICH THE ENTIRE CODE IS EMBBEDED. ALL SUBSEQUENT FUNCTIONS ARE INVOKED IN THIS SCRIPT; KIND REGARDS, GROUP 5
%% =======================================================================================================================================
clc;
clear;
legend('-DyanamicLegend');
clc
%%                                                  TO GOD BE THE GLORY
%% GROUP MEMBERS
%% 1.  Phelelani Eshmael Mamba (217058111)
%% 2.  Tashmira Ramjan         (217003070)
%% 3.  Taryn Gore              (217003065)
%% 4.  Prishalan Moodley       (217002787)
%% ========================================================================================================================================
%%                                                 THE DECLARATION OF VARIABLES
%% VARIABLE                                        VARIABLE DESCRIPTION
for k = 1:5                                        % Intiating for loop
colours = 'rgbmc';                                 % Defining colour of plot
P = [0.0476 0.0952 0.1428 0.1904 0.238];           % Generic vector P
used_color = colours(k);                           % Defining colour of plot
syms m_u K x y R d dy                              % SYMBOLIC DECLARATION OF VARIABLES
rho = 2500;                                        % density of particles
g   = 9.81;                                        % gravitational acceleration
R   = P(k);                                       % Radius of mill
P_crit = 0.6;                                      % 60 % of critical speed
crit_speed = 42.3/(sqrt(2*R));                     % Mill critical speed
speed = P_crit * crit_speed;                       % Mill actual speed in RPM
omega = (speed * 2*pi)/60;                         % omega od drum in rad/sec
eta  = 6e-3;                                       % Granular viscosity
d_p = 3e-3;                                        % Uniform particle diameter
h0 = 1.3*d_p;                                      % Depth of flowing layer
phi = 0.58;                                        % Average solids fraction
theta = deg2rad(40);                               % Angle of repose measured in the anticlockwise direction (Hence 180-40 notation)
p0 = rho*phi*g*h0*cos(theta);                      % Pressure at depth h0
alpha = 0.4;                                       % fill fraction  
m_u = 0.35;                                        % Mohr columb friction
K = (2*p0^3)/(3*eta * omega * rho^2*g^2);          % Bulk reighnological parameter
Tangent_vector = linspace(-R,R,100);               % Pre-defining domain for tangent line
robust = 0;                                        % Intersection function usage 
%% =========================================================================================================================
ode = K*((dy)^3 - m_u*dy^2 +dy - m_u) == R^2 - y^2 - x^2;
derivatives = solve(ode,dy,'MaxDegree',3);
root1 = derivatives(1);
root2 = derivatives(2);
root3 = derivatives(3);

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
    plot(xx,yy,used_color,'DisplayName',num2str(P(k)),'HandleVisibility','off'); 
    hold on
    plot(xx1,yy1,used_color,'DisplayName',num2str(P(k))); 
    plot(domain,range,'--','color','k','HandleVisibility','off')
    axis equal
    hold on
    plot(domain,range1,'--','color','k','HandleVisibility','off')
    axis equal
    hold on
    title(['Radius of drum parametric study for (\alpha) = ',num2str(alpha)],'FontSize',25); 
    legend('show')  
%% ===========================================================================================================================
%% ===========================================================================================================================
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
end


figure('Name','Viscosity parametric study')
%% THIS THE MAIN SCRIPT IN WHICH THE ENTIRE CODE IS EMBBEDED. ALL SUBSEQUENT FUNCTIONS ARE INVOKED IN THIS SCRIPT; KIND REGARDS, GROUP 5
%% =======================================================================================================================================
clc;
clear;
legend('-DyanamicLegend');
clc
%%                                                  TO GOD BE THE GLORY
%% GROUP MEMBERS
%% 1.  Phelelani Eshmael Mamba (217058111)
%% 2.  Tashmira Ramjan         (217003070)
%% 3.  Taryn Gore              (217003065)
%% 4.  Prishalan Moodley       (217002787)
%% ========================================================================================================================================
%%                                                 THE DECLARATION OF VARIABLES
%% VARIABLE                                        VARIABLE DESCRIPTION
for k = 1:5                                        % Intiating for loop
colours = 'rgbmc';                                 % Defining colour of plot
P = [0.2 0.4 0.6 0.8 1]*1e-2;                      % Generic vector P
used_color = colours(k);                           % Defining colour of plot
syms m_u K x y R d dy                              % SYMBOLIC DECLARATION OF VARIABLES
rho = 2500;                                        % density of particles
g   = 9.81;                                        % gravitational acceleration
R   = 0.238;                                       % Radius of mill
P_crit = 0.6;                                      % 60 % of critical speed
crit_speed = 42.3/(sqrt(2*R));                     % Mill critical speed
speed = P_crit * crit_speed;                       % Mill actual speed in RPM
omega = (speed * 2*pi)/60;                         % omega od drum in rad/sec
eta  = P(k);                                       % Granular viscosity
d_p = 3e-3;                                        % Uniform particle diameter
h0 = 1.3*d_p;                                      % Depth of flowing layer
phi = 0.58;                                        % Average solids fraction
theta = deg2rad(40);                               % Angle of repose measured in the anticlockwise direction (Hence 180-40 notation)
p0 = rho*phi*g*h0*cos(theta);                      % Pressure at depth h0
alpha = 0.4;                                       % fill fraction  
m_u = 0.35;                                        % Mohr columb friction
K = (2*p0^3)/(3*eta * omega * rho^2*g^2);          % Bulk reighnological parameter
Tangent_vector = linspace(-R,R,100);               % Pre-defining domain for tangent line
robust = 0;                                        % Intersection function usage 
%% =========================================================================================================================
ode = K*((dy)^3 - m_u*dy^2 +dy - m_u) == R^2 - y^2 - x^2;
derivatives = solve(ode,dy,'MaxDegree',3);
root1 = derivatives(1);
root2 = derivatives(2);
root3 = derivatives(3);

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
    plot(xx,yy,used_color,'DisplayName',num2str(P(k)),'HandleVisibility','off'); 
    hold on
    plot(xx1,yy1,used_color,'DisplayName',num2str(P(k))); 
    plot(domain,range,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
    plot(domain,range1,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
    title(['Viscosity parametric study for (\alpha) = ',num2str(alpha)],'FontSize',25); 
    legend('show')  
%% ===========================================================================================================================
%% ===========================================================================================================================
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
end
figure('Name','Omega parametric study')
%% THIS THE MAIN SCRIPT IN WHICH THE ENTIRE CODE IS EMBBEDED. ALL SUBSEQUENT FUNCTIONS ARE INVOKED IN THIS SCRIPT; KIND REGARDS, GROUP 5
%% =======================================================================================================================================
clc;
clear
legend('-DyanamicLegend');
clc
%%                                                  TO GOD BE THE GLORY
%% GROUP MEMBERS
%% 1.  Phelelani Eshmael Mamba (217058111)
%% 2.  Tashmira Ramjan         (217003070)
%% 3.  Taryn Gore              (217003065)
%% 4.  Prishalan Moodley       (217002787)
%% ========================================================================================================================================
%%                                                 THE DECLARATION OF VARIABLES
%% VARIABLE                                        VARIABLE DESCRIPTION

for k = 1:5
colours = 'rgbmc';
P = [0.2 0.4 0.6 0.8 1];
used_color = colours(k);
syms m_u K x y R d dy                              % SYMBOLIC DECLARATION OF VARIABLES
rho = 2500;                                        % density of particles
g   = 9.81;                                        % gravitational acceleration
R   = 0.238;                                       % Radius of mill
P_crit = P(k);                                      % 60 % of critical speed
crit_speed = 42.3/(sqrt(2*R));                     % Mill critical speed
speed = P_crit * crit_speed;                       % Mill actual speed in RPM
omega = (speed * 2*pi)/60;                         % omega od drum in rad/sec
eta  = 6e-3;                                       % Granular viscosity
d_p = 3e-3;                                        % Uniform particle diameter
h0 = 1.3*d_p;                                      % Depth of flowing layer
phi = 0.58;                                        % Average solids fraction
theta = deg2rad(40);                               % Angle of repose measured in the anticlockwise direction (Hence 180-40 notation)
p0 = rho*phi*g*h0*cos(theta);                      % Pressure at depth h0
alpha = 0.4;                                       % fill fraction  
m_u = 0.35;                                        % Mohr columb friction
K = (2*p0^3)/(3*eta * omega * rho^2*g^2);          % Bulk reighnological parameter
Tangent_vector = linspace(-R,R,100);               % Pre-defining domain for tangent line
robust = 0;                                        % Intersection function usage 
%% =========================================================================================================================
ode = K*((dy)^3 - m_u*dy^2 +dy - m_u) == R^2 - y^2 - x^2;
derivatives = solve(ode,dy,'MaxDegree',3);
root1 = derivatives(1);
root2 = derivatives(2);
root3 = derivatives(3);

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
    plot(xx,yy,used_color,'DisplayName',num2str(P(k)),'HandleVisibility','off'); 
    hold on
    plot(xx1,yy1,used_color,'DisplayName',num2str(P(k)*100)); 
    plot(domain,range,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
    plot(domain,range1,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
    title(['Omega parametric study for (\alpha) = ',num2str(alpha)],'FontSize',25); 
    legend('show')  
%% ===========================================================================================================================
%% ===========================================================================================================================
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
end

figure('Name','effect of mu below inflection point')
%% THIS THE MAIN SCRIPT IN WHICH THE ENTIRE CODE IS EMBBEDED. ALL SUBSEQUENT FUNCTIONS ARE INVOKED IN THIS SCRIPT; KIND REGARDS, GROUP 5
%% =======================================================================================================================================
clc;
clear
legend('-DyanamicLegend');
clc
%%                                                  TO GOD BE THE GLORY
%% GROUP MEMBERS
%% 1.  Phelelani Eshmael Mamba (217058111)
%% 2.  Tashmira Ramjan         (217003070)
%% 3.  Taryn Gore              (217003065)
%% 4.  Prishalan Moodley       (217002787)
%% ========================================================================================================================================
%%                                                 THE DECLARATION OF VARIABLES
%% VARIABLE                                        VARIABLE DESCRIPTION                                 
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
p0 = rho*phi*g*h0*cos(theta);                      % Pressure at depth h0
alpha = 0.4;                                       % fill fraction  
m_u = 0.35;                                        % Mohr columb friction
K = (2*p0^3)/(3*eta * omega * rho^2*g^2);          % Bulk reighnological parameter
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
    plot(xx,yy,'r','DisplayName',num2str(m_u),'HandleVisibility','off'); 
    hold on
    plot(xx1,yy1,'r','DisplayName',num2str(m_u)); 
    plot(domain,range,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
    plot(domain,range1,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
%% ===========================================================================================================================
%% ===========================================================================================================================
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%% REDIFINING VARIABLES                            VARIABLE DESCRIPTION
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
theta = deg2rad(40);                               %
p0 = rho*phi*g*h0*cos(theta);                      % Pressure at depth h0
alpha = 0.4;                                       % fill fraction
m_u = 0.35;                                        % Mohr columb friction
K = (2*p0^3)/(3*eta * omega * rho^2*g^2);          % Bulk reighnological parameter
robust = 0;                                        % Intersection function usage
%% ============================================================================================================================
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ode = K*((dy)^3 - m_u*dy^2 +dy - m_u) == R^2 - y^2 - x^2;
derivatives = solve(ode,dy,'MaxDegree',3);
root1 = derivatives(1);
% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%% *****************************************************************************************************************************
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
    X_intercept_upper = intersection_function(domain,range,xx,yy,robust);
    [xx,yy] = ode45(R1,[XX,X_intercept_upper(1)],YY);
    if length(X_intercept_upper) > 1
        disp('                                                THE DRUM IS CENTRIFUGING')
    end
    plot(xx,yy,'b','DisplayName',num2str(m_u),'HandleVisibility','off'); 
    hold on
    % plot(xx1,yy1,used_color,'DisplayName',num2str(P(k))); 
    plot(domain,range,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
    plot(domain,range1,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
%% REDIFINING VARIABLES                            VARIABLE DESCRIPTION
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
theta = deg2rad(40);                               %
p0 = rho*phi*g*h0*cos(theta);                      % Pressure at depth h0
alpha = 0.4;                                       % fill fraction
m_u = 0.35/1000;                                        % Mohr columb friction
K = (2*p0^3)/(3*eta * omega * rho^2*g^2);          % Bulk reighnological parameter
robust = 0;                                        % Intersection function usage
%% ============================================================================================================================
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ode = K*((dy)^3 - m_u*dy^2 +dy - m_u) == R^2 - y^2 - x^2;
derivatives = solve(ode,dy,'MaxDegree',3);
root1 = derivatives(1);
root2 = derivatives(2);
root3 = derivatives(3);
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
    % plot(xx,yy,used_color,'DisplayName',num2str(P(k)),'HandleVisibility','off'); 
    hold on
    plot(xx1,yy1,'b','DisplayName',num2str(m_u)); 
    plot(domain,range,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
    plot(domain,range1,'color','k','LineWidth',2,'HandleVisibility','off')
    axis equal
    hold on
    title(['Effect of changing (\mu) below inflection point study for (\alpha) = ',num2str(alpha)],'FontSize',25); 
    legend('show')  
%% ===========================================================================================================================
%% ===========================================================================================================================
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
