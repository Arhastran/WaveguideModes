%Adam Ignaciuk 
%Waveguide mode simulations 
clear
h = 15*1E-6; % core size width(x)
h2 = 15*1E-6; %core size (y)

dist = 0.01*1E-6;

%RI indices
n_e = 1.739;
n_o = 1.522;
n_core = 2.038; %1.639
n_clad = 1.41; 
ns = n_clad; %RI below surface, I assume
modes = 1;
lambda = 647*1E-9%589*1E-9;  
k = 2*pi/lambda;


a = (ns^2 - n_clad^2)/(n_core^2 - ns^2); %parametr nierówności 

%------------------------------------------------------------------------------------
% X mode propagation 
%------------------------------------------------------------------------------------

v = k*h*sqrt(n_core^2-n_clad^2);

b = (0:0.000001:1);

EigenEquation = @(b,m) m*pi+atan(sqrt(b./(1-b)))+atan(sqrt((b+a)./(1-b)));
Eigen2 = @(b) v.*sqrt(1-b);
%by = atan(sqrt(bx./(1 - bx))) +atan(sqrt((bx + gamma)./(1 - bx))) - 2.*vT.*sqrt(1 - bx); 
sol = [];
leg={};
figure('Color','w');
sol2 = Eigen2(b);
plot(b, sol2);
leg = [leg;"Left side of equation"];
hold on
for m = 1:modes
    sol = [sol;EigenEquation(b,m-1)];
    plot(b,sol(m,:));
    xlabel("b[normalised]");
    ylabel("EigenEquation values[a.u.]");
    leg = [leg;"Right side of equation"];
    %y_line = m*pi
    root = abs(EigenEquation(b,m-1)-Eigen2(b));
    index = find(root==min(root));
    intRoot = [0 1];
    Intersections = fzero(@(b) EigenEquation(b,m-1)-Eigen2(b),intRoot);
    Xintersect = Intersections;
    Yintersect = EigenEquation(Xintersect,m-1);
    plot(Xintersect,Yintersect, 'sg', 'MarkerSize',10);
    leg = [leg;"Intersection (err=0.0001)"];
    %print("Value for normalised b for "+m+"mode, :"+Xintersect)
    %plot([0 max(b)],[y_line y_line],'Color','red','LineStyle','--', HandleVisibility='off')
end
hold off
legend(leg);

%Effective index 
nTE = sqrt(Xintersect*(n_core^2 - ns^2) + ns^2);
%difference = nTE-n_core;

%some values (im not yet sure what they are about)
betaTE = nTE*k;
u = v*sqrt(1 - Xintersect);
w =  v*sqrt(Xintersect);
wP = v*sqrt(Xintersect + a);

% transverse propagation constants
ki = u/(h/2);       
sigma = w/(h/2);     
xi = wP/(h/2);    
% phase
phi = 0.5*atan(w/u) - 0.5*atan(wP/u);

% power confinment 
powerConfTE = (1 + (sin(u + phi))^2/(2*w) + (sin(u - phi))^2/(2*wP)) / (1 + 1/(2*w) + 1/(2*wP) );
   

% effective width
heTE = h + 1/sigma + 1/xi;
%diffWidth = h/heTE;

% numerical aperture
NA0Y = sqrt(n_core - n_clad);
NAsY = sqrt(n_core - ns);

%cutoff Wavelength
lambdaC0T = 2*pi/v*h/2*NA0Y;
lambdaCsT = 2*pi/v*h/2*NAsY;

x = (- h : dist : h);

    for xIdx = 1: length(x)
        % Electric field distribution: TE mode profile in the plane perpendicular
        % to the propagation direction
        ex(xIdx) = Efield(x(xIdx), h, ki, sigma, phi, xi);

    end

   figure("Color",'w');
   plot(x, ex, 'k');
   hold on
   plot( repmat(-(h/2),1, length(ex)), ex, 'r--');
   plot( repmat(+(h/2),1, length(ex)), ex, 'r--');
   xlabel('x[m]');
   ylabel('Efield/Emax');
   title('TE0 mode profile along X');
   hold off

%------------------------------------------------------------------------------------
% Y mode propagation 
%------------------------------------------------------------------------------------

v = k*h2*sqrt(n_core^2-n_clad^2);

b = (0:0.000001:1);

EigenEquation = @(b,m) m*pi+atan(sqrt(b./(1-b)))+atan(sqrt((b+a)./(1-b)));
Eigen2 = @(b) v.*sqrt(1-b);
%by = atan(sqrt(bx./(1 - bx))) +atan(sqrt((bx + gamma)./(1 - bx))) - 2.*vT.*sqrt(1 - bx); 
sol = [];
leg={};
figure('Color','w');
sol2 = Eigen2(b);
plot(b, sol2);
leg = [leg;"Left side of equation"];
hold on
for m = 1:modes
    sol = [sol;EigenEquation(b,m-1)];
    plot(b,sol(m,:));
    xlabel("b[normalised]");
    ylabel("EigenEquation values[a.u.]");
    leg = [leg;"Right side of equation"];
    %y_line = m*pi
    root = abs(EigenEquation(b,m-1)-Eigen2(b));
    index = find(root==min(root));
    Intersections = fzero(@(b) EigenEquation(b,m-1)-Eigen2(b),intRoot);
    intRoot = [0 1];
    Xintersect = Intersections;
    Yintersect = EigenEquation(Xintersect,m-1);
    plot(Xintersect,Yintersect, 'sg', 'MarkerSize',10);
    leg = [leg;"Intersection (err=0.0001)"];
    %print("Value for normalised b for "+m+"mode, :"+Xintersect)
    %plot([0 max(b)],[y_line y_line],'Color','red','LineStyle','--', HandleVisibility='off')
end
hold off
legend(leg);

%Effective index 
nTE = sqrt(Xintersect*(n_core^2 - ns^2) + ns^2);
%difference = nTE-n_core;

%some values (im not yet sure what they are about)
betaTE = nTE*k;
u = v*sqrt(1 - Xintersect);
w =  v*sqrt(Xintersect);
wP = v*sqrt(Xintersect + a);

% transverse propagation constants
ki = u/(h2/2);       
sigma = w/(h2/2);     
xi = wP/(h2/2);    
% phase
phi = 0.5*atan(w/u) - 0.5*atan(wP/u);

% power confinment 
powerConfTE = (1 + (sin(u + phi))^2/(2*w) + (sin(u - phi))^2/(2*wP)) / (1 + 1/(2*w) + 1/(2*wP) );
   

% effective width
heTE = h2 + 1/sigma + 1/xi;
%diffWidth = h/heTE;

% numerical aperture
NA0Y = sqrt(n_core - n_clad);
NAsY = sqrt(n_core - ns);

%cutoff Wavelength
lambdaC0T = 2*pi/v*h2/2*NA0Y;
lambdaCsT = 2*pi/v*h2/2*NAsY;

y = (- h2 : dist : h2);

    for yIdy = 1: length(y)
        % Electric field distribution: TE mode profile in the plane perpendicular
        % to the propagation direction
        ey(yIdy) = Efield(y(yIdy), h2, ki, sigma, phi, xi);

    end

   figure("Color",'w');
   plot(y, ey, 'k');
   hold on
   plot( repmat(-(h2/2),1, length(ey)), ey, 'r--');
   plot( repmat(+(h2/2),1, length(ey)), ey, 'r--');
   xlabel('y[m]');
   ylabel('Efield/Emax');
   title('TE0 mode profile along Y');
   hold off


if length(ex)>length(ey)

    N_max = max([length(ey), length(ex)]);
    n_v = linspace(0,1,N_max);
    new_X = interp1(linspace(0,1,length(ey)), ey, n_v);
    new_Y = interp1(linspace(0,1,length(ex)), ex, n_v);

    ez = new_X.*new_Y;

    N_max2 = max([length(y), length(x)]);
    n_v2 = linspace(0,1,N_max2);
    new_X2 = interp1(linspace(0,1,length(y)), y, n_v2);
    new_Y2 = interp1(linspace(0,1,length(x)), x, n_v2);


elseif length(ex)<length(ey)

    N_max = max([length(ex), length(ey)]);
    n_v = linspace(0,1,N_max);
    new_X = interp1(linspace(0,1,length(ex)), ex, n_v);
    new_Y = interp1(linspace(0,1,length(ey)), ey, n_v);

    ez = new_X.*new_Y;

    N_max2 = max([length(x), length(y)]);
    n_v2 = linspace(0,1,N_max2);
    new_X2 = interp1(linspace(0,1,length(x)), x, n_v2);
    new_Y2 = interp1(linspace(0,1,length(y)), y, n_v2);
else
    new_X=ex;
    new_Y=ey;
    ez = ex.*ey;
    new_X2 = x;
    new_Y2 = y;
    
end

        
figure("Color",'w');
Z = ex( 1:min(length(ex), length(ey)) )'*ey( 1:min(length(ex), length(ey)) );
% surf(Z);
image(flipud(Z),'CDataMapping','scaled'), colormap('jet')
axis equal
xlabel("x[um]");
ylabel("y[um]");
colorbar
%figure, mesh(flipud(Z))

xticklabels = ((-h)*1E6 : 5 : (h)*1E6);
set(gca, 'XTickLabel', xticklabels);
yticklabels = ((-h2)*1E6 : 5 : (h2)*1E6);
set(gca, 'YTickLabel', yticklabels);


   
