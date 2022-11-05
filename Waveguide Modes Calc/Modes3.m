
h = 20*10E-6;  % core size 

%RI indices
n_core = 1.477;
n_clad = 1.475; 
ns = 1.475; 

lambda = 0.647;  
k = 2*pi/lambda;

rootBInt = [0, 1];
a = 0
gamma = (ns^2 - n_clad^2)/(n_core^2 - ns^2); %parametr nierówności (XD nie wiem już)

V = (k*h)/2*sqrt(n_core^2-ns^2);
v = k*h*sqrt(n_core^2-n_clad^2)

b = (0:0.01:1);

EigenEquation = @(b,m) m*pi+atan(sqrt(b./(1-b)))+atan(sqrt((b+a)./(1-b))) - v.*sqrt(1-b)
%by = atan(sqrt(bx./(1 - bx))) +atan(sqrt((bx + gamma)./(1 - bx))) - 2.*vT.*sqrt(1 - bx); 
modes = 5
sol = []
leg={}
figure('Color','w')
for m = 0:modes
    sol = [sol;EigenEquation(b,m)];
    plot(b,sol(m+1,:));
    xlabel("b[normalised]")
    ylabel("EigenEquation values[a.u.]")
    leg = [leg;"Mode"+m]
    hold on
    y_line = m*pi
    plot([0 max(b)],[y_line y_line],'Color','red','LineStyle','--', HandleVisibility='off')
end
hold off
legend(leg)



