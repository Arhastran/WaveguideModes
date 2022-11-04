% Adam Ignaciuk 
% script for modes in waveguide slab 

%RI indices
n1 = 1.4513;
n2 = 1.4440;
neff_range = [n2 n1]; %nwm
h = 100; %width of waveguide
a = 0; %assymetry parameter
c=physconst('lightspeed'); %lightspeed
lambdas = (600:0.01:1000); %waves 
k = 2*pi./(lambdas); %wavenumbers 
lambda = lambdas(1); 
k0 = 2*pi/(lambda); %single wavenumebr 
V = h * k0*sqrt(n1^2-n2^2); %normalised frequency 

alpha = (pi:0.01*pi:2*pi);
N = n1*sin(alpha);
b = (N.^2-n1^2)./(n1^2-n2^2);

leftSide = @(b) V*sqrt(1-b);
left = leftSide(b);
rightSide = @(b,m) pi*m+atan(sqrt((b)/(1-b)))+atan(sqrt((b+a)/(1-b)));


plot(b,left)
xlabel("b[normalised]")
ylabel("crosssection")

rigthSol = []

for n=1:10
    rightSol(n) = abs(rightSide(b,n))
end





