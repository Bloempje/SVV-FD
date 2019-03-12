
flightdata = open('flightdata.mat')
%% series 1

% fixed input
W_op = 5000*g % N

% measured
t = [18*60+3, 20*60+4, 24*60+1, 27*60+26, 29*60+37, 37*60+0] % sec
h = [6990, 7000, 8010, 8000, 8010, 8980] % ft
V = [241, 219, 192, 161, 137, 119] % kts
TAT = [1.2, -1.2, -2.2, -6.5, -6.2, -7.0] % deg C
FFl = [720, 603, 492, 464, 387, 373] % lbs/hr
FFr = [756, 643, 533, 509, 426, 409] % lbs/hr
Wf = [370, 407, 481, 530, 559, 670] % lbs
alpha = [1.6, 2.4, 3.5, 5.3, 7.7, 10.2] % deg

for i=1:size(h,2)
    
    CL(i) = 2*(W_op-Wf(i))*0.45359237 / (rho0*(V(i)*0.5144)^2*S)
    CD(i) = 2*T / (rho0*(V(i)*0.5144)^2*S)
end

%% CD = CD0 + k* CL^2

% x vector
%CL = [0, 0.2, 0.4, 0.5, 0.8, 1.2];
% y vector
%CD = [0; 0.1; 0.1; 0.2; 0.5; 0.9]; %column vector, CD_flight

m = size(CL, 2); %nr of data points

X = zeros(m,2); % create X matrix with the correct dimentions
% fill X with the correct values
for i = 1:m
    X(i,1) = 1
    X(i,2) = CL(i)^2
end


% closed for optimal least squares solution
copt = ((transpose(X)*X)^(-1))*transpose(X) * CD

% cvx optimization
cvx_begin quiet
variable c(2); % variable is a vector with 4 elements
minimize norm((X*c-CD),2); %function to minimize. Note that we are not taking the
% norm squared. A norm is always positive so minimizing this is the same
% problem as optimizing the norm squared.
cvx_end
two_norm_regression = c

cvx_begin quiet
variable c(2); % variable is a vector with 4 elements
minimize norm((X*c-CD),1); %function to minimize. Note that we are not taking the
% norm squared. A norm is always positive so minimizing this is the same
% problem as optimizing the norm squared.
cvx_end
one_norm_regression = c

cvx_begin quiet
variable c(2); % variable is a vector with 4 elements
minimize norm((X*c-CD),inf); %function to minimize. Note that we are not taking the
% norm squared. A norm is always positive so minimizing this is the same
% problem as optimizing the norm squared.
cvx_end
inf_norm_regression = c


% coefficients optimized using the infinity-, 1- and 2-norm
cinf = inf_norm_regression; 
c1 = one_norm_regression;
c2 = two_norm_regression;


% plot data points
plot(CL,CD, 'o')
%xlim([x(1) - 1, x(m) + 1]) % set x-limits from 0 to 8
hold on % new plots will join in this window instead of overwriting 
x = 0 :0.01: CL(6)*1.2; % set a range of x values for plotting
% calculate range of y values for plotting of regression curves
yinf = cinf(1) + cinf(2)*x.^2;
y1 = c1(1) + c1(2)*x.^2;
y2 = c2(1) + c2(2)*x.^2;
ygiven = CD0 + x.^2 /(pi*A*e);
plot(x,y2) % plot regression curves
plot(x, y1);
plot(x, yinf);
plot(x, ygiven)
% add a legend
legend({'data points','2-norm optimization', '1-norm optimization', ...
    'infinity norm optimization'},'Location','northwest')