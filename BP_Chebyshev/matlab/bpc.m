%Band pass filter - Chebyshev
%Project creator :Ioannis Kessopoulos 9271

%Main specifications

a = [ 9 2 7 1];
m = 2;

%Freequencies

f_0 = 0.65*1000;
f_1 = 400 + 25*a(3);
f_2 = (f_0^2)/ f_1;
D = 2.3*((f_0^2 - f_1^2)/f_1);
f_3 = ((-D + sqrt((D)^2 + 4*(f_0^2)))/2);
f_4 = (f_0^2)/f_3;

a_min = 27.5 + a(4);
a_max = 0.5 + (a(3)-5)/10;

w_0 = 2*pi*f_0; 
w_1 = 2*pi*f_1;
w_2 = 2*pi*f_2; 
w_3 = 2*pi*f_3;
w_4 = 2*pi*f_4;

bw = w_2 - w_1;
W_p = 1;
W_s = (w_4 - w_3)/(w_2 - w_1);

n = acosh(sqrt(((10^(a_min/10) -1))/ (10^(a_max/10)- 1)))/acosh(W_s);
n=ceil(n);


%Starting the proscess in order to calculate the poles

epsilon = sqrt(10^(a_max/10) -1);
W_hp = cosh((1/n)*acosh(1/epsilon));

a = (1/n)*(asinh(1/epsilon));

%Butterworth angles-Definitionn=4 +- 22.5 +-67.5

s_1 = -sinh(a)*cos(deg2rad(22.5));
s_2 = -sinh(a)*cos(deg2rad(67.5));

W_1 = cosh(a)*sin(deg2rad(22.5));
W_2 = cosh(a)*sin(deg2rad(67.5));

%Poles-Calculation

%First pair

GeffeComplexPoles(-s_1 , W_1 , w_0, bw, 1);

%Second pair

GeffeComplexPoles(-s_2 , W_2 , w_0, bw, 2);

%Kathorismos timwn monadwn-Strathgikes sxediashs

%Monada 1 Q<5

w_o = Geffe_first_pair.w_01;
Q = Geffe_first_pair.Q;
DelyiannisFried1( w_o, Q, 1)

%Monada 2 Q<5

w_o = Geffe_first_pair.w_02;
Q = Geffe_first_pair.Q;
DelyiannisFried1( w_o, Q, 2)

%Monada 3 Q>5

w_o = Geffe_second_pair.w_01;
Q = Geffe_second_pair.Q;
DelyiannisFriedEnchancement( w_o, Q, 4, 1)

%Monada 4 Q>5

w_o = Geffe_second_pair.w_02;
Q = Geffe_second_pair.Q;
DelyiannisFriedEnchancement( w_o, Q, 4, 2)

%Transfer functions- Definitions -Calculations

tf_1 = tf([-Monada1.H*(Geffe_first_pair.w_01/Geffe_first_pair.Q) 0], [1 (Geffe_first_pair.w_01/Geffe_first_pair.Q) Geffe_first_pair.w_01^2])
tf_2 = tf([-Monada2.H*(Geffe_first_pair.w_02/Geffe_first_pair.Q) 0], [1 (Geffe_first_pair.w_02/Geffe_first_pair.Q) Geffe_first_pair.w_02^2])
tf_3 = tf([-Monada3.H*(Geffe_second_pair.w_01/Geffe_second_pair.Q) 0], [1 (Geffe_second_pair.w_01/Geffe_second_pair.Q) Geffe_second_pair.w_01^2])
tf_4 = tf([-Monada4.H*(Geffe_second_pair.w_02/Geffe_second_pair.Q) 0], [1 (Geffe_second_pair.w_02/Geffe_second_pair.Q) Geffe_second_pair.w_02^2])

%Transfer functions- Plots

plot_transfer_function( tf_1, [ f_3 f_1 f_0 f_2 f_4]);
plot_transfer_function( tf_2, [ f_3 f_1 f_0 f_2 f_4]);
plot_transfer_function( tf_3, [ f_3 f_1 f_0 f_2 f_4]);
plot_transfer_function( tf_4, [ f_3 f_1 f_0 f_2 f_4]);

%Gain
%Asked Gain
Gain_dB = 0;
Gain_num = 10^(Gain_dB/20);

%Recuired calculations for gain correction

w_o = Geffe_first_pair.w_01;
Q = Geffe_first_pair.Q;
k1 = (2*Q*w_o*w_0)/(sqrt((w_o^2 - w_0^2)^2 + (w_o*w_0/Q)^2));
 
w_o = Geffe_first_pair.w_02;
Q = Geffe_first_pair.Q;
k2 = (2*Q*w_o*w_0)/(sqrt((w_o^2 - w_0^2)^2 + (w_o*w_0/Q)^2));

w_o = Geffe_second_pair.w_01;
Q = Geffe_second_pair.Q;
k3 = ( ((Monada3.k*w_0)/((Monada3.k-1)*Monada3.R_1*Monada3.C_1)) / (sqrt((w_o^2 - w_0^2)^2 + (w_o*w_0/Q)^2)) );

w_o = Geffe_second_pair.w_02;
Q = Geffe_second_pair.Q;
k4 = ( ((Monada4.k*w_0)/((Monada4.k-1)*Monada4.R_1*Monada4.C_1)) / (sqrt((w_o^2 - w_0^2)^2 + (w_o*w_0/Q)^2)) );
a = Gain_num/(k1*k2*k3*k4);

%Voltage Dividers

a_2 = nthroot((a*(k3^2))/(k1^2),4);
a_1 = (k1/k3)*a_2;

%Monada 1
R_1A = Monada1.R_1/a_1;  
R_1B = Monada1.R_1/(1-a_1); 

%Monada 2
R_2A = Monada2.R_1/a_1;
R_2B = Monada2.R_1/(1-a_1);

%Monada 3
R_3A = Monada3.R_1/a_2;
R_3B = Monada3.R_1/(1-a_2);

%Monada 4
R_4A = Monada4.R_1/a_2;
R_4B = Monada4.R_1/(1-a_2);

%Final Results -Total transfer functions representations

Total_tf_1 = series(tf_1,tf_2);
Total_tf_2 = series(tf_3,tf_4);
Total_tf_3 = series(Total_tf_1, Total_tf_2);
Total_tf = series(a, Total_tf_3);

plot_transfer_function( Total_tf, [ f_3 f_1 f_0 f_2 f_4]);


%Inverted total transfer function representation

Total_inv_tf = inv(Total_tf);
plot_transfer_function(Total_inv_tf, [ f_3 f_1 f_0 f_2 f_4]);

%Fourier Analysis

f11 = (w_0 - (w_0-w_1)/2) / (2*pi);
f12 = (w_0 + (w_0+w_1)/3) / (2*pi);
f13 = 0.4*w_3 / (2*pi);
f14 = (2.5*w_4) / (2*pi);
f15 = (3*w_4) / (2*pi) ;
fs= 20000;
T=0.002;
dt=1/fs;
t=0:dt:(T);

u1 = cos(2*pi*f11*t)+0.8* cos(2*pi* f12*t)+0.8*cos(2*pi*0.8* f13*t)+0.6*cos(2*pi*f14*t)+0.5*cos(2*pi*f15*t);
figure
plot(u1)

N=T/dt;
figure
lsim(Total_tf,u1,t)
xt=lsim(Total_tf,u1,t);
figure
plot(t,xt)
n=2^nextpow2(N);
xfourier= fft(xt,n);
p2=abs(xfourier/n);
p1=p2(1:n/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=200*10^3*(0:(n/2))/n;
figure
plot(f,p1)

nfft=n;
y=fft(u1,nfft);
y = y(1:nfft/2); 
y_mag=abs(y);
f = (0:nfft/2-1)*fs/nfft; 
figure
plot(f,y_mag)

function plot_transfer_function( tf, frequency_markers )
%PLOT_TRANSFER_FUNCTION Plots bode of a transfer function with markers
%
%   tf                - The transfer function (created using tf)
%   frequency_markers - A matrix of frequencies in Hz
%
%   Example:
%       plot_transfer_function( tf([1000], [1 1000]), [10 1000 10000] );

figure;
x_space = logspace(1,5,5000); % 5000 points between 10^1 and 10^5
x_space = 2 * pi * x_space; % to rad / sec
[mag,~,wout] = bode(tf,x_space);
mag = squeeze(mag);
wout = squeeze(wout);
mag = 20*log10(mag);
wout = wout/2/pi;
semilogx(wout,mag,'-b');
axis([min(wout) max(wout) min(mag)-10 max(mag)+10]);
[num,den] = tfdata(tf);
syms s;
d1 = digits(5);
ltx = latex(vpa(poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s)));
digits(d1);
title(strcat('$',ltx,'$'), 'Interpreter','latex', 'FontSize', 24);
xlabel('Frequency (Hz)', 'FontSize', 18);
ylabel('Magnitude (dB)', 'FontSize', 18);
grid on;
hold on;
[dbMarks,~,frequency_markers] = bode(tf,2 * pi * frequency_markers);
dbMarks = squeeze(dbMarks);
frequency_markers = squeeze(frequency_markers);
dbMarks = 20*log10(dbMarks);
frequency_markers = frequency_markers/2/pi;
Aw = cell(size(frequency_markers, 1) + 1, 1);
Aw{1} = 'Transfer function';
for i = 1 : size(frequency_markers, 1)
    semilogx(frequency_markers(i),dbMarks(i),'o');
    Aw{i+1} = sprintf('Attenuation at %.2f Hz is %.2f dB', ...
        frequency_markers(i), dbMarks(i));
end
legend(Aw,'Location','best','FontSize',12);
set(gca,'FontSize',14);
end