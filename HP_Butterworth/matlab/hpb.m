%High pass filter - Butterworth
%Project creator :Ioannis Kessopoulos 9271

AEM = [9 2 7 1];

m=0;

% Pass & Stop frequency(Hz):
f_p = (4 + m)*1000;
f_s = f_p/2.6;

%Prodiagrafes bash AEM

a_min = 24 + AEM(3)*6/9;
a_max = 0.5 + AEM(4)/36;

%Gwniakes syxnothtes twn Stop & Pass Syxnotitwn

w_p = 2*pi*f_p;
w_s = 2*pi*f_s;

% Normalized circular frequencies based on w_p:
W_p = w_p/w_p;
W_s = w_p/w_s;


% Find the filter's order:
n = log10(((10^(a_min/10)-1)/(10^(a_max/10)-1)))/(2*log10(W_s/W_p));

n = ceil(n);

%Starting the proscess in order to calculate the poles

e = sqrt(10^(a_max/10)-1);

% sixnotita hmiseias isxios

%Circular

W_0 = W_p/((10^(a_max/10)-1)^(1/(2*n)));
w_0 = w_p/W_0;

%Normal

f_0 = w_0/(2*pi);

%Same n as in Low Pass Filter so Butterworth angles are the same (0, +36,
%+72, -36, -72

%Poles- Calculation

%Real part

s_1 = w_0*cos(deg2rad(0));
s_2 = w_0*cos(deg2rad(36));
s_3 = w_0*cos(deg2rad(72));

%Imaginary part

W_1 = w_0*sin(deg2rad(0));
W_2 = w_0*sin(deg2rad(36));
W_3 = w_0*sin(deg2rad(72));

%Poles- Q Part

Q_1 = 1/(2*cos(deg2rad(0)));
Q_2 = 1/(2*cos(deg2rad(36)));
Q_3 = 1/(2*cos(deg2rad(72)));

%Units- Normalized &Klimakopoihmena

%Unit 1

Q1 = Q_1;

%Normalized

C=10^7
R_1_norm = 1;
C_1_norm = 1;

%Unit's Gain

H1 = 1;

%Klimakopoihsh

k_f1 = w_0;
k_m1 = C_1_norm*(C /(k_f1));
R_1 = R_1_norm*k_m1;
C_1 = C_1_norm*(1/(k_m1*k_f1));

Unit1_normalized = struct ('k_f1' , k_f1 , 'k_m1' , k_m1 , 'R_1_norm', R_1_norm, 'C_1_norm', C_1_norm)
assignin('base', 'Unit1_normalized', Unit1_normalized);
Unit1 = struct ('k_f1', k_f1 ,'k_m1', k_m1 ,'R_1', R_1, 'C_1', C_1, 'H1', H1)
assignin('base', 'UNIT1', Unit1);

%Unit 2

Q2 = Q_2;

%Normalized

C_21_norm = 1;
C_22_norm = 1; 
R_21_norm = 1;
R_22_norm = 1;
r_21_norm = 1; 
r_22_norm = 2 - (1/Q2);

%Unit's Gain

H2 = 3 - (1/Q2);

%Klimakopoihsh

k_f2 = w_0; 
k_m2 = C_21_norm*(C / k_f2);
r_21 = r_21_norm*k_m2;
r_22 = r_22_norm*k_m2;
R_21 = R_21_norm*k_m2;
R_22 = R_22_norm*k_m2;
C_21 = C_21_norm*(1/(k_m2*k_f2));
C_22 = C_22_norm*(1/(k_m2*k_f2));

 Unit2_normalized = struct ('R_21_norm', R_21_norm, 'R_22_norm', R_22_norm, 'C_21_norm', C_21_norm, 'C_22_norm', C_22_norm, 'r_21_norm',r_21_norm, 'r_22_norm',r_22_norm)
   assignin('base', 'Unit2_normalized', Unit2_normalized);
 Unit2 = struct ('k_f2', k_f2 , 'k_m2', k_m2 ,'R_21', R_21, 'R_22', R_22, 'C_21', C_21, 'C_22', C_22,'r_21',r_21, 'r_22',r_22, 'H2', H2)
   assignin('base', 'Unit2', Unit2);

%Unit 3

Q3 = Q_3;

%Normalized

C_31_norm = 1;
C_32_norm = 1; 
R_31_norm = 1;
R_32_norm = 1;
r_31_norm = 1; 
r_32_norm = 2 - (1/Q3);

%Unit's Gain

H3 = 3 - (1/Q3);

%Klimakopoihsh

k_f3 = w_0; 
k_m3= C_1_norm*(C / k_f3);
r_31 = r_31_norm*k_m3;
r_32 = r_32_norm*k_m3;
R_31 = R_31_norm*k_m3;
R_32 = R_32_norm*k_m3;
C_31 = C_31_norm*(1/(k_m3*k_f3));
C_32 = C_32_norm*(1/(k_m3*k_f3));

 Unit3_normalized = struct ('R_31_norm', R_31_norm, 'R_32_norm', R_32_norm, 'C_31_norm', C_31_norm, 'C_32_norm', C_32_norm, 'r_31_norm',r_31_norm, 'r_32_norm',r_32_norm)
   assignin('base', 'Unit3_normalized', Unit3_normalized);
 Unit3 = struct ('k_f3', k_f3 , 'k_m3', k_m3 ,'R_31', R_31, 'R_32', R_32, 'C_31', C_31, 'C_32', C_32,'r_31',r_31, 'r_32',r_32, 'H3', H3)
   assignin('base', 'Unit33', Unit3);
   
%Transfer functions before gain adjacemnt

%Calculation

Tfhp1 = tf( [Unit1.H1 0] , [1 1/(Unit1.R_1*Unit1.C_1)] )
Tfhp2 = tf( [Unit2.H2 0 0] , [1 (w_0/Q2) w_0^2] )
Tfhp3 = tf( [Unit3.H3 0 0] , [1 (w_0/Q3) w_0^2] )

%Bodeplots

ltiview({'bodemag'}, Tfhp1)

ltiview({'bodemag'}, Tfhp2)

ltiview ({'bodemag'}, Tfhp3)

%Plots

plot_transfer_function( Tfhp1, [ f_s f_p ]);
plot_transfer_function( Tfhp2, [ f_s f_p ]);
plot_transfer_function( Tfhp3, [ f_s f_p ]);

%Gain Adjacement

%In dB

Gain_asked=10;

%In number

Gain_num=10^(0.5);

k_total=Unit1.H1*Unit2.H2*Unit3.H3;

a=Gain_num/k_total;

R_a=1000;
R_b=a*R_a;

%Transfer functions after gain adjacemnt

Total_Tfhp_first = series(a,Tfhp1);
Total_Tfhp_second = series(Tfhp2,Tfhp3);
Total_Tf = series(Total_Tfhp_first, Total_Tfhp_second)

%Bodeplots

ltiview({'bodemag'}, Total_Tf)

ltiview({'bodemag'}, Tfhp1, Tfhp2, Tfhp3, Total_Tf)

%Plots

plot_transfer_function( Total_Tf, [ f_s f_p ]);

%Inversed Transfer Function

Inv_Tfhp=inv(Total_Tf);

%Bodeplot

ltiview({'bodemag'}, Inv_Tfhp)

%Plot

plot_transfer_function( Inv_Tfhp, [ f_s f_p ]);

%Foyrier Ana;ysis

f1= (0.2*w_s) / (2*pi);
f2= (0.7*w_s) / (2*pi);
f3= (1.6*w_p) / (2*pi);
f4= (2.4*w_p) / (2*pi);
f5= (3.5*w_p) / (2*pi) ;

T = (1/100);
Fs = 10000;
dt = 1/Fs;
t = 0:dt:T_-dt;

uin= cos(2*pi*f1*t)+0.6*cos(2*pi*f2*t)+1.5*cos(2*pi*f3*t)+0.7*cos(2*pi*f4*t)+0.4*cos(2*pi*f5*t);
figure
plot(uin)

N=T/dt;
figure
lsim(Total_Tf,uin,t)
xt=lsim(Total_Tf,uin,t);
figure
plot(t,xt)

n=2^nextpow2(N);
xfourier= fft(xt,n);
p2=abs(xfourier/n);
p1=p2(1:n/2+1);
p1(2:end-1)=2*p1(2:end-1);
f=100000*(0:(n/2))/n;
figure
plot(f,p1)

nfft=n;
y=fft(uin,nfft);
y = y(1:nfft/2); 
y_mag=abs(y);
f = (0:nfft/2-1)*Fs/nfft; 
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