%Low pass filter - Chebyshev
%Project creator :Ioannis Kessopoulos 9271

a_1=9;
a_2=2;
a_3=7;
a_4=1;

m=1;
%Freequency
f_p=1.1*(3+m);
%Converting to kHz
f_p=f_p*1000;
f_s=1.6*f_p;
a_min=21.5+(max(1,a_3)-5)*0.5;
a_max=0.55+((max(1,a_4)-5)/10);

w_p=2*pi*f_p;

w_s=2*pi*f_s;

%In order to w_p=1,we adapt the frequencies
W_p=1;

W_s=w_s/w_p;

%Calculating the class of the filter
n = acosh(((10^(a_min/10)-1)/(10^(a_max/10)-1))^(1/2))/(acosh(W_s));
%Ceiling
n=ceil(n);

%Starting the proscess in order to calculate the poles
e=sqrt(10^(a_max/10)-1);

a=(1/n)*asinh(1/e);

w_hp=cosh((1/n)*acosh(1/e));

f_hp=w_hp*f_p;

%Butterworth angles-Calculation

k = 1:n;

yk = deg2rad(90) -(((2.*k-1)*pi)/(2*n));
ydeg = rad2deg(yk(1));

%Poles-Calculation

sk = sinh(a).*cos(yk);

wk = cosh(a).*sin(yk);

poles = complex(-sk,wk);

W_0 = sqrt(sk.^2+ wk.^2);

Q =  sqrt(sk.^2+wk.^2)./(2*k);

%Unit 1 - 1/RC

C=1;
R = 1/W_0(3);
kf1 = w_p;
km1 = (10^8)/kf1;
A1 = 1/(km1*kf1);
C1 = 10^(-6);
R1 = R*km1;

%Unit 2 -Sallen-Key 3rd strategy

%Normalized
R21 = 1;
C21 = 1;
r21 = 1;
r22 = 1;
R22 = Q(2);
C22 = 1/Q(2);

%Klimakopoihsh

kf2 = w_p * W_0(2);
km2 = (10^6)/kf2;
C21 = (1/(kf2*km2))*C21;
C22 = (1/(kf2*km2))*C22;
R21 = km2*R21;
R22 = km2*R22;
r21 = km2*r21;
r22 = km2*r22;

%Unit 3 -Sallen-Key 3rd strategy

%Normalized

R31 = 1;
C31 = 1;
r31 = 1;
r32 = 1;
R32 = Q(1);
C32 = 1/Q(1);

%Klimakopoihsh

kf3 = w_p * W_0(1);
km3 = (10^6)/kf3;
C31 = (1/(kf3*km3))*C31;
C32 = (1/(kf3*km3))*C32;
R31 = km3*R31;
R32 = km3*R32;
r31 = km3*r31;
r32 = km3*r32;

%Gain Adjacement

%Asked Gain is 0dB

%Gain in units 2 & 3  due to Sallen Key Strategy are =2

%Gain2=2;
%Gain3=2;

%totalGain = Gain2 * Gain3;

%Asked Gain

%aGain=1;

%a=aGain/totalGain;

%Resistors for gain adjacement unit

%RA=1000;
%RB= a*1000;

%Transfer functions before gain regulation

Tlp1 = tf(W_0(3)*w_p,[1 W_0(3)*w_p]);
plot_transfer_function(Tlp1,[w_p w_s]);

Tlp2 = tf(W_0(2)*w_p,[1 W_0(2)*w_p]);
plot_transfer_function(Tlp2,[w_p w_s]);

Tlp3 = tf(W_0(1)*w_p,[1 W_0(1)*w_p]);
plot_transfer_function(Tlp3,[w_p w_s]);

%TOTAL Transfer function

%First we calculate total transfer function between unit 1 &2

Total12 = series(Tlp1 , Tlp2);

TotalLP = series( Total12, Tlp3 );
plot_transfer_function(TotalLP,[w_p w_s]);

%Total transfer function before Gain Adjacement

invTLP = inv(TotalLP);
ltiview ({'bodemag'}, invTLP)
plot_transfer_function(invTLP,[w_p w_s]);

%Total transfer function after Gain Adjacement

%TotalLP = a * TotalLP;
%inverseTLP = inv(TotalLP);
%plot_transfer_function(inverseTLP,[w_p w_s]);

%Bodeplots

ltiview({'bodemag'}, Tlp1)

ltiview({'bodemag'}, Tlp2)

ltiview ({'bodemag'}, Tlp3)

ltiview ({'bodemag'}, TotalLP)

ltiview({'bodemag'}, Tlp1, Tlp2, Tlp3, TotalLP)

%Inverse

ltiview ({'bodemag'}, invTLP)

%ltiview ({'bodemag'}, inverseTLP)

%fourier Analysis - Spectrums

%Foureir Analysis

%Input and Output signals

t=0:10^(-5):0.001;
Uin = square(t*2*pi*2000,40);
Uout = lsim(TotalLP,Uin,t);
figure(1)
plot(Uout);
hold on
plot(Uin);

%Spectrums

Fs = 10^(5);
T=1/Fs;
L=100;
Y = fft(Uin);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(2)
plot(f,P1)
Y = fft(Uout);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(3)
plot(f,P1)


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