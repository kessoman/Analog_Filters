%Band Ellimination filter - Inverse Chebyshev
%Project creator :Ioannis Kessopoulos 9271

%Specifications-According to AEM
a=[9 2 7 1];
f_0 = 1800;
f_1 = 1200 + 25*a(4);
f_2 = (f_0^2)/ f_1;
D = (1/1.8)*((f_0^2 - f_1^2)/f_1);
f_3 = ((-D + sqrt((D)^2 + 4*(f_0^2)))/2);
f_4 = (f_0^2)/f_3;

amin= 30 - a(3);
amax= 0.5 + (a(4)/18);

w_0 = 2*pi*f_0; 
w_1 = 2*pi*f_1;
w_2 = 2*pi*f_2;
w_3 = 2*pi*f_3;
w_4 = 2*pi*f_4;

%Frequencies calculations -According to strategy

W_p =1;
W_s = (w_2 - w_1)/(w_4 - w_3);

bw = w_2 - w_1;
q_c = w_0 / bw;

%Taksi Filtroy

n = acosh(sqrt(((10^(amin/10) -1))/ (10^(amax/10)- 1)))/acosh(W_s);

%Ceiling n 

n = ceil(n);

%Starting the proscess in order to calculate the poles

e = (10^(amin/10)-1)^(-1/2);
a = asinh(1/e)/n;

whp = 1 / (cosh((1/n)*acosh(1/e)));

%Butterworth angles

y_k1 = 22.5;
y_k2 = -22.5;
y_k3 = 67.5;
y_k4 = -67.5;

%Poles - First Calculation

p_1 = -sinh(a)*cosd(y_k1)+1i*cosh(a)*sind(y_k1);
p_2 = -sinh(a)*cosd(y_k2)+1i*cosh(a)*sind(y_k2);
p_3 = -sinh(a)*cosd(y_k3)+1i*cosh(a)*sind(y_k3);
p_4 = -sinh(a)*cosd(y_k4)+1i*cosh(a)*sind(y_k4);

%Calculating factors based on pole calculation

W_12 = abs(p_1); 
Q_12 = W_12/(2*abs(real(p_1)));

W_34 = abs(p_3);
Q_34 = W_34/(2*abs(real(p_3)));

%Inversing them due to Inverse Chebyshev policy

IW_12 = 1/W_12;
IW_34 = 1/W_34;

%Klimakopoihsh twn polwn gia metafora se pedio syxnothtwn

I_W12 = IW_12*(1/W_p);
I_W34 = IW_34*(1/W_p);

%Finding the zeros of the Inverse Checy calculation

W_z12 = sec(0/(2*n));
W_z34 = sec(2*pi/(2*n));

%Klimakopoihsh mhdenikwn

W_z12 = W_z12 * (1/W_p);
W_z34 = W_z34 * (1/W_p);

%Reversin the InverseChebyshev poles and zeroes

IW_012 = 1/I_W12;
IW_034 = 1/I_W34;

IW_z12 = 1/W_z12;
IW_z34 = 1/W_z34;

%Final Real and imaginary parts of the Inverse Chebyshevs poles

S_12HP = -IW_012/(2*Q_12);
W_12HP = sqrt(IW_012^2 - S_12HP^2);

S_34HP = -IW_034/(2*Q_34);
W_34HP = sqrt(IW_034^2 - S_34HP^2);

%Using Geffe;s Algorithm to transform the poles

%First pair of poles

S_r1 = abs(S_12HP);
W_r1 = abs(W_12HP);

C_21 = S_r1^2 + W_r1^2;
D_1 = 2*S_r1/q_c;
E_1 = 4 + C_21/((q_c)^2);

G_1 = sqrt(E_1^2-4*(D_1^2));
Q_1 = sqrt((E_1+G_1)/2)/D_1;
k_21 = S_r1*Q_1/q_c;
W_1 = k_21 +sqrt((k_21^2) - 1);
w_012 = W_1*w_0;
w_011 = w_0/W_1;

%Second pait

S_r2 = abs(S_34HP);
W_r2 = abs(W_34HP);

C_2 = S_r2^2 + W_r2^2;
D_2 = 2*S_r2/q_c;
E_2 = 4 + C_2/((q_c)^2);

G_2 = sqrt(E_2^2-4*(D_2^2));
Q_2 = sqrt((E_2+G_2)/2)/D_2;
k_22 = S_r2*Q_2/q_c;
W_2 = k_22 +sqrt((k_22^2) - 1);
w_022 = W_2*w_0;
w_021 = w_0/W_2;

%From the above transformations occured : 2 pairs of zeroes so we need to
%transform them

K_3 = 2 + IW_z12^2/q_c^2;
x_3 = (K_3 +sqrt(K_3^2 -4))/2;
w_z11 = w_0 * sqrt(x_3);
w_z12 = w_0/sqrt(x_3);

%And for the second pair

K_4 = 2 + IW_z34^2/q_c^2;
x_4 = (K_4 +sqrt(K_4^2 -4))/2;
w_z21 = w_0 * sqrt(x_4);
w_z22 = w_0/sqrt(x_4);

%Next we design all 4 units follwing different strategies according to our
%AEM

%Unit 1 - Normalized
%w_z> w_0 for this unit so we follow the LPN strategy

W_o1_normalized = 1;
W_z1_normalized = w_z11/w_011;
R_11_normalized = 1;
C_normalized = 1/(2*Q_1);
R_12_normalizes = 4*(Q_1^2);
R_15_normalized = (4*Q_1^2)/(W_z1_normalized^2 -1);
H_high_1 = 1 / (1 + ((W_z1_normalized^2)/(2*Q_1^2)));
R_14_normalized = 1;
R_13_normalized = (W_z1_normalized^2)/(2*Q_1^2);

%Unit 1 -Klimakopoihsh

k_f1 = w_011; 
k_m1 = C_normalized*(10^8 / k_f1);  
R_11 = R_11_normalized*k_m1;
C11 = 0.01*(10^-6);
R_12 = R_12_normalizes*k_m1;
R_15 = R_15_normalized*k_m1;
H_low_1 = H_high_1*((w_z11/w_011)^2);
R_14 = R_14_normalized*k_m1;
R_13 = R_13_normalized*k_m1;
C12=C_normalized/(k_m1 * k_f1);

Monada1_normalized = struct ('W_z_norm', W_z1_normalized , 'W_o_norm', W_o1_normalized , 'R_1_norm', R_11_normalized, 'R_2_norm', R_12_normalizes, 'R_3_norm', R_13_normalized, 'R_4_norm', R_14_normalized, 'R_5_norm', R_15_normalized, 'C_norm', C_normalized)
assignin('base', 'Monada2_normalized', Monada1_normalized);
Monada1_results = struct ('k_f', k_f1 , 'k_m', k_m1 ,'R_1', R_11, 'R_2', R_12, 'R_3', R_13, 'R_4', R_14, 'R_5', R_15, 'C11', C11, 'C12' , C12 , 'H_high', H_high_1, 'H_low', H_low_1)
assignin('base', 'Monada2', Monada1_results);

%Unit 2 - Normalized
%w_z< w_0 for this unit so we follow the HPN strategy

W_02_normalized = 1;
W_z2_normalized = w_z12/w_012;
k_21= (W_02_normalized^2/W_z2_normalized^2) -1;
k_22 = ((2+k_21)*Q_1^2)/((2+k_21)*Q_1^2 +1);
C_2_normalized = 1/(Q_1*(2+k_21));
C_21_normalized = k_21*C_2_normalized;
R_21_normalized = 1;
R_22_normalized = (Q_1^2)*((k_21 + 2)^2);
R_23_normalized = 1;
R_24_normalized = (Q_1^2)*(k_21 + 2);
H_high_2 = k_22*(W_02_normalized^2/W_z2_normalized^2);
H_low_2 = H_high_2*((w_z12/w_012)^2);

%Unit 2 - Klimakopoihsh

k_f2 = w_012;
k_m2 = C_2_normalized*(10^8 / k_f2);   
C_21 =  0.01*(10^-6);
C_22 = C_2_normalized*(1/(k_m2*k_f2));
R_21 = R_21_normalized*k_m2;
R_22 = R_22_normalized*k_m2;
R_23 = R_23_normalized*k_m2;
R_24 = R_24_normalized*k_m2;

Monada2_normalized = struct ('W_z_norm', W_z2_normalized , 'W_o_norm', W_02_normalized , 'k_1',k_21,'k_2',k_22,'R_1_norm', R_21_normalized, 'R_2_norm', R_22_normalized, 'R_3_norm', R_23_normalized, 'R_4_norm', R_24_normalized, 'C_1_norm', C_21_normalized, 'C_norm', C_2_normalized)
assignin('base', 'Monada1_normalized', Monada1_normalized);
Monada2_results = struct ('k_f', k_f2 , 'k_m', k_m2 ,'R_1', R_21, 'R_2', R_22, 'R_3', R_23, 'R_4', R_24, 'C_21', C_21, 'C_22', C_22, 'H_high', H_high_2, 'H_low', H_low_2)
assignin('base', 'Monada1', Monada1_results);

%Unit 3 - Normalized
%w_z> w_0 for this unit so we follow the LPN strategy

W_03_normalized = 1;
W_z3_normalized = w_z21/w_021;
R_31_normalized = 1;
C_3_normalized = 1/(2*Q_2);
R_32_normalized = 4*(Q_2^2);
R_35_normalized = (4*Q_2^2)/(W_z3_normalized^2 -1);
H_high_3 = 1 / (1 + ((W_z3_normalized^2)/(2*Q_2^2)));
R_34_normalized = 1;
R_33_normalized = (W_z3_normalized^2)/(2*Q_2^2);

%Unit 3 - Klimakopoihsh

k_f3 = w_021; 
k_m3 = C_3_normalized*(10^8 / k_f3);  
R_31 = R_31_normalized*k_m3;
C31 = 0.01*(10^-6);
R_32 = R_32_normalized*k_m3;
R_35 = R_35_normalized*k_m3;
H_low_3 = H_high_3*((w_z21/w_021)^2);
R_34 = R_34_normalized*k_m3;
R_33 = R_33_normalized*k_m3;
C32=C_3_normalized/(k_m1 * k_f1);

Monada3_normalized = struct ('W_z_norm', W_z3_normalized , 'W_o_norm', W_03_normalized , 'R_1_norm', R_31_normalized, 'R_2_norm', R_32_normalized, 'R_3_norm', R_33_normalized, 'R_4_norm', R_34_normalized, 'R_5_norm', R_35_normalized, 'C_norm', C_3_normalized)
assignin('base', 'Monada2_normalized', Monada3_normalized);
Monada3_results = struct ('k_f', k_f3 , 'k_m', k_m3 ,'R_1', R_31, 'R_2', R_32, 'R_3', R_33, 'R_4', R_34, 'R_5', R_35, 'C31', C31, 'C32', C32 , 'H_high', H_high_3, 'H_low', H_low_3)
assignin('base', 'Monada2', Monada3_results);

%Unit 4 - Normalized
%w_z < w_0 for this unit so we follow the HPN strategy

W_04_normalized = 1;
W_z4_normalized = w_z22/w_022;
k_4 = (W_04_normalized^2/W_z4_normalized^2) -1;
k_42 = ((2+k_4)*Q_2^2)/((2+k_4)*Q_2^2 +1);
C_4_normalized = 1/(Q_2*(2+k_4));
C_41_normalized = k_4*C_4_normalized;
R_41_normalized = 1;
R_42_normalized = (Q_2^2)*((k_4 + 2)^2);
R_43_normalized = 1;
R_44_normalized = (Q_2^2)*(k_4 + 2);
H_high_4 = k_42*(W_04_normalized^2/W_z4_normalized^2);
H_low_4 = H_high_4*((w_z22/w_022)^2);

%Unit 4 Klimakopoihsh

k_f4 = w_022;
k_m4 = C_4_normalized*(10^8 / k_f4);   
C_41 =  0.01*(10^-6);
C_42 = C_4_normalized*(1/(k_m2*k_f2));
R_41 = R_41_normalized*k_m4;
R_42 = R_42_normalized*k_m4;
R_43 = R_43_normalized*k_m4;
R_44 = R_44_normalized*k_m4;

Monada4_normalized = struct ('W_z_norm', W_z4_normalized , 'W_o_norm', W_04_normalized , 'k_1',k_4,'k_2',k_42,'R_1_norm', R_41_normalized, 'R_2_norm', R_42_normalized, 'R_3_norm', R_43_normalized, 'R_4_norm', R_44_normalized, 'C_1_norm', C_41_normalized, 'C_norm', C_4_normalized)
assignin('base', 'Monada1_normalized', Monada4_normalized);
Monada4_results = struct ('k_f', k_f4 , 'k_m', k_m4 ,'R_1', R_41, 'R_2', R_42, 'R_3', R_43, 'R_4', R_44, 'C_1', C_41, 'C', C_42, 'H_high', H_high_4, 'H_low', H_low_4)
assignin('base', 'Monada1', Monada4_results);

%Synarthseis metaforas

Tfbe_1 = tf([H_high_1 0 H_high_1*(w_z12^2)],[1 w_011/Q_1 w_011^2]);

Tfbe_2 = tf([H_high_2 0 H_high_2*(w_z11^2)],[1 w_012/Q_1 w_012^2]);

Tfbe_3= tf([H_high_3 0 H_high_3*(w_z22^2)],[1 w_021/Q_2 w_021^2]);

Tfbe_4 = tf([H_high_4 0 H_high_4*(w_z21^2)],[1 w_022/Q_2 w_022^2]);

%Gain Adjacement

Tfbe_total = Tfbe_1*Tfbe_2*Tfbe_3*Tfbe_4;

H_high_total = H_high_1*H_high_2*H_high_3*H_high_4;

%20log(a*H_high_toal)=10

a = (10^0.5)/H_high_total;

%Resistors for the Gain Adjacement unit

R_A = 1 * 10^3;
R_B = R_A*a;

TBE = a*Tfbe_total;

ltiview({'bodemag'}, Tfbe_1)

ltiview({'bodemag'}, Tfbe_2)

ltiview({'bodemag'}, Tfbe_3)

ltiview({'bodemag'}, Tfbe_4)

ltiview({'bodemag'}, Tfbe_1, Tfbe_2, Tfbe_3, Tfbe_4)

ltiview({'bodemag'}, TBE)

ltiview({'bodemag'}, Tfbe_1, Tfbe_2, Tfbe_3, Tfbe_4, TBE)

plot_transfer_function(Tfbe_1, [f_1 f_2 f_3 f_4]);
 plot_transfer_function(Tfbe_2, [f_1 f_2 f_3 f_4]);
 plot_transfer_function(Tfbe_3, [f_1 f_2 f_3 f_4]);
 plot_transfer_function(Tfbe_4, [f_1 f_2 f_3 f_4]);
 plot_transfer_function(Tfbe_total, [f_1 f_2 f_3 f_4]);
 plot_transfer_function(TBE, [f_1 f_2 f_3 f_4]);


ltiview({'bodemag'}, Tfbe_total)

InvSys = inv (Tfbe_total);
ltiview({'bodemag'}, InvSys);
plot_transfer_function(InvSys, [f_1 f_2 f_3 f_4]);

InvSys_new = inv (TBE);
ltiview({'bodemag'}, InvSys_new);
plot_transfer_function(InvSys_new, [f_0 f_1 f_2 f_3 f_4]);

%Fourier Analysis
%Input%Output Signals- Spectrums

f1 = (w_0- (w_0-w_3)/2)/(2*pi);
f2 = (w_0+(w_0+w_3)/2)/(2*pi);
f3 = 0.5*w_1 /(2*pi);
f4 = (2.4*w_2) / (2*pi);
f5 = (3.5*w_2)/(2*pi) ;

fs = 200*10^3;
T = 0.002;
dt = 1/fs;
t = 0:dt:(T);

uin = 0.8* cos(2*pi*f1*t)+ cos(2*pi* f2*t)+ cos(2*pi* f3*t)+0.8*cos(2*pi*f4*t)+0.4*cos(2*pi*f5*t);
figure
plot(uin)

N=T/dt;
figure
lsim(TBE,uin,t)

xt=lsim(TBE,uin,t);
figure
plot(t,xt)

n=2^nextpow2(N);
xF= fft(xt,n);

psecond=abs(xF/n);
pfirst=psecond(1:n/2+1);
pfirst(2:end-1)=2*pfirst(2:end-1);

f=200*10^3*(0:(n/2))/n;
figure
plot(f,pfirst)

nfft=n;
y=fft(uin,nfft);
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