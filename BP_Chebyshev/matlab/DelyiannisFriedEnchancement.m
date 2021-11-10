%Band pass filter - Chebyshev
%Project creator :Ioannis Kessopoulos 9271

function DelyiannisFriedEnchancement( w_o, Q, beta, num)

C_1_norm = 1;
C_2_norm = 1;
W_o = 1; 
R_1_norm = 1/(sqrt(beta));  
R_2_norm = sqrt(beta); 
k = (Q*(beta + 2) - sqrt(beta))/(2*Q - sqrt(beta));
H = (k*beta)/(2*(k-1) - beta);
k_f = w_o;
k_m = C_1_norm*(10^8 / k_f);
R_1 = R_1_norm*k_m;
R_2 = R_2_norm*k_m;
C_1 = C_1_norm*(1/(k_m*k_f));
C_2 = C_2_norm*(1/(k_m*k_f));
R_A = 4700;
R_B = (k-1)*R_A;
if num ==1
   Monada3_normalized_results = struct ('R_1_norm', R_1_norm, 'R_2_norm', R_2_norm, 'C_1_norm', C_1_norm, 'C_2_norm', C_2_norm)
   assignin('base', 'Monada3_normalized', Monada3_normalized_results);
   Monada3_results = struct ('k_f', k_f , 'k_m', k_m ,'R_1', R_1, 'R_2', R_2, 'C_1', C_1, 'C_2', C_2, 'R_A', R_A, 'R_B', R_B, 'k', k, 'H', H)
   assignin('base', 'Monada3', Monada3_results);
else
   Monada4_normalized_results = struct ('R_1_norm', R_1_norm, 'R_2_norm', R_2_norm, 'C_1_norm', C_1_norm, 'C_2_norm', C_2_norm)
   assignin('base', 'Monada4_normalized', Monada4_normalized_results);
   Monada4_results = struct ('k_f', k_f , 'k_m', k_m ,'R_1', R_1, 'R_2', R_2, 'C_1', C_1, 'C_2', C_2, 'R_A', R_A, 'R_B', R_B, 'k', k, 'H', H)
   assignin('base', 'Monada4', Monada4_results);

end

