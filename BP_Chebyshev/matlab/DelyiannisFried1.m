%Band pass filter - Chebyshev
%Project creator :Ioannis Kessopoulos 9271
function DelyiannisFried1( w_o, Q, num)

C_1_norm = 1/(2*Q);
C_2_norm = 1/(2*Q);
R_1_norm = 1;   
R_2_norm = 4*Q^2; 
H = 2*Q^2; 
k_f = w_o;
k_m = C_1_norm*(10^7 / k_f);
R_1 = R_1_norm*k_m;
R_2 = R_2_norm*k_m;
C_1 = C_1_norm*(1/(k_m*k_f));
C_2 = C_2_norm*(1/(k_m*k_f));

if num ==1
   Monada1_normalized_results = struct ('R_1_norm', R_1_norm, 'R_2_norm', R_2_norm, 'C_1_norm', C_1_norm, 'C_2_norm', C_2_norm)
   assignin('base', 'Monada1_normalized', Monada1_normalized_results);
   Monada1_results = struct ('k_f', k_f , 'k_m', k_m ,'R_1', R_1, 'R_2', R_2, 'C_1', C_1, 'C_2', C_2, 'H', H)
   assignin('base', 'Monada1', Monada1_results);
elseif num==2 
    Monada2_normalized_results = struct ('R_1_norm', R_1_norm, 'R_2_norm', R_2_norm, 'C_1_norm', C_1_norm, 'C_2_norm', C_2_norm)
   assignin('base', 'Monada2_normalized', Monada2_normalized_results);
   Monada2_results = struct ('k_f', k_f , 'k_m', k_m ,'R_1', R_1, 'R_2', R_2, 'C_1', C_1, 'C_2', C_2, 'H', H)
   assignin('base', 'Monada2', Monada2_results); 
elseif num==3  
    Monada3_normalized_results = struct ('R_1_norm', R_1_norm, 'R_2_norm', R_2_norm, 'C_1_norm', C_1_norm, 'C_2_norm', C_2_norm)
   assignin('base', 'Monada3_normalized', Monada3_normalized_results);
   Monada3_results = struct ('k_f', k_f , 'k_m', k_m ,'R_1', R_1, 'R_2', R_2, 'C_1', C_1, 'C_2', C_2, 'H', H)
   assignin('base', 'Monada3', Monada3_results);
elseif num==4  
    Monada4_normalized_results = struct ('R_1_norm', R_1_norm, 'R_2_norm', R_2_norm, 'C_1_norm', C_1_norm, 'C_2_norm', C_2_norm)
   assignin('base', 'Monada4_normalized', Monada4_normalized_results);
   Monada4_results = struct ('k_f', k_f , 'k_m', k_m ,'R_1', R_1, 'R_2', R_2, 'C_1', C_1, 'C_2', C_2, 'H', H)
   assignin('base', 'Monada4', Monada4_results);
end



end

