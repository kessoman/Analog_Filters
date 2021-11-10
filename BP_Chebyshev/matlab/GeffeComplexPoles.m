%Band pass filter - Chebyshev
%Project creator :Ioannis Kessopoulos 9271

function GeffeComplexPoles( s, W, w_0, bw, pair)

q_c = w_0/bw;
C = s^2 + W^2;
D =(2*s)/ q_c;
E = 4 + (C/(q_c^2));
G = sqrt((E^2) - 4*(D^2));
Q = (1/D)*sqrt((1/2)*(E+G));
k = (s*Q)/q_c;
W = k + sqrt((k^2) -1);
w_02 = W * w_0;
w_01 = (1/W)* w_0;

if pair == 1
   Geffe_first_pair = struct ('q_c', q_c, 'C', C, 'D', D, 'E', E,'G', G,'Q', Q,'k', k,'W', W, 'w_02', w_02,'w_01', w_01)
   assignin('base', 'Geffe_first_pair', Geffe_first_pair);
else 
   Geffe_second_pair = struct ('q_c', q_c, 'C', C, 'D', D, 'E', E,'G', G,'Q', Q,'k', k,'W', W, 'w_02', w_02,'w_01', w_01)
   assignin('base', 'Geffe_second_pair', Geffe_second_pair);
end 

end

