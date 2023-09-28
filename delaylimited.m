Ns = 1 ;
Nd = 4 ;
T = 1 ;
I = 2 ;
Pp = 5 ;
alpha = 0.1 : 0.1 : 1 ;
eta = 0.8 ;
a = eta .* alpha .* Pp ./ (1 - alpha) ;
s_sjx = 1 ;
s_sjdk = 1 ;
s_psj = 1 ;
val1 = I ./ (a .* s_sjx) ;
R = 1 ;
l = 2 .^ (R ./ ((1 - alpha) .* T)) - 1 ;
No = 1 ;
s_pdk = 1 ;
k = No + Pp .* s_pdk ;
val2 = l .* k ./ ( s_sjdk .* a ) ;
term2 = 2 ./ factorial(Ns - 1) .* (val1 ./ s_psj ^ 2) .^ (Ns ./ 2) .* besselk(Ns, 2 .* sqrt(val1 ./ s_psj.^2));
term3 = 0 ;
for i = 0 : Nd - 1
    term3 = term3 + val2 .^ i ./ (factorial(Ns - 1) .* factorial(i) * s_psj ^ (2 * i)) .* 2 .* ( val2 ./ s_psj^2) .^ ((Ns - i)./ 2) .* besselk(Ns - i, 2 .* sqrt( val2 ./ s_psj^2 )) ;
end
term4 = 0 ;
for i = 0 : Nd - 1
    term4 = term4 + val2 .^ i ./ (factorial(Ns - 1) .* factorial(i) * s_psj ^ (2 * i)) .* 2 .* ( (val1 + val2) ./ s_psj^2 ) .^ ((Ns - i)./ 2) .* besselk(Ns - i, 2 .* sqrt((val1 + val2) / s_psj^2)) ; 
end
term5 = 0 ;
for i = 0 : Ns - 1
    term5 = term5 + val1 .^ i * s_sjx ^ (2 * (i - 1)) ./ (factorial(i) * s_psj ^ (2 * i)) .* 2 .* ( val1 * s_sjx ^ 4 / s_psj ^ 2) .^ ( (1 - i) ./ 2 ) .* besselk( 1 - i, 2 .* sqrt( val1 ./ s_psj^2 )) ;
end
term6 = 0 ;
for i = 0 : Nd - 1
    for j = 0 : Ns - 1
        term6 = term6 + val2 .^ i .* val1 .^ (j - i) * s_sjx ^ (2 * (j - i - 1)) ./ ( factorial(i) .* factorial(j) * s_psj ^ (2 * j)) .* 2 .* ( val1 * s_sjx^2 ./ s_psj^2 ./ ( val2 .* a ./ I + 1 / s_sjx^2 )) .^ ( (i - j + 1) ./ 2 ) .* besselk( i - j + 1, 2 .* sqrt( val1 * s_sjx^2 ./ s_psj^2 .* ( val2 .* a ./ I + 1 / s_sjx^2))) ;
    end
end
Pr = (1 - term2 - term3 + term4 + term5 - term6) .^ Ns ;
Thr = (1 - Pr) .* R .* (1 - alpha) ./ 2 ;
plot(alpha, Thr) ;
hold on ;
Ns = 1 ;
Nd = 4 ;
Samples = 10 ^ 6 ;
Pout = zeros( Ns, length(alpha) ) ;
for j = 1 : Ns 
h_sjx = 1/sqrt(2) * (randn( Samples, 1) + 1j*randn( Samples, 1)) ;
h_sjdk = 1/sqrt(2) * (randn( Samples, Nd) + 1j*randn( Samples, Nd)) ;
h_psj = 1/sqrt(2) * (randn( Samples, Ns) + 1j*randn( Samples, Ns)) ;
modh_sjx = (h_sjx) .* conj(h_sjx) ;
modh_sjdk = sum(h_sjdk .* conj(h_sjdk), 2) ;
modh_psj = sum(h_psj .* conj(h_psj), 2) ;
calc1 = a.*modh_psj .* (a .* modh_psj < I ./ modh_sjx) + I ./ modh_sjx .* (a .* modh_psj > I ./ modh_sjx) ;
Pout(j, :) = sum( calc1 .* modh_sjdk < l * k ) / Samples ;
Pout(j, :) = Pout(j, :) .^ Ns ;
end
if( Ns == 1 )
    Po = Pout ;
else
    Po = min(Pout) ;
end
Thr2 = (1 - Po) .* R .* (1 - alpha) ./ 2 ;
plot(alpha, Thr2, '--') ;
title('Throughput vs alpha for delay limited transmission for Ns=1,Nd=4') ;
xlabel('alpha') ;
ylabel('Throughput') ;
legend('Theory','Simulation') ;
xlim([0.1 0.9]) ;
grid on ;
hold off ;