Ns = 4 ;
Nd = 1 ;
T = 1 ;
I = 2 ;
Pp = 5 ;
alpha = 0.1 : 0.1 : 0.9 ;
eta = 0.8 ;
s_sjx = 1 ;
s_sjdk = 1 ;
s_psj = 1 ;
No = 1 ;
s_pdk = 1 ;
syms l ;
val = zeros( 1, length(alpha)) ;
for i = 1 : length(alpha)
fun = @(l) (1 - func(Ns, Nd, I, Pp, alpha(i), eta, s_sjx, s_sjdk, s_psj, No, s_pdk, l)) ./ (1 + l) ;
val(i) = 1 / ( 2 * log(2) ) * quadgk(fun, 0, 1000) ;
end
plot(alpha, val) ;
%plot(alpha, Pr) ;
hold on ;
Ns = 4 ;
Nd = 1 ;
Samples = 10 ^ 6 ;
%l = 0.1 : 0.01 : 40 ;
a = eta .* alpha .* Pp ./ (1 - alpha) ;
k = No + Pp .* s_pdk ;
Pout = zeros( Ns, length(alpha) ) ;
for j = 1 : Ns 
h_sjx = 1/sqrt(2) * (randn( Samples, 1) + 1j*randn( Samples, 1)) ;
h_sjdk = 1/sqrt(2) * (randn( Samples, Nd) + 1j*randn( Samples, Nd)) ;%
h_psj = 1/sqrt(2) * (randn( Samples, Ns) + 1j*randn( Samples, Ns)) ;
modh_sjx = (h_sjx) .* conj(h_sjx) ;
modh_sjdk = sum(h_sjdk .* conj(h_sjdk), 2) ;
modh_psj = sum(h_psj .* conj(h_psj), 2) ;
calc1 = a.*modh_psj .* (a .* modh_psj < I ./ modh_sjx) + I ./ modh_sjx .* (a .* modh_psj > I ./ modh_sjx) ;
Pout(j, :) = sum( log2( 1 + calc1 .* modh_sjdk ./ k ) ) / Samples ;
end
if( Ns == 1 )
    Po = Pout ;
else
    Po = max(Pout) ;
end
Thr2 = Pout .* (1 - alpha) ./ 2 ;
plot(alpha, Po, '--') ;
ylim([0 10]) ;
%legend('Theory Nd = 1', 'Simulation Nd = 1','Theory Nd = 2', 'Simulation Nd = 2','Theory Nd = 3', 'Simulation Nd = 3') ;
%xlabel('Alpha') ;
%ylabel('Outage Probability') ;
%ylabel('Throughput') ;
%xlim([0.1 0.9]) ;
%title('Outage probability vs alpha for delay limited case, with Ns=4,Ith=2,Pp=5') ;
hold off ;

