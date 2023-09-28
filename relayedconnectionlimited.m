Ns = 2 ;
Nd = 2 ;
L = 3 ;
Nr = 4 ;
s_prik = [0.8 0.7 0.6] ;
T = 1 ;
I = 2 ;
Pp = 5 ;
alpha = 0.1 : 0.02 : 0.9 ;
eta = 0.8 ;
s_sjx = 1 ;
s_sjrik = [1.2 1.1 0.9] ;
s_rilx = [0.4 0.9 0.7] ;
s_psj = 1 ;
s_rildn = [1.1 1.0 0.9] ;
No = 1 ;
s_pdk = 1 ;
R = 1 ;
l = 2 .^ (R ./ ((1 - alpha) .* T)) - 1 ;
val = zeros( 1, length(alpha)) ;
Po = zeros( 1, length(alpha)) ;
for i = 1 : length(alpha)
    fun = 1 ;
    for j = 1 : length(L)
        fun1 = func(Ns, Nr, I, Pp, alpha(i), eta, s_sjx, s_sjrik(j), s_psj, No, s_prik(j), l(i)) ;
        fun2 = func(Nr, Nd, I, Pp, alpha(i), eta, s_rilx(j), s_rildn(j), s_prik(j), No, s_pdk, l(i)) ; 
        fun = fun .* (fun1 + fun2 - fun1 .* fun2);
    end
    Po(i) = fun ;
    Thr(i) = (1 - fun) .* R .* (1 - alpha(i)) ./ 2 ;
end
plot(alpha, Po) ;
hold on ;
%Simulation
k = No + Pp .* s_pdk ^ 2 ;
Samples = 10 ^ 6 ;
g1 = zeros( L, Samples, length(alpha) ) ;
h_psj = s_psj ./ sqrt(2) .* (randn(Samples, Ns) + 1j*randn(Samples, Ns)) ;
modh_psj = h_psj .* conj(h_psj) ;
h_sjx = s_sjx ./ sqrt(2) .* (randn(Samples, Ns) + 1j*randn(Samples, Ns)) ;
modh_sjx = h_sjx .* conj(h_sjx) ;
Psmax = zeros(Samples, 1) ;
Psj = zeros(Samples, 1) ;
Primax = zeros(Samples, 1) ;
Prik = zeros(Samples, 1) ;
gamma_1 = zeros(Samples, Ns) ;
gamma_2 = zeros(Samples, Nr) ;
Ui = zeros(Samples, 1) ;
Vi = zeros(Samples, 1) ;
Pout = zeros(1, length(alpha)) ;
for m = 1 : length(alpha)
for i = 1 : L
    for j = 1 : Ns      
        h_sjrik = s_sjrik(i) ./ sqrt(2) .* (randn(Samples, Nr) + 1j*randn(Samples, Nr)) ;
        modh_sjrik = h_sjrik .* conj(h_sjrik) ;
        Psmax = 2 .* eta .* alpha(m) .* Pp ./ ( 1 - alpha(m) ) .* sum(modh_psj, 2) ;
        Psj = min(Psmax, I ./ (modh_sjx(:, j))) ;
        modh_sjri = sum(modh_sjrik, 2) ;
        gamma_1(:, j) = Psj .* modh_sjri ./ k ;
    end
    Ui = max(gamma_1, [], 2) ;
    for j = 1 : Nr
        k = No + Pp * s_prik(i) .^ 2 ;
        h_prik = s_prik(i) / sqrt(2) .* (randn(Samples, Nr) + 1j*randn(Samples, Nr)) ;
        modh_prik = h_prik .* conj(h_prik) ;
        h_rilx = s_rilx(i) / sqrt(2) .* (randn(Samples, Nr) + 1j*randn(Samples, Nr)) ;
        modh_rilx = h_rilx .* conj(h_rilx) ;
        h_rildn = s_rildn(i) ./ sqrt(2) .* (randn(Samples, Nd) + 1j*randn(Samples, Nd)) ;
        modh_rildn = h_rildn .* conj(h_rildn) ;
        Primax = 2 .* eta .* alpha(m) .* Pp ./ ( 1 - alpha(m) ) .* sum(modh_prik, 2) ;
        Prik = min(Primax, I ./ (modh_rilx(:, j))) ;
        modh_ril = sum(modh_rildn, 2) ;
        gamma_2(:, j) = Prik .* modh_ril ./ k ;
    end
    Vi = max(gamma_2, [], 2) ;
    Pvalue = min(Ui, Vi) ;
    Pout(m) = sum(Pvalue < l(m)) / Samples ;
    Pout(m) = Pout(m) .^ Ns ;
end
end
Thr2 = (1 - Pout) .* R .* (1 - alpha) ./ 2 ;
plot(alpha, Pout, '--') ;
ylim([0 1]) ;
hold off ;