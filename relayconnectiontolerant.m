Ns = 4 ;
Nd = 1 ;
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
syms l ;
val = zeros( 1, length(alpha)) ;
for i = 1 : length(alpha)
    fun = @(l) 1 ;
    for j = 1 : length(L)
        fun1 = @(l) func(Ns, Nr, I, Pp, alpha(i), eta, s_sjx, s_sjrik(j), s_psj, No, s_prik(j), l) ;
        fun2 = @(l) func(Nr, Nd, I, Pp, alpha(i), eta, s_rilx(j), s_rildn(j), s_prik(j), No, s_pdk, l) ; 
        fun = @(l) fun(l) .* (fun1(l) + fun2(l) - fun1(l) .* fun2(l));
    end
    funval = @(l) (1 - fun(l)) ./ ( 1 + l );
    val(i) = 1 / ( 2 * log(2) ) * quadgk(funval, 0, 100)  ;
end
plot(alpha, val) ;