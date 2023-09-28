clc; clear all;
close all;
N=1e5;

Ns=2; %No of antenna at S
Nd=3; %No of antenna at D

vardB = 0;
var = 10^(vardB/10);  %variance of awgn noise + interf from PTx to S and relay


R=1; %bits/sec/Hz
gamma_th = (2^(R))-1;

eta_eff = .8; %efficiency of energy conversion

beta = 4;  %path loss exponent


% S_loc = [0,0];
% D_loc = [1,0];
% X_loc = [1,1];
% T_loc = [0,1];



%channel from source to primary Rx
h_sx = zeros(N,Ns);

%channel from primary Tx to S
h_ts = zeros(N,Ns);

h_sd = zeros(N,Ns,Nd);


dis_SX = sqrt(2); %sqrt((1-0)^2 + (1-0)^2)
mu_SX = dis_SX^(-beta);


dis_TS = 1; %sqrt((0-0)^2 + (1-0)^2)
mu_TS = dis_TS^(-beta);

dis_SD = 1; %sqrt((0-0)^2 + (0-1)^2)
mu_SD = dis_SD^(-beta);


%from source to primary Rx
for j=1:Ns
    h_sx(:,j) = sqrt(mu_SX/2)*(randn(N,1)+1j*randn(N,1));
end


%from primary Tx to S
for j=1:Ns
    h_ts(:,j) = sqrt(mu_TS/2)*(randn(N,1)+1j*randn(N,1));
end

%from S to D
for j=1:Ns
        for m=1:Nd
            h_sd(:,j,m)= sqrt(mu_SD/2)*(randn(N,1)+1j*randn(N,1));
        end
end

%----------------  end channel --------------------


alpha= 0.6;   %time switch factor
b=2*eta_eff*alpha/(1-alpha);
Ith_dB=5;
Ith=10^(Ith_dB/10);

count=0;
for Pmax_dB = -20:2:15 % Primary Tx power
    count = count+1
    Pmax = 10^(Pmax_dB/10);
    gamma_P_db = Pmax_dB - vardB;
       
    %Ith = eta * Pmax;
        
    
    %---start power allocation----
       
    sum_TS = zeros(N,1);
    %for T to S
       for j=1:Ns
          sum_TS(:,1) = sum_TS(:,1) + abs(h_ts(:,j)).^2;           
       end
    Ps_max = b*Pmax*sum_TS;
    
    
        
Ps=zeros(N,Ns);

    for j=1:Ns
        Ps(:,j) = min(Ps_max , Ith./(abs(h_sx(:,j)).^2));
    end
    
     %----------end power allocation-------
  
    %MRC combined SNR
    sum_SD = zeros(N,Ns);
        
    %for S to D
    
        for j=1:Ns
            for m=1:Nd
                gamma_SD(:,j,m) = Ps(:,j) .* (abs(h_sd(:,j,m)).^2)./var;
            end
        end
    
    
    %MRC for SD link
    
        for j=1:Ns
           for m=1:Nd
                sum_SD(:,j) = sum_SD(:,j)+ gamma_SD(:,j,m);
            end
        end
   
    
    
    %---------E2E SNR calculation for best option
    max_gammaSD = zeros(N,1);
    
    %Find the best transmit antenna at S 
  
        for n=1:N
            max_gammaSD(n,1) = max(sum_SD(n,:));            
        end
       
    Pout_sim = sum(max_gammaSD < gamma_th)/N;
    Tput_sim = (1-Pout_sim)*R*(1-alpha)/2;
    %---------------end simulation-----
    
    
    
    
%-----start analysis-----         
Pbar=Pmax/var;
Ibar=Ith/var;
eta=Ibar/Pbar;


%------------check integrals for a given y ----
y=.6;

T1_int= quadgk(@(x)x.^(Nd-1).*exp(-x./mu_SD)./(mu_SD^Nd.*factorial(Nd-1)), 0, gamma_th./(b.*Pbar.*y))
 

sum1=0;
for t=0:Nd-1
sum1=sum1+(gamma_th./(b.*Pbar.*y.*mu_SD))^t*(1/factorial(t));
end
T1_sum=1-exp(-(gamma_th./(b.*Pbar.*y.*mu_SD))).*sum1
   
T1_incompgma=gammainc(gamma_th./(b.*Pbar.*y.*mu_SD), Nd, 'lower').*...
    gamma(gamma_th./(b.*Pbar.*y.*mu_SD))./factorial(Nd-1)
  


   T2_int=   quadgk(@(x)x.^(Nd-1).*exp(-x.*(1./mu_SD+Ibar./(gamma_th.*mu_SX)))./...
       (mu_SD.^Nd.*factorial(Nd-1)), gamma_th./(b.*Pbar.*y), inf)
 

sum2=0;
for t=0:Nd-1
sum2=sum2+(gamma_th./(b.*Pbar.*y))^t*(1/factorial(t))./...
    (1./mu_SD+Ibar./(gamma_th*mu_SX))^(Nd-t);
end
T2_sum= exp(-(gamma_th./(b.*Pbar.*y.*mu_SD) + eta./(b.*y.*mu_SX))).*sum2
   
T2_incompgma=  (1+Ibar*mu_SD./(gamma_th.*mu_SX))^(-Nd).*...
        gammainc((gamma_th./(b.*Pbar.*y.*mu_SD) + (eta./(b.*y.*mu_SX))), Nd, 'upper')...
        .*gamma((gamma_th./(b.*Pbar.*y.*mu_SD) + (eta./(b.*y.*mu_SX))))./gamma(Nd)
    
%------------end check integrals----

%---------Pout ananlytical expression----
I1=zeros(N,1);
for n=1:N
    y=sum_TS(n);
    
sum1=0;
for t=0:Nd-1
sum1=sum1+(gamma_th./(b.*Pbar.*y.*mu_SD))^t*(1/factorial(t));
end
T1_sum=1-exp(-(gamma_th./(b.*Pbar.*y.*mu_SD))).*sum1;

sum2=0;
for t=0:Nd-1
sum2=sum2+(gamma_th./(b.*Pbar.*y))^t*(1/factorial(t))./...
    (1./mu_SD+Ibar./(gamma_th*mu_SX))^(Nd-t);
end
T2_sum= exp(-(gamma_th./(b.*Pbar.*y.*mu_SD) + eta./(b.*y.*mu_SX))).*sum2;
   
I1(n)=(T1_sum+T2_sum)^Ns;
   
end

Pout_analysis = mean(I1);
Tput_analysis = (1-Pout_analysis)*R*(1-alpha)/2;

result(count,:) = [gamma_P_db, Pout_sim, Pout_analysis,  Tput_sim , Tput_analysis]; 

% Y1=( 1./gamma(Nd).*( gammainc(gamma_th./(b.*Pbar.*y.*mu_SD), Nd, 'lower')...
%         .*gamma(gamma_th./(b.*Pbar.*y.*mu_SD))+...
%         (1+Ibar*mu_SD./(gamma_th.*mu_SX))^(-Nd).*...
%         gammainc((gamma_th./(b.*Pbar.*y.*mu_SD) + (eta./mu_SX)), Nd, 'upper')...
%         .*gamma((gamma_th./(b.*Pbar.*y.*mu_SD) + (eta./mu_SX))))).^Ns;

%-------------
%     func=@(y) ( 1./gamma(Nd).*( gammainc(gamma_th./(b.*Pbar.*y.*mu_SD), Nd, 'lower')...
%         .*gamma(gamma_th./(b.*Pbar.*y.*mu_SD))+...
%         (1+Ibar*mu_SD./(gamma_th.*mu_SX))^(-Nd).*...
%         gammainc((gamma_th./(b.*Pbar.*y.*mu_SD) + (eta./mu_SX)), Nd, 'upper')...
%         .*gamma((gamma_th./(b.*Pbar.*y.*mu_SD) + (eta./mu_SX))))).^Ns.*...
%         y.^(Ns-1).*exp(-y./mu_TS)./(mu_TS.^Ns.*gamma(Ns));
%     
  %Pout_analysis = quadgk( func, 0.01, 1e2)   ;

%   Tput_analysis = (1-Pout_analysis)*R*(1-alpha)/2;
%     
% 

end



figure(1);
semilogy(result(:,1),result(:,2),'bo','LineWidth',1);
 hold on;
semilogy(result(:,1),result(:,3),'k-','LineWidth',1);
legend('Simulation', 'Analysis');
title('Pout vs SNR') ;



