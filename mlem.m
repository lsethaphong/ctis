function maxE_ntF = mart(ntF_hat, ntG_hat, ntG_prior, ntF_prior, Hmat, psinv_Hmat,cycle_stop)
nn=length(ntF_hat);
mm=length(ntG_hat);
ntF_k1 = zeros(nn,1);
ntG_k = zeros(mm,1);
denomiator=zeros(nn,1);
prob_quot = float(0.0);
numerator = zeros(nn,1);
for cs = 1:nn
    for hs = 1:mm
        denominator(cs) += Hmat(hs,cs);
    end
end
cyclic = 1
ntG_k = psinv_Hmat*ntF_prior;
ntF_k1 = Hmat*ntG_k;
if cycle_stop == 1
   while (cyclic ~= cycle_stop)&&(prob_quot < 1.0)
       % calculate Hmn*Gm/G_hat_m
       numerator = psinv_Hmat*ntG_prior;
       for i = 1:nn
           for j = 1:mm
                 numerator(i)=psinv_Hmat(j,i)*ntG_prior(j);
           end
           prob_quot(i) = denominator(i)/numerator(i);
           ntF_hat(i) = ntF_prior(i)*ntG_scale(i)*prob_quot(i);
       end
       ntG_hat = psinv_Hmat*ntF_hat;
       ntG_prior = ;
   end
else
    message = 'cycle_stop must be 1 or greater'
end
maxE_ntF=ntF_hat;