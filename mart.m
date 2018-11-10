function maxE_ntF = mart(ntG_prior, ntF_prior, Ho, Hp,cycle_stop);
nn=length(ntF_prior);
mm=length(ntG_prior);
ntF_k1 = zeros(nn,1); % k+1 estimate
ntF_k = zeros(nn,1); % kth estimate
ntG_k = zeros(mm,1); % kth estimate
denomiator=zeros(nn,1);
prob_quot = 0.0; % was float
numerator = zeros(nn,1);
ntG_k = Ho*ntF_prior; % initial forward estimate
ntF_k = Hp*ntG_k; % initial back projection
numerator = Hp*ntG_prior; % what you should see
cyclic = 0;
%   while (cyclic ~= cycle_stop)&&(prob_quot < 1.0)
    while (cyclic < cycle_stop)
    denominator = Hp*ntG_k;
       for i = 1:nn
           if denominator(i) > 0.0
               prob_quot=numerator(i)/denominator(i);
           else
               prob_quot= 0.0;
           end
           ntF_k1(i)= ntF_k(i)*prob_quot;
       end
       % least squares calculation
       ntF_k = ntF_k1;
       ntG_k = Ho*ntF_k1;
       cyclic= cyclic + 1;
   end
%cnt
%message = strcat('iterations: ',char(cyclic))
cyclic
maxE_ntF=ntF_k1; % returns only a voxel representation