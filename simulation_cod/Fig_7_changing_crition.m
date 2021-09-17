clear all
d=10;

M=5; %number of users
hall = ones(1,M);%[1 2 3 4 5 ]; %all the channels
dvec = [1:1:10];

%Dall = 2*ones(1,M-1);
Dall = [5 ]; 
N=10 ;%for case I %N=10 for case II
cond=zeros(M,1);

for i = 1 : length(dvec) 
    Dall(1) = 10;
    for m =2: M
        Dall(m) = Dall(m-1) + dvec(i);
    end
    %OMA
    Poma=0; % max power needed by OMA
    for m =1 : M
        if m >1
            Poma(m) = (exp(N/(Dall(m) - Dall(m-1)))-1)/hall(m);
        else
            Poma(m) = (exp(N/(Dall(m) ))-1)/hall(m);
        end
    end

    %hybrid NOMA
    tm(1) = Dall(1);% User 1's deadline, t1
    P = zeros(M,M); 
    P(1,1) = (exp(N/tm(1))-1)/hall(1); %User 1's transmit power
    am = zeros(M,M);
    am(1) = 1/(1+P(1,1)*hall(1)); % a21 is decided by User 1's parameters
    for m = 2: M
        tm(m) = Dall(m) - Dall(m-1);
            %there are two ways. Method 1: the following is based on Lemma 1
%             tempx1 = exp((N-sum(tm(1:m-1).*log(am(m-1,1:m-1))) )/Dall(m));
%             if min(am(m-1,1:m-1))> exp(-(N-sum(tm(1:m-1).*log(am(m-1,1:m-1))) )/Dall(m))%tm(m)<=-N/log(am1(m-1))  %hybrid NOMA
%                 Pmm = (tempx1-1)/hall(m);
%                 P(m,1:m-1) =  (am(m-1,1:m-1)*tempx1-1)./am(m-1,1:m-1)/hall(m);
%                 P(m,m) = Pmm;
%             else %OMA
%                 P(m,m) =  (exp(N/tm(m))-1)/hall(m);
%             end
            
%             %%% method 2: using fmincon as the alternative 
            nonlcon = @mycons;%(x,N,tm,m,P,hall);
            options = optimoptions('fmincon','Display', 'off','MaxFunctionEvaluations', 300000); %display off
            x0 = zeros(m,1);
            A = []; % No other constraints
            b = [];
            Aeq = [];
            beq = [];
            lb = [];
            ub = [];
            x=[];
            x = fmincon(@(x) sum(x'.*tm(1:m)),x0,A,b,Aeq,beq,lb,ub,@(x) mycons(x,N,tm,m,P,hall),options);
            P(m,1:m)=x';
%             %%%%
            for n = 1 : m
                am(m,n) = 1/(1 + sum(P(n:m,n).*hall(n:m).'));
            end
    end 

    ojbe(:,i) = tm(1)*P(:,1);
    oma(:,i) = tm'.*Poma';
    for m = 2: M
        for n = 2 : m
            ojbe(m,i) = ojbe(m,i) + tm(n)*P(m,n);            
        end
    end
    
    cond(:,i) = -N./(log(am(:,1)));
end

plot(dvec,oma(2,:),dvec, ojbe(1:M,:)) 

function [c,ceq] = mycons(x,N,tm,m,P,hall)
hm = hall(m);
c(1) = N;
for i = 1: m        
    c(1) = c(1) - tm(i)*log(1 + hm*x(i)/(1+sum(hall(i:m-1).*P(i:m-1,i)'))) ;
    c(i+1) = -x(i);
end
    ceq = [];
 
end