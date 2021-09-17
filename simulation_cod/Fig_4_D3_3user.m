clear all
close all

M = 3;
hall = ones(1,M);%[1 2 3 4 5 ]; %all the channels
Dall = [5 0 0 ]; 
N = 10;
Dvec = [1 :0.1: 5];
 for i = 1 : length(Dvec)
     for j = 1 : length(Dvec)
         Dall(2) = Dall(1)+Dvec(i);
         Dall(3) = Dall(2)+Dvec(j);


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

        ojbe(i,j) = tm(1)*P(1,1);
        oma(i,j) = tm(1)*Poma(1,1);
        for m = 2: M
            oma(i,j) = oma(i,j) + tm(m)*Poma(m);
            for n = 1 : m
                ojbe(i,j) = ojbe(i,j) + tm(n)*P(m,n);            
            end
        end
        
        %just for plotting, we need to force some values to inf
        if oma(i,j)>=700
            oma(i,j) = 700;
        end
        if ojbe(i,j)>=700
            ojbe(i,j) = 700;
        end
        
     end
 end

surf(Dvec,Dvec, oma )
hold

surf(Dvec,Dvec, ojbe )

function [c,ceq] = mycons(x,N,tm,m,P,hall)
hm = hall(m);
c(1) = N;
for i = 1: m        
    c(1) = c(1) - tm(i)*log(1 + hm*x(i)/(1+sum(hall(i:m-1).*P(i:m-1,i)'))) ;
    c(i+1) = -x(i);
end
    ceq = [];
 
end

 