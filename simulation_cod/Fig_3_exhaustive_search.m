clear all
close all
M=3; %number of users
hall = ones(1,M); %all the channels
Dall = [8 12 16]; %all the deadlines
N = 10; %number of bits

Poma=0; % max power needed by OMA
for m =1 : M
    if m >1
        Poma(m) = (exp(N/(Dall(m) - Dall(m-1)))-1)/hall(m);
    else
        Poma(m) = (exp(N/(Dall(m) ))-1)/hall(m);
    end
end

tm(1) = Dall(1);% User 1's deadline, t1
P = zeros(M,M); 
P(1,1) = (exp(N/tm(1))-1)/hall(1); %User 1's transmit power
am1(1) = 1/(1+P(1,1)*hall(1)); % a21 is decided by User 1's parameters
for m = 2: M
    tm(m) = Dall(m) - Dall(m-1);
    if tm(m)<=-N/log(am1(m-1))  %hybrid NOMA
        Pm1 = (am1(m-1)*exp( (N-Dall(m-1)*log(am1(m-1)))/(Dall(m))  )-1)/am1(m-1)/hall(m);
        Pmm = ( exp((N-Dall(m-1)*log(am1(m-1)))/(Dall(m))  )-1)/hall(m);
        P(m,1:m-1) = Pm1;
        P(m,m) = Pmm;
    else %OMA
        P(m,m) =  (exp(N/tm(m))-1)/hall(m);
    end
    am1(m) = 1/(1 + sum(P(1:m,1).*hall(1:m).'));
end 

oma = tm(1)*P(1,1);;
ojbe = tm(1)*P(1,1);
for m = 2: M
    oma = oma + tm(m)*Poma(m);
    for n = 1 : m
        ojbe = ojbe + tm(n)*P(m,n);
    end
end


P11 = P(1,1); P21 = P(2,1); P22 = P(2,2);
t1 = Dall(1); t2 = Dall(2) - Dall(1); t3 = Dall(3) - Dall(2);
h1 = hall(1); h2 = hall(2); h3 = hall(3);
Pvec = [0:0.1: 5];
for i = 1 : length(Pvec)
    for j = 1: length(Pvec)
            P31 = Pvec(i); P32 = Pvec(j);  
            nbit_2slot = t1*log(1+P31*h3/(1+P11*h1+P21*h2))...
                +t2*log(1+P32*h3/(1+P22*h2));
            P33x = (exp( (N-nbit_2slot)/t3 )-1)/h3;
            if P33x>0 & P33x<Poma(3)
                es(i,j) = t1*P31 + t2*P32 + t3*P33x;                
            else
                es(i,j) = inf;           
            end
    end
end

obj3 = t1*P(3,1) + t2*P(3,2) + t3*P(3,3);

surf(Pvec,Pvec,es) 

hold

surf(Pvec, Pvec, obj3*ones(length(Pvec),length(Pvec)))
%alpha 0.3
            

