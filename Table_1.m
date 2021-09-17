clear all
close all
M=3; %number of users
hall = ones(1,M); %all the channels
Dall = [8 12 16]; %all the deadlines
N = 12; %number of bits
w = [1 1 1]/3;%[1 2 3]/6;   you need to change this for different weighting coefficients

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

oma = tm(1)*P(1,1);
ojbe = tm(1)*P(1,1);
for m = 2: M
    oma = oma + tm(m)*Poma(m);
    for n = 1 : m
        ojbe = ojbe + tm(n)*P(m,n);
    end
end

%three user
P11 = P(1,1); P21 = P(2,1);P22 = P(2,2); P31 = P(3,1);P32=P(3,2);P33=P(3,3);
t1 = Dall(1); t2 = Dall(2) - Dall(1); t3 = Dall(3) - Dall(2);
ojbe =  w(1)*t1*P11 + w(2)*(t1*P21+t2*P22) +w(3)*(t1*P31+t2*P32+t3*P33);

P11 = P(1,1); P21 = P(2,1); P22 = P(2,2);
t1 = Dall(1); t2x = Dall(2) - Dall(1); t3x = Dall(3) - Dall(2);
h1 = hall(1); h2 = hall(2); h3 = hall(3);
Pvec = [0:0.1: 10];
t3vec = [3.01:0.5:Dall(3) - Dall(2) Dall(3) - Dall(2)];
t2vec = [3.01:0.5:Dall(2) - Dall(1) Dall(2) - Dall(1)];

all_vec=[0 0 0 0 0 0 0 0];
energy=100000;
index1=1;
for i2 = 1 : length(Pvec)
    for j2 = 1 : length(t2vec)
        P21 = Pvec(i2); t2 = t2vec(j2);
        nbit_2slot = t1*log(1+P21*h2/(1+P11*h1));
        P22 = (exp( (N-nbit_2slot)/t2 )-1)/h2;
        for i3 = 1 : length(Pvec)
            for j3 = 1: length(Pvec)
                for k3 = 1: length(t3vec)
                    P31 = Pvec(i3); P32 = Pvec(j3);  t3 = t3vec(k3);
                    nbit_3slot = t1*log(1+P31*h3/(1+P11*h1+P21*h2))...
                        +t2*log(1+P32*h3/(1+P22*h2));% + t3*log(1+P33*h3);
                    P33 = (exp( (N-nbit_3slot)/t3 )-1)/h3;
                    x=w(1)*t1*P11 + w(2)*(t1*P21+t2*P22) +w(3)*(t1*P31+t2*P32+t3*P33);
                    if P22>=0 & P33>=0 & P22<Poma(2) & P33<Poma(3) & x<=energy
                        all_vec = [t2 P21 P22 t3 P31 P32 P33 x];
                        energy = x;
                        %all_vec(index1,end) = t1*P11 + t1*P21+t2*P22+t1*P31+t2*P32+t3*P33;
                        index1 = index1+1;
                    end
                end
            end
        end
    end
end


[energy ojbe energy-ojbe]'

