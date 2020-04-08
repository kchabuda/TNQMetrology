alfa = 1; gamma = 2; beta = 0.1;
tioptinf = 1.17334774934831;
foptinf = 0.510405826179561;
mlist = [3:20].';
talist = logspace(-1,2,100).';
avar = zeros(length(talist),1);
qavar = zeros(length(talist),length(mlist));
qavarinf = zeros(length(talist),1);
for j = 1:length(mlist)
    for i = 1:length(talist)
        avar(i,1) = sigma0(talist(i),alfa,gamma,beta);
        if (mlist(j) == 2 || 0.5*tioptinf*mlist(j) < talist(i)) && talist(i) < 2*tioptinf*mlist(j)
            [~,qavar(i,j),~,~] = MPO_A_opt(mlist(j),talist(i));
        else
            qavar(i,j) = avar(i,1);
        end
        qavarinf(i,1) = (2*alfa/gamma+beta-foptinf)/talist(i);
    end
end

function y = sigma0(t,alfa,gamma,beta)
    y = alfa*(2*gamma*t+4*exp(-gamma*t)-exp(-2*gamma*t)-3)/(gamma*t)^2+beta/t;
end