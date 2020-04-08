imprecision = 10^-2;
mlist = [2:9,unique(round(logspace(1,3,20)))].';
mlen = length(mlist);
tilist = zeros([mlen,1]);
flist = zeros([mlen,1]);
for i = 1:mlen
    fun = @(ti)-MPO_A_opt(mlist(i),ti);
    options = optimset('Display','iter','TolX',0.1*imprecision);
    if i == 1
        [tiopt,fopt] = fminbnd(fun,0,10,options);
        tilist(i) = tiopt;
        flist(i) = -fopt;
        [tiopt,fopt] = fminbnd(fun,0.5*tilist(i),1.1*tilist(i),options);
        tilist(i) = tiopt;
        flist(i) = -fopt;
    else
        [tiopt,fopt] = fminbnd(fun,0.5*tilist(i-1),1.1*tilist(i-1),options);
        tilist(i) = tiopt;
        flist(i) = -fopt;
    end
end