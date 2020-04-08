function [tilist,flist] = Exact_A_Inf_N00N
alfa = 1; gamma = 2; beta = 0.1;
function y = c0(ti)
    y = 2*alfa*(gamma*ti+exp(-gamma*ti)-1)/gamma^2+beta*ti;
end
function y = corr(dn,ti)
    y = c0(ti)^2*(dn-1)^2*exp(-c0(ti)*(dn-1)^2)/ti;
end
imprecision = 10^-3;
dnlist = 2:11;
dnlist = dnlist.';
dnlen = length(dnlist);
tilist = zeros([dnlen,1]);
flist = zeros([dnlen,1]);
for i = 1:dnlen
    fun = @(ti)-corr(dnlist(i),ti);
    options = optimset('Display','iter','TolX',0.1*imprecision);
    if i == 1
        [tiopt,fopt] = fminbnd(fun,0,10);
        tilist(i) = tiopt;
        flist(i) = -fopt;
        [tiopt,fopt] = fminbnd(fun,0.1*tilist(i),1.1*tilist(i),options);
        tilist(i) = tiopt;
        flist(i) = -fopt;
    else
        [tiopt,fopt] = fminbnd(fun,0.1*tilist(i-1),1.1*tilist(i-1),options);
        tilist(i) = tiopt;
        flist(i) = -fopt;
    end
end
end