imprecision = 10^-3;
mlist = 2:7;
mlist = mlist.';
mlen = length(mlist);
tilist = zeros([mlen,1]);
flist = zeros([mlen,1]);
for i = 1:mlen
    fun = @(ti)-Exact_A_OBC_prod(mlist(i),ti);
    %fun = @(ti)-Exact_A_OBC(mlist(i),ti);
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