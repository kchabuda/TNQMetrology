imprecision = 10^(-5/2);
dlist = 2:11;
dlist = dlist.';
dlen = length(dlist);
tilist = zeros([dlen,1]);
flist = zeros([dlen,1]);
for i = 1:dlen
    fun = @(ti)-MPO_A_PBC_Inf_opt(dlist(i),ti);
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