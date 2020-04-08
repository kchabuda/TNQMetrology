numof2box = 1;
d = 2;
int0 = 1;
int1 = -0.49:0.01:0.49;
tab = zeros([length(int1),2]);
for i = 1:length(int1)
    fun = @(int0p)Bound_F_general(numof2box,d,int0,int1(i),int0p);
    [tab(i,1),tab(i,2)] = fminbnd(fun,0,int0/2);
end