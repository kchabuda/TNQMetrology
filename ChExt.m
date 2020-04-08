function c = ChExt(kraus, dkraus)
%returns the asymptotic bound coeeficient c, (c/sqrt(N)) given 
%the set of linearly independent Kraus operators and their parameter
%derivatives. kraus, dkraus should be a cell arrays with kraus{i} being 
%Kraus matrix, etc
n = length(kraus); %number of Kraus operators
dims=size(kraus{1}); %dimensions of the kraus matrix (number of rows, number of columns),
% (dimension of the output system, dimension of the input system
din=dims(2);
dout=dims(1); %dimensions of the output and input systems
kc_dk = zeros(din,din);
for j = 1:n
   kc_dk=kc_dk+kraus{j}'*dkraus{j};
end
cvx_begin quiet
    variable h(n,n) hermitian;
    variable t;
    kc_h_k = zeros(din,din);
    for j = 1:n
        for k = 1:n
            kc_h_k=kc_h_k+kraus{j}'*kraus{k}*h(j,k);
        end
    end
    dd=din+n*dout;
    a = t*eye(dd,dd);
    for j = 1:n
        w{j}=-i*dkraus{j};
        for k = 1:n
            w{j} = w{j} + h(j,k)*kraus{k};
        end        
    end
    for j = 1:n
        for k = 1 : dout
            for l = 1 : din
                a(din + (j-1)*dout+k,l)=w{j}(k,l);
                a(l,din+(j-1)*dout+k)=w{j}(k,l)';
            end
        end    
    end
    minimize(t);
    subject to
        i*kc_dk == kc_h_k;
        a == semidefinite(dd,dd);  
cvx_end
h;
%opnorm = norm(a);
c=1/(2*t);
end