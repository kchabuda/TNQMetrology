function [result,qAvar,resultm,A0,c] = MPO_A_opt
approx = 0;  % 1 = yes, /else/ = no
figureofmerit = 1; % 1 = 2Tr(r'L)-Tr(rLL), 2 = ||2*r'-rL-Lr||^2/||2*r'||^2
m = 5;
if approx == 1
    n = m-1;
else
    n = 2*m-1;
end
d = 2;
ti = 1;
noiserange = 2;
alfa = 1; gamma = 2; beta = 0.1;
lherm = 1; % 1 = yes, /else/ = no
imprecision = 10^-2;
bdpsimax = 1;
bdlmax = 100;
%rng('default');
resultm = zeros(100);
function y = sigma0(t,alfa,gamma,beta)
    y = alfa*(2*gamma*t+4*exp(-gamma*t)-exp(-2*gamma*t)-3)/(gamma*t)^2+beta/t;
end
integral0 = 2*alfa*(gamma*ti+exp(-gamma*ti)-1)/gamma^2+beta*ti;
integral1 = 2*alfa*(cosh(gamma*ti)-1)*exp(-1*gamma*ti)/gamma^2;
integral2 = 2*alfa*(cosh(gamma*ti)-1)*exp(-2*gamma*ti)/gamma^2;
bdpsi = 1;
bdl = 2;
a0 = sqrt(2/(d+1))*sin((1:d)*pi/(d+1));
a0 = permute(a0,[1,3,2]);
A0 = cell([1,n]);
[A0{:}] = deal(a0);
store = cell([3,2]);
for rep = 1:2
    c = (rand([bdl,bdl,d,d,n])+1i*rand([bdl,bdl,d,d,n]))/bdl;
    c = (c+conj(permute(c,[1,2,4,3,5])))/2;
    c(2:bdl,:,:,:,1) = zeros([bdl-1,bdl,d,d]);
    c(:,2:bdl,:,:,n) = zeros([bdl,bdl-1,d,d]);
    [f,A0new,cnew] = MPO_A_OBC(approx,figureofmerit,m,d,bdl,ti,noiserange,integral0,integral1,integral2,lherm,imprecision,A0,c);
    store{1,rep} = f;
    store{2,rep} = A0new;
    store{3,rep} = cnew;
end
[resultm(bdl,bdpsi),position] = max(cell2mat(store(1,:)));
A0 = store{2,position};
c = store{3,position};
while 1
    while 1
        if bdpsi == bdpsimax
            break
        else
            A0old = A0;
            bdpsi = bdpsi+1;
            factor = 1/2;
            while 1
                A0 = cell([1,n]);
                [s1,s2,~] = size(A0old{1});
                for nx = 1:d
                    meanrea0old = sum(sum(abs(real(A0old{1}(:,:,nx)))))/(s1*s2);
                    meanima0old = sum(sum(abs(imag(A0old{1}(:,:,nx)))))/(s1*s2);
                    A0{1}(:,:,nx) = (meanrea0old*rand([s1,s2+1])+1i*meanima0old*rand([s1,s2+1]))*factor;
                end
                A0{1}(1:s1,1:s2,:) = A0old{1};
                for x = 2:n-1
                    [s1,s2,~] = size(A0old{x});
                    for nx = 1:d
                        meanrea0old = sum(sum(abs(real(A0old{x}(:,:,nx)))))/(s1*s2);
                        meanima0old = sum(sum(abs(imag(A0old{x}(:,:,nx)))))/(s1*s2);
                        A0{x}(:,:,nx) = (meanrea0old*rand([s1+1,s2+1])+1i*meanima0old*rand([s1+1,s2+1]))*factor;
                    end
                    A0{x}(1:s1,1:s2,:) = A0old{x};
                end
                [s1,s2,~] = size(A0old{n});
                for nx = 1:d
                    meanrea0old = sum(sum(abs(real(A0old{n}(:,:,nx)))))/(s1*s2);
                    meanima0old = sum(sum(abs(imag(A0old{n}(:,:,nx)))))/(s1*s2);
                    A0{n}(:,:,nx) = (meanrea0old*rand([s1+1,s2])+1i*meanima0old*rand([s1+1,s2]))*factor;
                end
                A0{n}(1:s1,1:s2,:) = A0old{n};
                tensors = {conj(A0{n}),A0{n}};
                legs = {[-1,-3,1],[-2,-4,1]};
                r1 = ncon(tensors,legs);
                A0{n} = A0{n}/norm(r1(:))^(1/2);
                tensors = {conj(A0{n}),A0{n}};
                legs = {[-1,-3,1],[-2,-4,1]};
                r2 = ncon(tensors,legs);
                for x = n-1:-1:2
                    tensors = {conj(A0{x}),A0{x},r2};
                    legs = {[-1,2,1],[-2,3,1],[2,3]};
                    r1 = ncon(tensors,legs);
                    A0{x} = A0{x}/norm(r1(:))^(1/2);
                    tensors = {conj(A0{x}),A0{x},r2};
                    legs = {[-1,2,1],[-2,3,1],[2,3]};
                    r2 = ncon(tensors,legs);
                end
                tensors = {conj(A0{1}),A0{1},r2};
                legs = {[-1,2,1],[-2,3,1],[2,3]};
                r1 = ncon(tensors,legs);
                A0{1} = A0{1}/norm(r1(:))^(1/2);
                [resultm(bdl,bdpsi),A0new,cnew] = MPO_A_OBC(approx,figureofmerit,m,d,bdl,ti,noiserange,integral0,integral1,integral2,lherm,imprecision,A0,c);
                if resultm(bdl,bdpsi) > resultm(bdl,bdpsi-1)
                    break
                end
                factor = factor/2;
            end
            if resultm(bdl,bdpsi) < (1+imprecision)*resultm(bdl,bdpsi-1)
                bdpsi = bdpsi-1;
                A0 = A0old;
                A0copy = A0new;
                ccopy = cnew;
                break
            else
                A0 = A0new;
                c = cnew;
            end
        end
    end
    if bdl == bdlmax
        if bdpsi == bdpsimax
            resultm = resultm(1:bdl,1:bdpsi);
            result = resultm(bdl,bdpsi);
        else
            A0 = A0copy;
            c = ccopy;
            resultm = resultm(1:bdl,1:bdpsi+1);
            result = resultm(bdl,bdpsi+1);
        end
        break
    else
        cold = c;
        bdl = bdl+1;
        factor = 1/10;
        while 1
            c = zeros([bdl,bdl,d,d,n]);
            for nx = 1:d
                for nxp = 1:d
                    meanrecold = sum(sum(abs(real(cold(:,:,nx,nxp,1)))))/(bdl-1);
                    meanimcold = sum(sum(abs(imag(cold(:,:,nx,nxp,1)))))/(bdl-1);
                    c(1,:,nx,nxp,1) = (meanrecold*rand([1,bdl])+1i*meanimcold*rand([1,bdl]))*factor;
                    for x = 2:n-1
                        meanrecold = sum(sum(abs(real(cold(:,:,nx,nxp,x)))))/(bdl-1)^2;
                        meanimcold = sum(sum(abs(imag(cold(:,:,nx,nxp,x)))))/(bdl-1)^2;
                        c(:,:,nx,nxp,x) = (meanrecold*rand(bdl)+1i*meanimcold*rand(bdl))*factor;
                    end
                    meanrecold = sum(sum(abs(real(cold(:,:,nx,nxp,n)))))/(bdl-1);
                    meanimcold = sum(sum(abs(imag(cold(:,:,nx,nxp,n)))))/(bdl-1);
                    c(:,1,nx,nxp,n) = (meanrecold*rand([bdl,1])+1i*meanimcold*rand([bdl,1]))*factor;
                end
            end
            c = (c+conj(permute(c,[1,2,4,3,5])))/2;
            c(1:bdl-1,1:bdl-1,:,:,:) = cold;
            [resultm(bdl,bdpsi),A0new,cnew] = MPO_A_OBC(approx,figureofmerit,m,d,bdl,ti,noiserange,integral0,integral1,integral2,lherm,imprecision,A0,c);
            if resultm(bdl,bdpsi) > resultm(bdl-1,bdpsi)
                A0 = A0new;
                c = cnew;
                break
            end
            factor = factor/2;
        end
        if resultm(bdl,bdpsi) < (1+imprecision)*resultm(bdl-1,bdpsi)
            if bdpsi == bdpsimax
                resultm = resultm(1:bdl,1:bdpsi);
                result = resultm(bdl,bdpsi);
            else
                if resultm(bdl,bdpsi) < resultm(bdl-1,bdpsi+1)
                    A0 = A0copy;
                    c = ccopy;
                    resultm = resultm(1:bdl,1:bdpsi+1);
                    result = resultm(bdl-1,bdpsi+1);
                else
                    resultm = resultm(1:bdl,1:bdpsi+1);
                    result = resultm(bdl,bdpsi);
                end
            end
            break
        end
    end
end
%{
fullnoise = 0; % 1 = yes, /else/ = no
if bdpsimax == 1
    [exact,exactqAvar] = Exact_A_OBC_prod(m,d,ti,fullnoise,noiserange,alfa,gamma,beta);
else
    [exact,exactqAvar] = Exact_A_OBC(m,d,ti,fullnoise,noiserange,alfa,gamma,beta,imprecision);
end
%}
qAvar = sigma0(ti*m,alfa,gamma,beta)-result/(ti*m);
end