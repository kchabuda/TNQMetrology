function [result,resultm,a0,c] = MPO_F_PBC_Inf_opt
figureofmerit = 1; % 1 = 2Tr(r'L)-Tr(rLL), 2 = ||2*r'-rL-Lr||^2/||2*r'||^2
d = 2;
noiserange = 2;
integral0 = 1;
integral1 = 0.1;
integral2 = 0.01;
lherm = 1; % 1 = yes, /else/ = no
phi = 10^-4;
imprecision = 10^(-5/2);
ratio = 10^(-5/2);
bdpsimax = 100;
bdlmax = 100;
%rng('default')
resultm = zeros([1,100]);
bdpsi = 1;
bdl = 1;
a0 = sqrt(2/(d+1))*sin((1:d)*pi/(d+1));
a0 = permute(a0,[1,3,2]);
c = triu(ones(d)-eye(d));
c = 1i*phi*exp(-integral0/2)*c/(1+2*exp(-integral0)*sinh(integral1));
c = c+c';
c = eye(d)+c;
c = reshape(c,[bdl,bdl,d,d]);
[resultm(bdl,bdpsi),a0,c] = MPO_F_PBC_Inf(figureofmerit,d,bdpsi,bdl,noiserange,integral0,integral1,integral2,lherm,phi,imprecision,a0,c);
while 1
    while 1
        if bdpsi == bdpsimax
            break
        else
            a0old = a0;
            bdpsi = bdpsi+1;
            a0 = zeros([bdpsi,bdpsi,d]);
            for i = 1:d
                if i <= ceil(d/2)
                    a0oldihalf = triu(rot90(a0old(:,:,i),-1));
                    a0(1:bdpsi-1,2:bdpsi,i) = a0oldihalf;
                    a0(:,:,i) = a0(:,:,i)+a0(:,:,i).';
                    a0(:,:,i) = a0(:,:,i)+diag([0;diag(a0(:,:,i),2);0]);
                    a0(:,:,i) = rot90(a0(:,:,i),1);
                    a0(1,bdpsi,i) = ratio*(1+1i)*abs(a0(1,bdpsi-1,i));
                    a0(bdpsi,1,i) = conj(a0(1,bdpsi,i));
                    if i == ceil(d/2) && mod(d,2) == 1
                        a0(:,:,i) = (a0(:,:,i)+a0(:,:,i).')/2;
                    end
                else
                    a0(:,:,i) = a0(:,:,d+1-i).';
                end
            end
            tensors = {conj(a0),a0};
            legs = {[-1,-3,1],[-2,-4,1]};
            tm = ncon(tensors,legs);
            tm = reshape(tm,[bdpsi*bdpsi,bdpsi*bdpsi]);
            a0norm = eigs(tm,1);
            a0norm = abs(a0norm)^(1/2);
            a0 = a0/a0norm;
            [resultm(bdl,bdpsi),a0new,cnew] = MPO_F_PBC_Inf(figureofmerit,d,bdpsi,bdl,noiserange,integral0,integral1,integral2,lherm,phi,imprecision,a0,c);
            if resultm(bdl,bdpsi) < (1+imprecision)*resultm(bdl,bdpsi-1)
                bdpsi = bdpsi-1;
                a0 = a0old;
                a0copy = a0new;
                ccopy = cnew;
                break
            else
                a0 = a0new;
                c = cnew;
            end
        end
    end
    if bdl == bdlmax
        if bdpsi == bdpsimax
            resultm = resultm(1:bdl,1:bdpsi);
            result = resultm(bdl,bdpsi);
        else
            a0 = a0copy;
            c = ccopy;
            resultm = resultm(1:bdl,1:bdpsi+1);
            result = resultm(bdl,bdpsi+1);
        end
        break
    else
        cold = c;
        bdl = bdl+1;
        factor = 1/2;
        while 1
            c = zeros([bdl,bdl,d,d]);
            for nx = 1:d
                for nxp = 1:d
                    if nx ~= nxp
                        meanrecold = sum(sum(abs(real(cold(:,:,nx,nxp)))))/(bdl-1)^2;
                        meanimcold = sum(sum(abs(imag(cold(:,:,nx,nxp)))))/(bdl-1)^2;
                        c(:,:,nx,nxp) = (meanrecold*rand(bdl)+1i*meanimcold*rand(bdl))*factor;
                    end
                end
            end
            c = (c+conj(permute(c,[1,2,4,3])))/2;
            c(1:bdl-1,1:bdl-1,:,:,:) = cold;
            tensors = {c};
            legs = {[-1,-2,1,1]};
            tm = ncon(tensors,legs);
            ctr = eigs(tm,1);
            ctr = real(ctr);
            c = d*c/ctr;
            [resultm(bdl,bdpsi),a0new,cnew] = MPO_F_PBC_Inf(figureofmerit,d,bdpsi,bdl,noiserange,integral0,integral1,integral2,lherm,phi,imprecision,a0,c);
            if resultm(bdl,bdpsi) > resultm(bdl-1,bdpsi)
                a0 = a0new;
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
                    a0 = a0copy;
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
if noiserange == 1
    eta = exp(-integral0/2);
    boundloc = eta^2/(1-eta^2);
elseif noiserange == 2
    bound05 = Bound_F(1/2,d,integral0,integral1);
    bound1 = Bound_F(1,d,integral0,integral1);
    bound2 = Bound_F(2,d,integral0,integral1);
    bound3 = Bound_F(3,d,integral0,integral1);
end
%}
%{
bdpsi = size(a0,1);
tensors = {conj(a0),a0};
legs = {[-1,-3,1],[-2,-4,1]};
tm = ncon(tensors,legs);
tm = reshape(tm,[bdpsi*bdpsi,bdpsi*bdpsi]);
[tmvr,~] = eigs(tm,1);
[tmvl,~] = eigs(tm.',1);
tmvnorm = tmvl.'*tmvr;
r = reshape(tmvr/tmvnorm^(1/2),[bdpsi,bdpsi]);
r = r+r';
l = reshape(tmvl.'/tmvnorm^(1/2),[bdpsi,bdpsi]);
l = l+l';
rhoreduced = sqrtm(l)*r*sqrtm(l);
rhoreduced = rhoreduced+rhoreduced';
rhoreduced = rhoreduced/trace(rhoreduced);
rhoreducedeigval = eig(rhoreduced);
entropy = -rhoreducedeigval.'*log2(rhoreducedeigval);
%}
end