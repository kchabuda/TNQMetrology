function [result,resultm,c] = MPO_Fid_opt(boundary,n,d,A,B,phi)
%boundary = 'O'; % 'O', 'P'
figureofmerit = 1; % 1 = 2Tr(r'L)-Tr(rLL), 2 = ||2*r'-rL-Lr||^2/||2*r'||^2
%n = 10;
%d = 2;
lherm = 1; % 1 = yes, /else/ = no
imprecision = 10^-2;
bdlmax = 100;
%rng('default');
resultm = zeros([100,1]);
if boundary == 'O'
    tensors = {A{n}};
    legs = {[-1,-2,1,1]};
    Atr = ncon(tensors,legs);
    for x = n-1:-1:2
        tensors = {A{x},Atr};
        legs = {[-1,2,1,1],[2]};
        Atr = ncon(tensors,legs);
    end
    tensors = {A{1},Atr};
    legs = {[-1,2,1,1],[2]};
    Atr = ncon(tensors,legs);
    Atr = abs(Atr)^(1/n);
    A = cellfun(@(x) x/Atr,A,'UniformOutput',false);
    tensors = {B{n}};
    legs = {[-1,-2,1,1]};
    Btr = ncon(tensors,legs);
    for x = n-1:-1:2
        tensors = {B{x},Btr};
        legs = {[-1,2,1,1],[2]};
        Btr = ncon(tensors,legs);
    end
    tensors = {B{1},Btr};
    legs = {[-1,2,1,1],[2]};
    Btr = ncon(tensors,legs);
    Btr = abs(Btr)^(1/n);
    B = cellfun(@(x) x/Btr,B,'UniformOutput',false);
    bdl = 2;
    store = cell([2,2]);
    for rep = 1:2
        c = (rand([bdl,bdl,d,d,n])+1i*rand([bdl,bdl,d,d,n]))/bdl;
        c = (c+conj(permute(c,[1,2,4,3,5])))/2;
        c(2:bdl,:,:,:,1) = zeros([bdl-1,bdl,d,d]);
        c(:,2:bdl,:,:,n) = zeros([bdl,bdl-1,d,d]);
        [f,cnew] = MPO_Fid_OBC(figureofmerit,n,d,bdl,lherm,imprecision,A,B,phi,c);
        store{1,rep} = f;
        store{2,rep} = cnew;
    end
    [resultm(bdl),position] = max(cell2mat(store(1,:)));
    c = store{2,position};
    while 1
        if bdl == bdlmax
            resultm = resultm(1:bdl);
            result = resultm(bdl);
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
                [resultm(bdl),cnew] = MPO_Fid_OBC(figureofmerit,n,d,bdl,lherm,imprecision,A,B,phi,c);
                if resultm(bdl) > resultm(bdl-1)
                    c = cnew;
                    break
                end
                factor = factor/2;
            end
            if resultm(bdl) < (1+imprecision)*resultm(bdl-1)
                resultm = resultm(1:bdl);
                result = resultm(bdl);
                break
            end
        end
    end
    %exact = Exact_Fid_OBC(n,d,A,B,phi);
elseif boundary == 'P'
    tensors = {A(:,:,:,:,n)};
    legs = {[-1,-2,1,1]};
    Atr = ncon(tensors,legs);
    for x = n-1:-1:2
        tensors = {A(:,:,:,:,x),Atr};
        legs = {[-1,2,1,1],[2,-2]};
        Atr = ncon(tensors,legs);
    end
    tensors = {A(:,:,:,:,1),Atr};
    legs = {[3,2,1,1],[2,3]};
    Atr = ncon(tensors,legs);
    Atr = abs(Atr)^(1/n);
    A = A/Atr;
    tensors = {B(:,:,:,:,n)};
    legs = {[-1,-2,1,1]};
    Btr = ncon(tensors,legs);
    for x = n-1:-1:2
        tensors = {B(:,:,:,:,x),Btr};
        legs = {[-1,2,1,1],[2,-2]};
        Btr = ncon(tensors,legs);
    end
    tensors = {B(:,:,:,:,1),Btr};
    legs = {[3,2,1,1],[2,3]};
    Btr = ncon(tensors,legs);
    Btr = abs(Btr)^(1/n);
    B = B/Btr;
    bdl = 2;
    store = cell([2,2]);
    for rep = 1:2
        c = (rand([bdl,bdl,d,d,n])+1i*rand([bdl,bdl,d,d,n]))/bdl;
        c = (c+conj(permute(c,[1,2,4,3,5])))/2;
        [f,cnew] = MPO_Fid_PBC(figureofmerit,n,d,bdl,lherm,imprecision,A,B,phi,c);
        store{1,rep} = f;
        store{2,rep} = cnew;
    end
    [resultm(bdl),position] = max(cell2mat(store(1,:)));
    c = store{2,position};
    while 1
        if bdl == bdlmax
            resultm = resultm(1:bdl);
            result = resultm(bdl);
            break
        else
            cold = c;
            bdl = bdl+1;
            factor = 1/10;
            while 1
                c = zeros([bdl,bdl,d,d,n]);
                for nx = 1:d
                    for nxp = 1:d
                        for x = 1:n
                            meanrecold = sum(sum(abs(real(cold(:,:,nx,nxp,x)))))/(bdl-1)^2;
                            meanimcold = sum(sum(abs(imag(cold(:,:,nx,nxp,x)))))/(bdl-1)^2;
                            c(:,:,nx,nxp,x) = (meanrecold*rand(bdl)+1i*meanimcold*rand(bdl))*factor;
                        end
                    end
                end
                c = (c+conj(permute(c,[1,2,4,3,5])))/2;
                c(1:bdl-1,1:bdl-1,:,:,:) = cold;
                [resultm(bdl),cnew] = MPO_Fid_PBC(figureofmerit,n,d,bdl,lherm,imprecision,A,B,phi,c);
                if resultm(bdl) > resultm(bdl-1)
                    c = cnew;
                    break
                end
                factor = factor/2;
            end
            if resultm(bdl,bdpsi) < (1+imprecision)*resultm(bdl-1,bdpsi)
                resultm = resultm(1:bdl);
                result = resultm(bdl);
                break
            end
        end
    end
    %exact = Exact_Fid_PBC(n,d,A,B,phi);
end
fids = 1/4*result;
fid = 1-1/8*result*phi^2;
end