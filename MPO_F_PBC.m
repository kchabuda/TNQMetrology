function [result,a0,c] = MPO_F_PBC(figureofmerit,n,d,bdpsi,bdl,noiserange,integral0,integral1,integral2,lherm,imprecision,a0,c)
% Parameters
tol1 = 0.1*imprecision/n^2;
tol2 = 0.1*imprecision/n^2;
told = 0.1*imprecision/n^2;
relunc = 0.1*imprecision;
relunc1 = 0.1*imprecision;
relunc2 = imprecision;
reluncd = 0.1*imprecision;
global ncon_skipCheckInputs;
ncon_skipCheckInputs = true;
% Functions
function y = pinv2(x,tol)
    [uf,sf,vf] = svd(x);
    sf = diag(sf);
    positionf = sf > sf(1)*tol;
    sf2 = sf(positionf).^-1;
    sf2(end+1:numel(sf)) = 0;
    sf2 = diag(sf2);
    y = vf*sf2*uf';
end
% Calculations
%  MPO
if noiserange == 0
    bdnoise = 1;
	an = ones([1,1,d,d]);
elseif noiserange == 1
    nxpm = repmat(1:d,[d,1]);
    nxm = nxpm.';
    nxminusnxp = nxm-nxpm;
    bdnoise = 1;
    an = exp(-nxminusnxp.^2*integral0/2);
    an = permute(an,[3,4,1,2]);
elseif noiserange == 2
    nxpm = repmat(1:d,[d,1]);
    nxm = nxpm.';
    nxminusnxp = nxm-nxpm;
    bdnoise1 = 1;
    an1 = exp(-nxminusnxp.^2*integral0/2);
    an1 = permute(an1,[3,4,1,2]);
    bdnoise2 = 2*d-1;
    if n == 2
        corr12 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral1/2);
    else
        corr12 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral1);
    end
    [u,s,v] = svds(corr12,bdnoise2);
    s = diag(s);
    us = u*diag(s.^(1/2));
    sv = diag(s.^(1/2))*v';
    an2 = zeros([bdnoise2,bdnoise2,d,d]);
    for nx = 1:d
        for nxp = 1:d
            an2(:,:,nx,nxp) = sv(:,d*(nx-1)+nxp)*us(d*(nx-1)+nxp,:);
        end
    end
    bdnoise = bdnoise1*bdnoise2;
    an = zeros([bdnoise,bdnoise,d,d]);
    for nx = 1:d
        for nxp = 1:d
            an(:,:,nx,nxp) = an1(:,:,nx,nxp)*an2(:,:,nx,nxp);
        end
    end
elseif noiserange == 3
    nxpm = repmat(1:d,[d,1]);
    nxm = nxpm.';
    nxminusnxp = nxm-nxpm;
    bdnoise1 = 1;
    an1 = exp(-nxminusnxp.^2*integral0/2);
    an1 = permute(an1,[3,4,1,2]);
    bdnoise2 = 2*d-1;
    if n == 2
        corr12 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral1/2);
    else
        corr12 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral1);
    end
    [u,s,v] = svds(corr12,bdnoise2);
    s = diag(s);
    us = u*diag(s.^(1/2));
    sv = diag(s.^(1/2))*v';
    an2 = zeros([bdnoise2,bdnoise2,d,d]);
    for nx = 1:d
        for nxp = 1:d
            an2(:,:,nx,nxp) = sv(:,d*(nx-1)+nxp)*us(d*(nx-1)+nxp,:);
        end
    end
    if n == 2 || n == 3
        bdnoise3 = 1;
        an3 = ones([1,1,d,d]);
    else
        bdnoise3ini = 2*d-1;
        bdnoise3 = bdnoise3ini^2;
        if n == 4
            corr13 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral2/2);
        else
            corr13 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral2);
        end
        [u,s,v] = svds(corr13,bdnoise3ini);
        s = diag(s);
        us = u*diag(s.^(1/2));
        sv = diag(s.^(1/2))*v';
        an3 = zeros([bdnoise3ini,bdnoise3ini,d,d]);
        for nx = 1:d
            for nxp = 1:d
                an3(:,:,nx,nxp) = sv(:,d*(nx-1)+nxp)*us(d*(nx-1)+nxp,:);
            end
        end
        tensors = {an3,eye(bdnoise3ini)};
        legs = {[-1,-4,-5,-6],[-2,-3]};
        an3 = ncon(tensors,legs);
        an3 = reshape(an3,[bdnoise3,bdnoise3,d,d]);
    end
    bdnoise = bdnoise1*bdnoise2*bdnoise3;
    an = zeros([bdnoise,bdnoise,d,d]);
    for nx = 1:d
        for nxp = 1:d
            an(:,:,nx,nxp) = an1(:,:,nx,nxp)*kron(an2(:,:,nx,nxp),an3(:,:,nx,nxp));
        end
    end
end
f_iter = zeros([2*100,1]);
iter = 0;
while 1
    iter = iter+1;
    a = zeros([bdpsi^2*bdnoise,bdpsi^2*bdnoise,d,d,n]);
    for nx = 1:d
        for nxp = 1:d
            for x = 1:n
                a(:,:,nx,nxp,x) = kron(kron(a0(:,:,nx,x),conj(a0(:,:,nxp,x))),an(:,:,nx,nxp));
            end
        end
    end
    b = zeros([2*bdpsi^2*bdnoise,2*bdpsi^2*bdnoise,d,d,n]);
    for nx = 1:d
        for nxp = 1:d
            b(:,:,nx,nxp,1) = kron([-1i*(nx-nxp),1;0,0],a(:,:,nx,nxp,1));
            for x = 2:n-1
                b(:,:,nx,nxp,x) = kron([1,0;-1i*(nx-nxp),1],a(:,:,nx,nxp,x));
            end
            b(:,:,nx,nxp,n) = kron([1,0;-1i*(nx-nxp),0],a(:,:,nx,nxp,n));
        end
    end
%  2Tr(r'L)-Tr(rLL) maximization
    if figureofmerit == 1
        l1f = zeros([bdl,2*bdpsi^2*bdnoise,bdl,2*bdpsi^2*bdnoise,n-1]);
        l2f = zeros([bdl,bdpsi^2*bdnoise,bdl,bdl,bdpsi^2*bdnoise,bdl,n-1]);
        fom1 = zeros([n*100,1]);
        index = 0;
        iter1 = 0;
        while 1
            iter1 = iter1+1;
            tensors = {c(:,:,:,:,n),b(:,:,:,:,n)};
            legs = {[-1,-3,1,2],[-2,-4,2,1]};
            l1f(:,:,:,:,n-1) = ncon(tensors,legs);
            tensors = {c(:,:,:,:,n),a(:,:,:,:,n),c(:,:,:,:,n)};
            legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
            l2f(:,:,:,:,:,:,n-1) = ncon(tensors,legs);
            for x = n-1:-1:2
                tensors = {c(:,:,:,:,x),b(:,:,:,:,x),l1f(:,:,:,:,x)};
                legs = {[-1,3,1,2],[-2,4,2,1],[3,4,-3,-4]};
                l1f(:,:,:,:,x-1) = ncon(tensors,legs);
                tensors = {c(:,:,:,:,x),a(:,:,:,:,x),c(:,:,:,:,x),l2f(:,:,:,:,:,:,x)};
                legs = {[-1,4,1,2],[-2,5,2,3],[-3,6,3,1],[4,5,6,-4,-5,-6]};
                l2f(:,:,:,:,:,:,x-1) = ncon(tensors,legs);
            end
            tensors = {b(:,:,:,:,1),l1f(:,:,:,:,1)};
            legs = {[2,1,-4,-3],[-2,1,-1,2]};
            l1 = ncon(tensors,legs);
            l1 = reshape(l1,[bdl*bdl*d*d,1]);
            tensors = {a(:,:,:,:,1),eye(d),l2f(:,:,:,:,:,:,1)};
            legs = {[2,1,-4,-7],[-8,-3],[-2,1,-6,-1,2,-5]};
            l2 = ncon(tensors,legs);
            l2 = reshape(l2,[bdl*bdl*d*d,bdl*bdl*d*d]);
            dl2 = l2+l2.';
            dl1 = 2*l1;
            dl2pinv = pinv2(dl2,tol1);
            dl2pinv = (dl2pinv+dl2pinv.')/2;
            cv = dl2pinv*dl1;
            c(:,:,:,:,1) = reshape(cv,[bdl,bdl,d,d]);
            if lherm == 1
                c(:,:,:,:,1) = (c(:,:,:,:,1)+conj(permute(c(:,:,:,:,1),[1,2,4,3])))/2;
                cv = reshape(c(:,:,:,:,1),[bdl*bdl*d*d,1]);
            end
            index = index+1;
            fom1(index) = real(2*cv.'*l1-cv.'*l2*cv);
            tensors = {c(:,:,:,:,1),b(:,:,:,:,1)};
            legs = {[-1,-3,1,2],[-2,-4,2,1]};
            l1c = ncon(tensors,legs);
            tensors = {c(:,:,:,:,1),a(:,:,:,:,1),c(:,:,:,:,1)};
            legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
            l2c = ncon(tensors,legs);
            for x = 2:n-1
                tensors = {l1c,b(:,:,:,:,x),l1f(:,:,:,:,x)};
                legs = {[3,4,-1,1],[1,2,-4,-3],[-2,2,3,4]};
                l1 = ncon(tensors,legs);
                l1 = reshape(l1,[bdl*bdl*d*d,1]);
                tensors = {l2c,a(:,:,:,:,x),eye(d),l2f(:,:,:,:,:,:,x)};
                legs = {[3,4,5,-1,1,-5],[1,2,-4,-7],[-8,-3],[-2,2,-6,3,4,5]};
                l2 = ncon(tensors,legs);
                l2 = reshape(l2,[bdl*bdl*d*d,bdl*bdl*d*d]);
                dl2 = l2+l2.';
                dl1 = 2*l1;
                dl2pinv = pinv2(dl2,tol1);
                dl2pinv = (dl2pinv+dl2pinv.')/2;
                cv = dl2pinv*dl1;
                c(:,:,:,:,x) = reshape(cv,[bdl,bdl,d,d]);
                if lherm == 1
                    c(:,:,:,:,x) = (c(:,:,:,:,x)+conj(permute(c(:,:,:,:,x),[1,2,4,3])))/2;
                    cv = reshape(c(:,:,:,:,x),[bdl*bdl*d*d,1]);
                end
                index = index+1;
                fom1(index) = real(2*cv.'*l1-cv.'*l2*cv);
                tensors = {l1c,c(:,:,:,:,x),b(:,:,:,:,x)};
                legs = {[-1,-2,3,4],[3,-3,1,2],[4,-4,2,1]};
                l1c = ncon(tensors,legs);
                tensors = {l2c,c(:,:,:,:,x),a(:,:,:,:,x),c(:,:,:,:,x)};
                legs = {[-1,-2,-3,4,5,6],[4,-4,1,2],[5,-5,2,3],[6,-6,3,1]};
                l2c = ncon(tensors,legs);
            end
            tensors = {l1c,b(:,:,:,:,n)};
            legs = {[-2,2,-1,1],[1,2,-4,-3]};
            l1 = ncon(tensors,legs);
            l1 = reshape(l1,[bdl*bdl*d*d,1]);
            tensors = {l2c,a(:,:,:,:,n),eye(d)};
            legs = {[-2,2,-6,-1,1,-5],[1,2,-4,-7],[-8,-3]};
            l2 = ncon(tensors,legs);
            l2 = reshape(l2,[bdl*bdl*d*d,bdl*bdl*d*d]);
            dl2 = l2+l2.';
            dl1 = 2*l1;
            dl2pinv = pinv2(dl2,tol1);
            dl2pinv = (dl2pinv+dl2pinv.')/2;
            cv = dl2pinv*dl1;
            c(:,:,:,:,n) = reshape(cv,[bdl,bdl,d,d]);
            if lherm == 1
                c(:,:,:,:,n) = (c(:,:,:,:,n)+conj(permute(c(:,:,:,:,n),[1,2,4,3])))/2;
                cv = reshape(c(:,:,:,:,n),[bdl*bdl*d*d,1]);
            end
            index = index+1;
            fom1(index) = real(2*cv.'*l1-cv.'*l2*cv);
            if iter1 >= 2 && all(fom1((iter1-2)*n+1:iter1*n) > 0) && std(fom1((iter1-2)*n+1:iter1*n),1)/mean(fom1((iter1-2)*n+1:iter1*n)) <= relunc1
                break
            end
        end
        fom1 = fom1(1:index);
        %plot(fom1)
        f = fom1(index);
%  ||2*r'-rL-Lr||^2/||2*r'||^2 minimization
    elseif figureofmerit == 2
        tensors = {b(:,:,:,:,n),b(:,:,:,:,n)};
        legs = {[-1,-3,1,2],[-2,-4,2,1]};
        l0 = ncon(tensors,legs);
        for x = n-1:-1:2
            tensors = {b(:,:,:,:,x),b(:,:,:,:,x),l0};
            legs = {[-1,3,1,2],[-2,4,2,1],[3,4,-3,-4]};
            l0 = ncon(tensors,legs);
        end
        tensors = {b(:,:,:,:,1),b(:,:,:,:,1),l0};
        legs = {[5,3,1,2],[6,4,2,1],[3,4,5,6]};
        l0 = ncon(tensors,legs);
        l0 = abs(l0);
        if lherm == 1
            l1_1f = zeros([bdl,2*bdpsi^2*bdnoise,bdpsi^2*bdnoise,bdl,2*bdpsi^2*bdnoise,bdpsi^2*bdnoise,n-1]);
            l1_2f = zeros([bdl,bdpsi^2*bdnoise,2*bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise,2*bdpsi^2*bdnoise,n-1]);
            l2_1f = zeros([bdl,bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise,n-1]);
            l2_2f = zeros([bdl,bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise,bdl,bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise,n-1]);
            fom2 = zeros([n*100,1]);
            index = 0;
            iter2 = 0;
            while 1
                iter2 = iter2+1;
                tensors = {c(:,:,:,:,n),b(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                l1_1f(:,:,:,:,:,:,n-1) = ncon(tensors,legs);
                tensors = {c(:,:,:,:,n),a(:,:,:,:,n),b(:,:,:,:,n)};
                legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                l1_2f(:,:,:,:,:,:,n-1) = ncon(tensors,legs);
                tensors = {c(:,:,:,:,n),a(:,:,:,:,n),c(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                l2_1f(:,:,:,:,:,:,:,:,n-1) = ncon(tensors,legs);
                tensors = {c(:,:,:,:,n),c(:,:,:,:,n),a(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                l2_2f(:,:,:,:,:,:,:,:,n-1) = ncon(tensors,legs);
                for x = n-1:-1:2
                    tensors = {c(:,:,:,:,x),b(:,:,:,:,x),a(:,:,:,:,x),l1_1f(:,:,:,:,:,:,x)};
                    legs = {[-1,4,1,2],[-2,5,2,3],[-3,6,3,1],[4,5,6,-4,-5,-6]};
                    l1_1f(:,:,:,:,:,:,x-1) = ncon(tensors,legs);
                    tensors = {c(:,:,:,:,x),a(:,:,:,:,x),b(:,:,:,:,x),l1_2f(:,:,:,:,:,:,x)};
                    legs = {[-1,4,1,2],[-2,5,2,3],[-3,6,3,1],[4,5,6,-4,-5,-6]};
                    l1_2f(:,:,:,:,:,:,x-1) = ncon(tensors,legs);
                    tensors = {c(:,:,:,:,x),a(:,:,:,:,x),c(:,:,:,:,x),a(:,:,:,:,x),l2_1f(:,:,:,:,:,:,:,:,x)};
                    legs = {[-1,5,1,2],[-2,6,2,3],[-3,7,3,4],[-4,8,4,1],[5,6,7,8,-5,-6,-7,-8]};
                    l2_1f(:,:,:,:,:,:,:,:,x-1) = ncon(tensors,legs);
                    tensors = {c(:,:,:,:,x),c(:,:,:,:,x),a(:,:,:,:,x),a(:,:,:,:,x),l2_2f(:,:,:,:,:,:,:,:,x)};
                    legs = {[-1,5,1,2],[-2,6,2,3],[-3,7,3,4],[-4,8,4,1],[5,6,7,8,-5,-6,-7,-8]};
                    l2_2f(:,:,:,:,:,:,:,:,x-1) = ncon(tensors,legs);
                end
                tensors = {b(:,:,:,:,1),a(:,:,:,:,1),l1_1f(:,:,:,:,:,:,1)};
                legs = {[4,2,-4,1],[5,3,1,-3],[-2,2,3,-1,4,5]};
                l1_1 = ncon(tensors,legs);
                l1_1 = reshape(l1_1,[bdl*bdl*d*d,1]);
                tensors = {a(:,:,:,:,1),b(:,:,:,:,1),l1_2f(:,:,:,:,:,:,1)};
                legs = {[4,2,-4,1],[5,3,1,-3],[-2,2,3,-1,4,5]};
                l1_2 = ncon(tensors,legs);
                l1_2 = reshape(l1_2,[bdl*bdl*d*d,1]);
                tensors = {a(:,:,:,:,1),a(:,:,:,:,1),l2_1f(:,:,:,:,:,:,:,:,1)};
                legs = {[2,1,-4,-7],[4,3,-8,-3],[-2,1,-6,3,-1,2,-5,4]};
                l2_1 = ncon(tensors,legs);
                l2_1 = reshape(l2_1,[bdl*bdl*d*d,bdl*bdl*d*d]);
                tensors = {eye(d),a(:,:,:,:,1),a(:,:,:,:,1),l2_2f(:,:,:,:,:,:,:,:,1)};
                legs = {[-4,-7],[4,2,-8,1],[5,3,1,-3],[-2,-6,2,3,-1,-5,4,5]};
                l2_2 = ncon(tensors,legs);
                l2_2 = reshape(l2_2,[bdl*bdl*d*d,bdl*bdl*d*d]);
                dl2 = l2_1+l2_1.'+l2_2+l2_2.';
                dl1 = 2*(l1_1+l1_2);
                dl2pinv = pinv2(dl2,tol2);
                dl2pinv = (dl2pinv+dl2pinv.')/2;
                cv = dl2pinv*dl1;
                c(:,:,:,:,1) = reshape(cv,[bdl,bdl,d,d]);
                c(:,:,:,:,1) = (c(:,:,:,:,1)+conj(permute(c(:,:,:,:,1),[1,2,4,3])))/2;
                cv = reshape(c(:,:,:,:,1),[bdl*bdl*d*d,1]);
                index = index+1;
                fom2(index) = abs(1+(-4*(cv.'*l1_1+cv.'*l1_2)+2*(cv.'*l2_1*cv+cv.'*l2_2*cv))/(4*l0));
                tensors = {c(:,:,:,:,1),b(:,:,:,:,1),a(:,:,:,:,1)};
                legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                l1_1c = ncon(tensors,legs);
                tensors = {c(:,:,:,:,1),a(:,:,:,:,1),b(:,:,:,:,1)};
                legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                l1_2c = ncon(tensors,legs);
                tensors = {c(:,:,:,:,1),a(:,:,:,:,1),c(:,:,:,:,1),a(:,:,:,:,1)};
                legs = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                l2_1c = ncon(tensors,legs);
                tensors = {c(:,:,:,:,1),c(:,:,:,:,1),a(:,:,:,:,1),a(:,:,:,:,1)};
                legs = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                l2_2c = ncon(tensors,legs);
                for x = 2:n-1
                    tensors = {l1_1c,b(:,:,:,:,x),a(:,:,:,:,x),l1_1f(:,:,:,:,:,:,x)};
                    legs = {[6,7,8,-1,2,3],[2,4,-4,1],[3,5,1,-3],[-2,4,5,6,7,8]};
                    l1_1 = ncon(tensors,legs);
                    l1_1 = reshape(l1_1,[bdl*bdl*d*d,1]);
                    tensors = {l1_2c,a(:,:,:,:,x),b(:,:,:,:,x),l1_2f(:,:,:,:,:,:,x)};
                    legs = {[6,7,8,-1,2,3],[2,4,-4,1],[3,5,1,-3],[-2,4,5,6,7,8]};
                    l1_2 = ncon(tensors,legs);
                    l1_2 = reshape(l1_2,[bdl*bdl*d*d,1]);
                    tensors = {l2_1c,a(:,:,:,:,x),a(:,:,:,:,x),l2_1f(:,:,:,:,:,:,:,:,x)};
                    legs = {[5,6,7,8,-1,1,-5,2],[1,3,-4,-7],[2,4,-8,-3],[-2,3,-6,4,5,6,7,8]};
                    l2_1 = ncon(tensors,legs);
                    l2_1 = reshape(l2_1,[bdl*bdl*d*d,bdl*bdl*d*d]);
                    tensors = {l2_2c,eye(d),a(:,:,:,:,x),a(:,:,:,:,x),l2_2f(:,:,:,:,:,:,:,:,x)};
                    legs = {[6,7,8,9,-1,-5,2,3],[-4,-7],[2,4,-8,1],[3,5,1,-3],[-2,-6,4,5,6,7,8,9]};
                    l2_2 = ncon(tensors,legs);
                    l2_2 = reshape(l2_2,[bdl*bdl*d*d,bdl*bdl*d*d]);
                    dl2 = l2_1+l2_1.'+l2_2+l2_2.';
                    dl1 = 2*(l1_1+l1_2);
                    dl2pinv = pinv2(dl2,tol2);
                    dl2pinv = (dl2pinv+dl2pinv.')/2;
                    cv = dl2pinv*dl1;
                    c(:,:,:,:,x) = reshape(cv,[bdl,bdl,d,d]);
                    c(:,:,:,:,x) = (c(:,:,:,:,x)+conj(permute(c(:,:,:,:,x),[1,2,4,3])))/2;
                    cv = reshape(c(:,:,:,:,x),[bdl*bdl*d*d,1]);
                    index = index+1;
                    fom2(index) = abs(1+(-4*(cv.'*l1_1+cv.'*l1_2)+2*(cv.'*l2_1*cv+cv.'*l2_2*cv))/(4*l0));
                    tensors = {l1_1c,c(:,:,:,:,x),b(:,:,:,:,x),a(:,:,:,:,x)};
                    legs = {[-1,-2,-3,4,5,6],[4,-4,1,2],[5,-5,2,3],[6,-6,3,1]};
                    l1_1c = ncon(tensors,legs);
                    tensors = {l1_2c,c(:,:,:,:,x),a(:,:,:,:,x),b(:,:,:,:,x)};
                    legs = {[-1,-2,-3,4,5,6],[4,-4,1,2],[5,-5,2,3],[6,-6,3,1]};
                    l1_2c = ncon(tensors,legs);
                    tensors = {l2_1c,c(:,:,:,:,x),a(:,:,:,:,x),c(:,:,:,:,x),a(:,:,:,:,x)};
                    legs = {[-1,-2,-3,-4,5,6,7,8],[5,-5,1,2],[6,-6,2,3],[7,-7,3,4],[8,-8,4,1]};
                    l2_1c = ncon(tensors,legs);
                    tensors = {l2_2c,c(:,:,:,:,x),c(:,:,:,:,x),a(:,:,:,:,x),a(:,:,:,:,x)};
                    legs = {[-1,-2,-3,-4,5,6,7,8],[5,-5,1,2],[6,-6,2,3],[7,-7,3,4],[8,-8,4,1]};
                    l2_2c = ncon(tensors,legs);
                end
                tensors = {l1_1c,b(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-2,4,5,-1,2,3],[2,4,-4,1],[3,5,1,-3]};
                l1_1 = ncon(tensors,legs);
                l1_1 = reshape(l1_1,[bdl*bdl*d*d,1]);
                tensors = {l1_2c,a(:,:,:,:,n),b(:,:,:,:,n)};
                legs = {[-2,4,5,-1,2,3],[2,4,-4,1],[3,5,1,-3]};
                l1_2 = ncon(tensors,legs);
                l1_2 = reshape(l1_2,[bdl*bdl*d*d,1]);
                tensors = {l2_1c,a(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-2,2,-6,4,-1,1,-5,3],[1,2,-4,-7],[3,4,-8,-3]};
                l2_1 = ncon(tensors,legs);
                l2_1 = reshape(l2_1,[bdl*bdl*d*d,bdl*bdl*d*d]);
                tensors = {l2_2c,eye(d),a(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-2,-6,4,5,-1,-5,2,3],[-4,-7],[2,4,-8,1],[3,5,1,-3]};
                l2_2 = ncon(tensors,legs);
                l2_2 = reshape(l2_2,[bdl*bdl*d*d,bdl*bdl*d*d]);
                dl2 = l2_1+l2_1.'+l2_2+l2_2.';
                dl1 = 2*(l1_1+l1_2);
                dl2pinv = pinv2(dl2,tol2);
                dl2pinv = (dl2pinv+dl2pinv.')/2;
                cv = dl2pinv*dl1;
                c(:,:,:,:,n) = reshape(cv,[bdl,bdl,d,d]);
                c(:,:,:,:,n) = (c(:,:,:,:,n)+conj(permute(c(:,:,:,:,n),[1,2,4,3])))/2;
                cv = reshape(c(:,:,:,:,n),[bdl*bdl*d*d,1]);
                index = index+1;
                fom2(index) = abs(1+(-4*(cv.'*l1_1+cv.'*l1_2)+2*(cv.'*l2_1*cv+cv.'*l2_2*cv))/(4*l0));
                if fom2(index) < 10^-10 || iter2 >= 20 || (iter2 >= 2 && std(fom2((iter2-2)*n+1:iter2*n),1)/mean(fom2((iter2-2)*n+1:iter2*n)) <= relunc2)
                    break
                end
            end
        else
            l1_1f = zeros([bdl,2*bdpsi^2*bdnoise,bdpsi^2*bdnoise,bdl,2*bdpsi^2*bdnoise,bdpsi^2*bdnoise,n-1]);
            l1_2f = zeros([bdl,bdpsi^2*bdnoise,2*bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise,2*bdpsi^2*bdnoise,n-1]);
            l1_3f = zeros([bdl,2*bdpsi^2*bdnoise,bdpsi^2*bdnoise,bdl,2*bdpsi^2*bdnoise,bdpsi^2*bdnoise,n-1]);
            l1_4f = zeros([bdl,bdpsi^2*bdnoise,2*bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise,2*bdpsi^2*bdnoise,n-1]);
            l2_1f = zeros([bdl,bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise,n-1]);
            l2_2f = zeros([bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise,bdl,bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise,bdl,n-1]);
            l2_3f = zeros([bdl,bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise,bdl,bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise,n-1]);
            fom2 = zeros([n*100,1]);
            index = 0;
            iter2 = 0;
            while 1
                iter2 = iter2+1;
                tensors = {conj(c(:,:,:,:,n)),b(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
                l1_1f(:,:,:,:,:,:,n-1) = ncon(tensors,legs);
                tensors = {conj(c(:,:,:,:,n)),a(:,:,:,:,n),b(:,:,:,:,n)};
                legs = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
                l1_2f(:,:,:,:,:,:,n-1) = ncon(tensors,legs);
                tensors = {c(:,:,:,:,n),b(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                l1_3f(:,:,:,:,:,:,n-1) = ncon(tensors,legs);
                tensors = {c(:,:,:,:,n),a(:,:,:,:,n),b(:,:,:,:,n)};
                legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                l1_4f(:,:,:,:,:,:,n-1) = ncon(tensors,legs);
                tensors = {conj(c(:,:,:,:,n)),a(:,:,:,:,n),c(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                l2_1f(:,:,:,:,:,:,:,:,n-1) = ncon(tensors,legs);
                tensors = {conj(c(:,:,:,:,n)),a(:,:,:,:,n),a(:,:,:,:,n),c(:,:,:,:,n)};
                legs = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                l2_2f(:,:,:,:,:,:,:,:,n-1) = ncon(tensors,legs);
                tensors = {conj(c(:,:,:,:,n)),c(:,:,:,:,n),a(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                l2_3f(:,:,:,:,:,:,:,:,n-1) = ncon(tensors,legs);
                for x = n-1:-1:2
                    tensors = {conj(c(:,:,:,:,x)),b(:,:,:,:,x),a(:,:,:,:,x),l1_1f(:,:,:,:,:,:,x)};
                    legs = {[-1,4,2,1],[-2,5,2,3],[-3,6,3,1],[4,5,6,-4,-5,-6]};
                    l1_1f(:,:,:,:,:,:,x-1) = ncon(tensors,legs);
                    tensors = {conj(c(:,:,:,:,x)),a(:,:,:,:,x),b(:,:,:,:,x),l1_2f(:,:,:,:,:,:,x)};
                    legs = {[-1,4,2,1],[-2,5,2,3],[-3,6,3,1],[4,5,6,-4,-5,-6]};
                    l1_2f(:,:,:,:,:,:,x-1) = ncon(tensors,legs);
                    tensors = {c(:,:,:,:,x),b(:,:,:,:,x),a(:,:,:,:,x),l1_3f(:,:,:,:,:,:,x)};
                    legs = {[-1,4,1,2],[-2,5,2,3],[-3,6,3,1],[4,5,6,-4,-5,-6]};
                    l1_3f(:,:,:,:,:,:,x-1) = ncon(tensors,legs);
                    tensors = {c(:,:,:,:,x),a(:,:,:,:,x),b(:,:,:,:,x),l1_4f(:,:,:,:,:,:,x)};
                    legs = {[-1,4,1,2],[-2,5,2,3],[-3,6,3,1],[4,5,6,-4,-5,-6]};
                    l1_4f(:,:,:,:,:,:,x-1) = ncon(tensors,legs);
                    tensors = {conj(c(:,:,:,:,x)),a(:,:,:,:,x),c(:,:,:,:,x),a(:,:,:,:,x),l2_1f(:,:,:,:,:,:,:,:,x)};
                    legs = {[-1,5,2,1],[-2,6,2,3],[-3,7,3,4],[-4,8,4,1],[5,6,7,8,-5,-6,-7,-8]};
                    l2_1f(:,:,:,:,:,:,:,:,x-1) = ncon(tensors,legs);
                    tensors = {conj(c(:,:,:,:,x)),a(:,:,:,:,x),a(:,:,:,:,x),c(:,:,:,:,x),l2_2f(:,:,:,:,:,:,:,:,x)};
                    legs = {[-1,5,2,1],[-2,6,2,3],[-3,7,3,4],[-4,8,4,1],[5,6,7,8,-5,-6,-7,-8]};
                    l2_2f(:,:,:,:,:,:,:,:,x-1) = ncon(tensors,legs);
                    tensors = {conj(c(:,:,:,:,x)),c(:,:,:,:,x),a(:,:,:,:,x),a(:,:,:,:,x),l2_3f(:,:,:,:,:,:,:,:,x)};
                    legs = {[-1,5,2,1],[-2,6,2,3],[-3,7,3,4],[-4,8,4,1],[5,6,7,8,-5,-6,-7,-8]};
                    l2_3f(:,:,:,:,:,:,:,:,x-1) = ncon(tensors,legs);
                end
                tensors = {b(:,:,:,:,1),a(:,:,:,:,1),l1_1f(:,:,:,:,:,:,1)};
                legs = {[4,2,-3,1],[5,3,1,-4],[-2,2,3,-1,4,5]};
                l1_1 = ncon(tensors,legs);
                l1_1 = reshape(l1_1,[bdl*bdl*d*d,1]);
                tensors = {a(:,:,:,:,1),b(:,:,:,:,1),l1_2f(:,:,:,:,:,:,1)};
                legs = {[4,2,-3,1],[5,3,1,-4],[-2,2,3,-1,4,5]};
                l1_2 = ncon(tensors,legs);
                l1_2 = reshape(l1_2,[bdl*bdl*d*d,1]);
                tensors = {b(:,:,:,:,1),a(:,:,:,:,1),l1_3f(:,:,:,:,:,:,1)};
                legs = {[4,2,-4,1],[5,3,1,-3],[-2,2,3,-1,4,5]};
                l1_3 = ncon(tensors,legs);
                l1_3 = reshape(l1_3,[bdl*bdl*d*d,1]);
                tensors = {a(:,:,:,:,1),b(:,:,:,:,1),l1_4f(:,:,:,:,:,:,1)};
                legs = {[4,2,-4,1],[5,3,1,-3],[-2,2,3,-1,4,5]};
                l1_4 = ncon(tensors,legs);
                l1_4 = reshape(l1_4,[bdl*bdl*d*d,1]);
                tensors = {a(:,:,:,:,1),a(:,:,:,:,1),l2_1f(:,:,:,:,:,:,:,:,1)};
                legs = {[2,1,-3,-7],[4,3,-8,-4],[-2,1,-6,3,-1,2,-5,4]};
                l2_1 = ncon(tensors,legs);
                l2_1 = reshape(l2_1,[bdl*bdl*d*d,bdl*bdl*d*d]);
                tensors = {a(:,:,:,:,1),a(:,:,:,:,1),eye(d),l2_2f(:,:,:,:,:,:,:,:,1)};
                legs = {[4,2,-3,1],[5,3,1,-7],[-8,-4],[-2,2,3,-6,-1,4,5,-5]};
                l2_2 = ncon(tensors,legs);
                l2_2 = reshape(l2_2,[bdl*bdl*d*d,bdl*bdl*d*d]);
                tensors = {eye(d),a(:,:,:,:,1),a(:,:,:,:,1),l2_3f(:,:,:,:,:,:,:,:,1)};
                legs = {[-3,-7],[4,2,-8,1],[5,3,1,-4],[-2,-6,2,3,-1,-5,4,5]};
                l2_3 = ncon(tensors,legs);
                l2_3 = reshape(l2_3,[bdl*bdl*d*d,bdl*bdl*d*d]);
                dl2 = 2*l2_1+l2_2+l2_3;
                dl2 = (dl2+dl2')/2;
                dl1 = 2*(l1_1+l1_2);
                dl2pinv = pinv2(dl2,tol2);
                dl2pinv = (dl2pinv+dl2pinv')/2;
                cv = dl2pinv*dl1;
                c(:,:,:,:,1) = reshape(cv,[bdl,bdl,d,d]);
                index = index+1;
                fom2(index) = abs(1+(-2*(cv'*l1_1+cv'*l1_2+cv.'*l1_3+cv.'*l1_4)+2*cv'*l2_1*cv+cv'*l2_2*cv+cv'*l2_3*cv)/(4*l0));
                tensors = {conj(c(:,:,:,:,1)),b(:,:,:,:,1),a(:,:,:,:,1)};
                legs = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
                l1_1c = ncon(tensors,legs);
                tensors = {conj(c(:,:,:,:,1)),a(:,:,:,:,1),b(:,:,:,:,1)};
                legs = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
                l1_2c = ncon(tensors,legs);
                tensors = {c(:,:,:,:,1),b(:,:,:,:,1),a(:,:,:,:,1)};
                legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                l1_3c = ncon(tensors,legs);
                tensors = {c(:,:,:,:,1),a(:,:,:,:,1),b(:,:,:,:,1)};
                legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                l1_4c = ncon(tensors,legs);
                tensors = {conj(c(:,:,:,:,1)),a(:,:,:,:,1),c(:,:,:,:,1),a(:,:,:,:,1)};
                legs = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                l2_1c = ncon(tensors,legs);
                tensors = {conj(c(:,:,:,:,1)),a(:,:,:,:,1),a(:,:,:,:,1),c(:,:,:,:,1)};
                legs = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                l2_2c = ncon(tensors,legs);
                tensors = {conj(c(:,:,:,:,1)),c(:,:,:,:,1),a(:,:,:,:,1),a(:,:,:,:,1)};
                legs = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                l2_3c = ncon(tensors,legs);
                for x = 2:n-1
                    tensors = {l1_1c,b(:,:,:,:,x),a(:,:,:,:,x),l1_1f(:,:,:,:,:,:,x)};
                    legs = {[6,7,8,-1,2,3],[2,4,-3,1],[3,5,1,-4],[-2,4,5,6,7,8]};
                    l1_1 = ncon(tensors,legs);
                    l1_1 = reshape(l1_1,[bdl*bdl*d*d,1]);
                    tensors = {l1_2c,a(:,:,:,:,x),b(:,:,:,:,x),l1_2f(:,:,:,:,:,:,x)};
                    legs = {[6,7,8,-1,2,3],[2,4,-3,1],[3,5,1,-4],[-2,4,5,6,7,8]};
                    l1_2 = ncon(tensors,legs);
                    l1_2 = reshape(l1_2,[bdl*bdl*d*d,1]);
                    tensors = {l1_3c,b(:,:,:,:,x),a(:,:,:,:,x),l1_3f(:,:,:,:,:,:,x)};
                    legs = {[6,7,8,-1,2,3],[2,4,-4,1],[3,5,1,-3],[-2,4,5,6,7,8]};
                    l1_3 = ncon(tensors,legs);
                    l1_3 = reshape(l1_3,[bdl*bdl*d*d,1]);
                    tensors = {l1_4c,a(:,:,:,:,x),b(:,:,:,:,x),l1_4f(:,:,:,:,:,:,x)};
                    legs = {[6,7,8,-1,2,3],[2,4,-4,1],[3,5,1,-3],[-2,4,5,6,7,8]};
                    l1_4 = ncon(tensors,legs);
                    l1_4 = reshape(l1_4,[bdl*bdl*d*d,1]);
                    tensors = {l2_1c,a(:,:,:,:,x),a(:,:,:,:,x),l2_1f(:,:,:,:,:,:,:,:,x)};
                    legs = {[5,6,7,8,-1,1,-5,2],[1,3,-3,-7],[2,4,-8,-4],[-2,3,-6,4,5,6,7,8]};
                    l2_1 = ncon(tensors,legs);
                    l2_1 = reshape(l2_1,[bdl*bdl*d*d,bdl*bdl*d*d]);
                    tensors = {l2_2c,a(:,:,:,:,x),a(:,:,:,:,x),eye(d),l2_2f(:,:,:,:,:,:,:,:,x)};
                    legs = {[6,7,8,9,-1,2,3,-5],[2,4,-3,1],[3,5,1,-7],[-8,-4],[-2,4,5,-6,6,7,8,9]};
                    l2_2 = ncon(tensors,legs);
                    l2_2 = reshape(l2_2,[bdl*bdl*d*d,bdl*bdl*d*d]);
                    tensors = {l2_3c,eye(d),a(:,:,:,:,x),a(:,:,:,:,x),l2_3f(:,:,:,:,:,:,:,:,x)};
                    legs = {[6,7,8,9,-1,-5,2,3],[-3,-7],[2,4,-8,1],[3,5,1,-4],[-2,-6,4,5,6,7,8,9]};
                    l2_3 = ncon(tensors,legs);
                    l2_3 = reshape(l2_3,[bdl*bdl*d*d,bdl*bdl*d*d]);
                    dl2 = 2*l2_1+l2_2+l2_3;
                    dl2 = (dl2+dl2')/2;
                    dl1 = 2*(l1_1+l1_2);
                    dl2pinv = pinv2(dl2,tol2);
                    dl2pinv = (dl2pinv+dl2pinv')/2;
                    cv = dl2pinv*dl1;
                    c(:,:,:,:,x) = reshape(cv,[bdl,bdl,d,d]);
                    index = index+1;
                    fom2(index) = abs(1+(-2*(cv'*l1_1+cv'*l1_2+cv.'*l1_3+cv.'*l1_4)+2*cv'*l2_1*cv+cv'*l2_2*cv+cv'*l2_3*cv)/(4*l0));
                    tensors = {l1_1c,conj(c(:,:,:,:,x)),b(:,:,:,:,x),a(:,:,:,:,x)};
                    legs = {[-1,-2,-3,4,5,6],[4,-4,2,1],[5,-5,2,3],[6,-6,3,1]};
                    l1_1c = ncon(tensors,legs);
                    tensors = {l1_2c,conj(c(:,:,:,:,x)),a(:,:,:,:,x),b(:,:,:,:,x)};
                    legs = {[-1,-2,-3,4,5,6],[4,-4,2,1],[5,-5,2,3],[6,-6,3,1]};
                    l1_2c = ncon(tensors,legs);
                    tensors = {l1_3c,c(:,:,:,:,x),b(:,:,:,:,x),a(:,:,:,:,x)};
                    legs = {[-1,-2,-3,4,5,6],[4,-4,1,2],[5,-5,2,3],[6,-6,3,1]};
                    l1_3c = ncon(tensors,legs);
                    tensors = {l1_4c,c(:,:,:,:,x),a(:,:,:,:,x),b(:,:,:,:,x)};
                    legs = {[-1,-2,-3,4,5,6],[4,-4,1,2],[5,-5,2,3],[6,-6,3,1]};
                    l1_4c = ncon(tensors,legs);
                    tensors = {l2_1c,conj(c(:,:,:,:,x)),a(:,:,:,:,x),c(:,:,:,:,x),a(:,:,:,:,x)};
                    legs = {[-1,-2,-3,-4,5,6,7,8],[5,-5,2,1],[6,-6,2,3],[7,-7,3,4],[8,-8,4,1]};
                    l2_1c = ncon(tensors,legs);
                    tensors = {l2_2c,conj(c(:,:,:,:,x)),a(:,:,:,:,x),a(:,:,:,:,x),c(:,:,:,:,x)};
                    legs = {[-1,-2,-3,-4,5,6,7,8],[5,-5,2,1],[6,-6,2,3],[7,-7,3,4],[8,-8,4,1]};
                    l2_2c = ncon(tensors,legs);
                    tensors = {l2_3c,conj(c(:,:,:,:,x)),c(:,:,:,:,x),a(:,:,:,:,x),a(:,:,:,:,x)};
                    legs = {[-1,-2,-3,-4,5,6,7,8],[5,-5,2,1],[6,-6,2,3],[7,-7,3,4],[8,-8,4,1]};
                    l2_3c = ncon(tensors,legs);
                end
                tensors = {l1_1c,b(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-2,4,5,-1,2,3],[2,4,-3,1],[3,5,1,-4]};
                l1_1 = ncon(tensors,legs);
                l1_1 = reshape(l1_1,[bdl*bdl*d*d,1]);
                tensors = {l1_2c,a(:,:,:,:,n),b(:,:,:,:,n)};
                legs = {[-2,4,5,-1,2,3],[2,4,-3,1],[3,5,1,-4]};
                l1_2 = ncon(tensors,legs);
                l1_2 = reshape(l1_2,[bdl*bdl*d*d,1]);
                tensors = {l1_3c,b(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-2,4,5,-1,2,3],[2,4,-4,1],[3,5,1,-3]};
                l1_3 = ncon(tensors,legs);
                l1_3 = reshape(l1_3,[bdl*bdl*d*d,1]);
                tensors = {l1_4c,a(:,:,:,:,n),b(:,:,:,:,n)};
                legs = {[-2,4,5,-1,2,3],[2,4,-4,1],[3,5,1,-3]};
                l1_4 = ncon(tensors,legs);
                l1_4 = reshape(l1_4,[bdl*bdl*d*d,1]);
                tensors = {l2_1c,a(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-2,2,-6,4,-1,1,-5,3],[1,2,-3,-7],[3,4,-8,-4]};
                l2_1 = ncon(tensors,legs);
                l2_1 = reshape(l2_1,[bdl*bdl*d*d,bdl*bdl*d*d]);
                tensors = {l2_2c,a(:,:,:,:,n),a(:,:,:,:,n),eye(d)};
                legs = {[-2,4,5,-6,-1,2,3,-5],[2,4,-3,1],[3,5,1,-7],[-8,-4]};
                l2_2 = ncon(tensors,legs);
                l2_2 = reshape(l2_2,[bdl*bdl*d*d,bdl*bdl*d*d]);
                tensors = {l2_3c,eye(d),a(:,:,:,:,n),a(:,:,:,:,n)};
                legs = {[-2,-6,4,5,-1,-5,2,3],[-3,-7],[2,4,-8,1],[3,5,1,-4]};
                l2_3 = ncon(tensors,legs);
                l2_3 = reshape(l2_3,[bdl*bdl*d*d,bdl*bdl*d*d]);
                dl2 = 2*l2_1+l2_2+l2_3;
                dl2 = (dl2+dl2')/2;
                dl1 = 2*(l1_1+l1_2);
                dl2pinv = pinv2(dl2,tol2);
                dl2pinv = (dl2pinv+dl2pinv')/2;
                cv = dl2pinv*dl1;
                c(:,:,:,:,n) = reshape(cv,[bdl,bdl,d,d]);
                index = index+1;
                fom2(index) = abs(1+(-2*(cv'*l1_1+cv'*l1_2+cv.'*l1_3+cv.'*l1_4)+2*cv'*l2_1*cv+cv'*l2_2*cv+cv'*l2_3*cv)/(4*l0));
                if fom2(index) < 10^-10 || iter2 >= 20 || (iter2 >= 2 && std(fom2((iter2-2)*n+1:iter2*n),1)/mean(fom2((iter2-2)*n+1:iter2*n)) <= relunc2)
                    break
                end
            end
        end
        fom2 = fom2(1:index);
        %plot(-log10(fom2))
        tensors = {b(:,:,:,:,n),c(:,:,:,:,n)};
        legs = {[-1,-3,1,2],[-2,-4,2,1]};
        f = ncon(tensors,legs);
        for x = n-1:-1:2
            tensors = {b(:,:,:,:,x),c(:,:,:,:,x),f};
            legs = {[-1,3,1,2],[-2,4,2,1],[3,4,-3,-4]};
            f = ncon(tensors,legs);
        end
        tensors = {b(:,:,:,:,1),c(:,:,:,:,1),f};
        legs = {[5,3,1,2],[6,4,2,1],[3,4,5,6]};
        f = ncon(tensors,legs);
        f = real(f);
    end
    f_iter(2*iter-1) = f;
	if iter >= 3 && std(f_iter(2*iter-4:2*iter-1),1)/mean(f_iter(2*iter-4:2*iter-1)) <= relunc
        break
	end
    if figureofmerit == 2 && iter >= 20
        break
    end
%  MPO (dual)
    if lherm == 1
        bdherm = 1;
        ch = c;
    else
        bdherm = 2;
        ch = zeros([bdherm*bdl,bdherm*bdl,d,d,n]);
        for nx = 1:d
            for nxp = 1:d
                ch(:,:,nx,nxp,:) = [c(:,:,nx,nxp,:),zeros([bdl,bdl,1,1,n]);zeros([bdl,bdl,1,1,n]),conj(c(:,:,nxp,nx,:))];
            end
        end
        ch = ch/2^(1/n);
    end
    c2 = zeros([(bdherm*bdl)^2,(bdherm*bdl)^2,d,d,n]);
    for nx = 1:d
        for nxp = 1:d
            for nxpp = 1:d
                for x = 1:n
                    c2(:,:,nx,nxp,x) = c2(:,:,nx,nxp,x)+kron(ch(:,:,nx,nxpp,x),ch(:,:,nxpp,nxp,x));
                end
            end
        end
    end
    cp = zeros([2*bdherm*bdl,2*bdherm*bdl,d,d,n]);
    for nx = 1:d
        for nxp = 1:d
            cp(:,:,nx,nxp,1) = kron([1i*(nx-nxp),1;0,0],ch(:,:,nx,nxp,1));
            for x = 2:n-1
                cp(:,:,nx,nxp,x) = kron([1,0;1i*(nx-nxp),1],ch(:,:,nx,nxp,x));
            end
            cp(:,:,nx,nxp,n) = kron([1,0;1i*(nx-nxp),0],ch(:,:,nx,nxp,n));
        end
    end
    c2d = zeros([(bdherm*bdl)^2*bdnoise,(bdherm*bdl)^2*bdnoise,d,d,n]);
    cpd = zeros([2*bdherm*bdl*bdnoise,2*bdherm*bdl*bdnoise,d,d,n]);
    for nx = 1:d
        for nxp = 1:d
            for x = 1:n
                c2d(:,:,nx,nxp,x) = kron(an(:,:,nx,nxp),c2(:,:,nx,nxp,x));
                cpd(:,:,nx,nxp,x) = kron(an(:,:,nx,nxp),cp(:,:,nx,nxp,x));
            end
        end
    end
%  Tr[r0*X] maximization
    l2df = zeros([bdpsi,(bdherm*bdl)^2*bdnoise,bdpsi,bdpsi,(bdherm*bdl)^2*bdnoise,bdpsi,n-1]);
    lpdf = zeros([bdpsi,2*bdherm*bdl*bdnoise,bdpsi,bdpsi,2*bdherm*bdl*bdnoise,bdpsi,n-1]);
    psinormf = zeros(bdpsi,bdpsi,bdpsi,bdpsi,n-1);
    fomd = zeros([n*100,1]);
    indexd = 0;
    iterd = 0;
    while 1
        iterd = iterd+1;
        tensors = {conj(a0(:,:,:,n)),c2d(:,:,:,:,n),a0(:,:,:,n)};
        legs = {[-1,-4,1],[-2,-5,1,2],[-3,-6,2]};
        l2df(:,:,:,:,:,:,n-1) = ncon(tensors,legs);
        tensors = {conj(a0(:,:,:,n)),cpd(:,:,:,:,n),a0(:,:,:,n)};
        legs = {[-1,-4,1],[-2,-5,1,2],[-3,-6,2]};
        lpdf(:,:,:,:,:,:,n-1) = ncon(tensors,legs);
        tensors = {conj(a0(:,:,:,n)),a0(:,:,:,n)};
        legs = {[-1,-3,1],[-2,-4,1]};
        psinormf(:,:,:,:,n-1) = ncon(tensors,legs);
        for x = n-1:-1:2
            tensors = {conj(a0(:,:,:,x)),c2d(:,:,:,:,x),a0(:,:,:,x),l2df(:,:,:,:,:,:,x)};
            legs = {[-1,4,1],[-2,5,1,2],[-3,6,2],[4,5,6,-4,-5,-6]};
            l2df(:,:,:,:,:,:,x-1) = ncon(tensors,legs);
            tensors = {conj(a0(:,:,:,x)),cpd(:,:,:,:,x),a0(:,:,:,x),lpdf(:,:,:,:,:,:,x)};
            legs = {[-1,4,1],[-2,5,1,2],[-3,6,2],[4,5,6,-4,-5,-6]};
            lpdf(:,:,:,:,:,:,x-1) = ncon(tensors,legs);
            tensors = {conj(a0(:,:,:,x)),a0(:,:,:,x),psinormf(:,:,:,:,x)};
            legs = {[-1,2,1],[-2,3,1],[2,3,-3,-4]};
            psinormf(:,:,:,:,x-1) = ncon(tensors,legs);
        end
        tensors = {c2d(:,:,:,:,1),l2df(:,:,:,:,:,:,1)};
        legs = {[2,1,-3,-6],[-2,1,-5,-1,2,-4]};
        l2d = ncon(tensors,legs);
        l2d = reshape(l2d,[bdpsi*bdpsi*d,bdpsi*bdpsi*d]);
        tensors = {cpd(:,:,:,:,1),lpdf(:,:,:,:,:,:,1)};
        legs = {[2,1,-3,-6],[-2,1,-5,-1,2,-4]};
        lpd = ncon(tensors,legs);
        lpd = reshape(lpd,[bdpsi*bdpsi*d,bdpsi*bdpsi*d]);
        tensors = {psinormf(:,:,:,:,1)};
        legs = {[-2,-4,-1,-3]};
        psinorm = ncon(tensors,legs);
        psinorm = reshape(psinorm,[bdpsi*bdpsi,bdpsi*bdpsi]);
        psinorm = (psinorm+psinorm')/2;
        psinormpinv = pinv2(psinorm,told);
        psinormpinv = (psinormpinv+psinormpinv')/2;
        psinormpinv = kron(eye(d),psinormpinv);
        eiginput = 2*lpd-l2d;
        eiginput = (eiginput+eiginput')/2;
        eiginput = psinormpinv*eiginput;
        [a0v,fomdval] = eigs(eiginput,1,'lr');
        a0v = a0v/abs(a0v'*kron(eye(d),psinorm)*a0v)^(1/2);
        a0(:,:,:,1) = reshape(a0v,[bdpsi,bdpsi,d]);
        indexd = indexd+1;
        fomd(indexd) = real(fomdval);
        tensors = {conj(a0(:,:,:,1)),c2d(:,:,:,:,1),a0(:,:,:,1)};
        legs = {[-1,-4,1],[-2,-5,1,2],[-3,-6,2]};
        l2dc = ncon(tensors,legs);
        tensors = {conj(a0(:,:,:,1)),cpd(:,:,:,:,1),a0(:,:,:,1)};
        legs = {[-1,-4,1],[-2,-5,1,2],[-3,-6,2]};
        lpdc = ncon(tensors,legs);
        tensors = {conj(a0(:,:,:,1)),a0(:,:,:,1)};
        legs = {[-1,-3,1],[-2,-4,1]};
        psinormc = ncon(tensors,legs);
        for x = 2:n-1
            tensors = {l2dc,c2d(:,:,:,:,x),l2df(:,:,:,:,:,:,x)};
            legs = {[3,4,5,-1,1,-4],[1,2,-3,-6],[-2,2,-5,3,4,5]};
            l2d = ncon(tensors,legs);
            l2d = reshape(l2d,[bdpsi*bdpsi*d,bdpsi*bdpsi*d]);
            tensors = {lpdc,cpd(:,:,:,:,x),lpdf(:,:,:,:,:,:,x)};
            legs = {[3,4,5,-1,1,-4],[1,2,-3,-6],[-2,2,-5,3,4,5]};
            lpd = ncon(tensors,legs);
            lpd = reshape(lpd,[bdpsi*bdpsi*d,bdpsi*bdpsi*d]);
            tensors = {psinormc,psinormf(:,:,:,:,x)};
            legs = {[1,2,-1,-3],[-2,-4,1,2]};
            psinorm = ncon(tensors,legs);
            psinorm = reshape(psinorm,[bdpsi*bdpsi,bdpsi*bdpsi]);
            psinorm = (psinorm+psinorm')/2;
            psinormpinv = pinv2(psinorm,told);
            psinormpinv = (psinormpinv+psinormpinv')/2;
            psinormpinv = kron(eye(d),psinormpinv);
            eiginput = 2*lpd-l2d;
            eiginput = (eiginput+eiginput')/2;
            eiginput = psinormpinv*eiginput;
            [a0v,fomdval] = eigs(eiginput,1,'lr');
            a0v = a0v/abs(a0v'*kron(eye(d),psinorm)*a0v)^(1/2);
            a0(:,:,:,x) = reshape(a0v,[bdpsi,bdpsi,d]);
            indexd = indexd+1;
            fomd(indexd) = real(fomdval);
            tensors = {l2dc,conj(a0(:,:,:,x)),c2d(:,:,:,:,x),a0(:,:,:,x)};
            legs = {[-1,-2,-3,3,4,5],[3,-4,1],[4,-5,1,2],[5,-6,2]};
            l2dc = ncon(tensors,legs);
            tensors = {lpdc,conj(a0(:,:,:,x)),cpd(:,:,:,:,x),a0(:,:,:,x)};
            legs = {[-1,-2,-3,3,4,5],[3,-4,1],[4,-5,1,2],[5,-6,2]};
            lpdc = ncon(tensors,legs);
            tensors = {psinormc,conj(a0(:,:,:,x)),a0(:,:,:,x)};
            legs = {[-1,-2,2,3],[2,-3,1],[3,-4,1]};
            psinormc = ncon(tensors,legs);
        end
        tensors = {l2dc,c2d(:,:,:,:,n)};
        legs = {[-2,2,-5,-1,1,-4],[1,2,-3,-6]};
        l2d = ncon(tensors,legs);
        l2d = reshape(l2d,[bdpsi*bdpsi*d,bdpsi*bdpsi*d]);
        tensors = {lpdc,cpd(:,:,:,:,n)};
        legs = {[-2,2,-5,-1,1,-4],[1,2,-3,-6]};
        lpd = ncon(tensors,legs);
        lpd = reshape(lpd,[bdpsi*bdpsi*d,bdpsi*bdpsi*d]);
        tensors = {psinormc};
        legs = {[-2,-4,-1,-3]};
        psinorm = ncon(tensors,legs);
        psinorm = reshape(psinorm,[bdpsi*bdpsi,bdpsi*bdpsi]);
        psinorm = (psinorm+psinorm')/2;
        psinormpinv = pinv2(psinorm,told);
        psinormpinv = (psinormpinv+psinormpinv')/2;
        psinormpinv = kron(eye(d),psinormpinv);
        eiginput = 2*lpd-l2d;
        eiginput = (eiginput+eiginput')/2;
        eiginput = psinormpinv*eiginput;
        [a0v,fomdval] = eigs(eiginput,1,'lr');
        a0v = a0v/abs(a0v'*kron(eye(d),psinorm)*a0v)^(1/2);
        a0(:,:,:,n) = reshape(a0v,[bdpsi,bdpsi,d]);
        indexd = indexd+1;
        fomd(indexd) = real(fomdval);
        if iterd >= 2 && all(fomd((iterd-2)*n+1:iterd*n) > 0) && std(fomd((iterd-2)*n+1:iterd*n),1)/mean(fomd((iterd-2)*n+1:iterd*n)) <= reluncd
            break
        end
    end
    fomd = fomd(1:indexd);
    %plot(fomd)
    f_iter(2*iter) = fomd(iterd);
end
f_iter = f_iter(1:2*iter-1);
%plot(f_iter)
result = f_iter(2*iter-1);
end