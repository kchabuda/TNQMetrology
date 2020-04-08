function [result,A0,c] = MPO_A_OBC(approx,figureofmerit,m,d,bdl,ti,noiserange,integral0,integral1,integral2,lherm,imprecision,A0,c)
% Parameters
if approx == 1
    n = m-1;
else
    n = 2*m-1;
end
tol1 = 0.1*imprecision/n^2;
tol2 = 0.1*imprecision/n^2;
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
function [u,s,v] = svd2(x)
    [s1,s2] = size(x);
    if s1 >= s2
        [u,s,v] = svd(x,0);
    else
        [v,s,u] = svd(x',0);
    end
    v = v';
end
% Calculations
%  MPO
for x = n:-1:2
    [d1,d2,d3] = size(A0{x});
    A0{x} = permute(A0{x},[1,3,2]);
    A0{x} = reshape(A0{x},[d1,d3*d2]);
    [u,s,v] = svd2(A0{x});
    A0{x} = reshape(v,[size(s,1),d3,d2]);
    A0{x} = permute(A0{x},[1,3,2]);
    tensors = {A0{x-1},u*s};
    legs = {[-1,1,-3],[1,-2]};
    A0{x-1} = ncon(tensors,legs);
end
if noiserange == 0
    AN = cell([1,5]);
    [AN{:}] = deal(ones([1,1,d,d]));
elseif noiserange == 1
    nxpm = repmat(1:d,[d,1]);
    nxm = nxpm.';
    nxminusnxp = nxm-nxpm;
    AN = cell([1,5]);
	AN{1} = exp(-nxminusnxp.^2*integral0/2);
    AN{1} = permute(AN{1},[3,4,1,2]);
    [AN{:}] = deal(AN{1});
elseif noiserange == 2
    nxpm = repmat(1:d,[d,1]);
    nxm = nxpm.';
    nxminusnxp = nxm-nxpm;
    AN1 = cell(1);
    AN1{1} = exp(-nxminusnxp.^2*integral0/2);
    AN1{1} = permute(AN1{1},[3,4,1,2]);
    corr12 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral1);
    [u,s,v] = svds(corr12,2*d-1);
    s = diag(s);
    us = u*diag(s.^(1/2));
    sv = diag(s.^(1/2))*v';
    AN2 = cell([1,3]);
    for nx = 1:d
        for nxp = 1:d
            AN2{1}(:,:,nx,nxp) = us(d*(nx-1)+nxp,:);
            AN2{2}(:,:,nx,nxp) = sv(:,d*(nx-1)+nxp)*us(d*(nx-1)+nxp,:);
            AN2{3}(:,:,nx,nxp) = sv(:,d*(nx-1)+nxp);
        end
    end
    AN = cell([1,5]);
    for nx = 1:d
        for nxp = 1:d
            AN{1}(:,:,nx,nxp) = AN1{1}(:,:,nx,nxp)*AN2{1}(:,:,nx,nxp);
            AN{2}(:,:,nx,nxp) = AN1{1}(:,:,nx,nxp)*AN2{2}(:,:,nx,nxp);
            AN{3}(:,:,nx,nxp) = AN1{1}(:,:,nx,nxp)*AN2{2}(:,:,nx,nxp);
            AN{4}(:,:,nx,nxp) = AN1{1}(:,:,nx,nxp)*AN2{2}(:,:,nx,nxp);
            AN{5}(:,:,nx,nxp) = AN1{1}(:,:,nx,nxp)*AN2{3}(:,:,nx,nxp);
        end
    end
elseif noiserange == 3
    nxpm = repmat(1:d,[d,1]);
    nxm = nxpm.';
    nxminusnxp = nxm-nxpm;
    AN1 = cell(1);
    AN1{1} = exp(-nxminusnxp.^2*integral0/2);
    AN1{1} = permute(AN1{1},[3,4,1,2]);
    corr12 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral1);
    [u,s,v] = svds(corr12,2*d-1);
    s = diag(s);
    us = u*diag(s.^(1/2));
    sv = diag(s.^(1/2))*v';
    AN2 = cell([1,3]);
    for nx = 1:d
        for nxp = 1:d
            AN2{1}(:,:,nx,nxp) = us(d*(nx-1)+nxp,:);
            AN2{2}(:,:,nx,nxp) = sv(:,d*(nx-1)+nxp)*us(d*(nx-1)+nxp,:);
            AN2{3}(:,:,nx,nxp) = sv(:,d*(nx-1)+nxp);
        end
    end
    if n == 2
        AN3 = cell([1,5]);
        [AN3{:}] = deal(ones([1,1,d,d]));
    elseif n == 3
        corr13 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral2);
        [u,s,v] = svds(corr13,2*d-1);
        s = diag(s);
        us = u*diag(s.^(1/2));
        sv = diag(s.^(1/2))*v';
        AN3 = cell([1,5]);
        for nx = 1:d
            for nxp = 1:d
                AN3{1}(:,:,nx,nxp) = us(d*(nx-1)+nxp,:);
                AN3{3}(:,:,nx,nxp) = eye(2*d-1);
                AN3{5}(:,:,nx,nxp) = sv(:,d*(nx-1)+nxp);
            end
        end
        AN3{2} = AN3{3};
        AN3{4} = AN3{3};
    else
        corr13 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral2);
        [u,s,v] = svds(corr13,2*d-1);
        s = diag(s);
        us = u*diag(s.^(1/2));
        sv = diag(s.^(1/2))*v';
        AN3 = cell([1,5]);
        for nx = 1:d
            for nxp = 1:d
                AN3{1}(:,:,nx,nxp) = us(d*(nx-1)+nxp,:);
                AN3{3}(:,:,nx,nxp) = sv(:,d*(nx-1)+nxp)*us(d*(nx-1)+nxp,:);
                AN3{5}(:,:,nx,nxp) = sv(:,d*(nx-1)+nxp);
            end
        end
        tensors = {AN3{1},eye(2*d-1)};
        legs = {[-6,-3,-4,-5],[-1,-2]};
        AN3{2} = ncon(tensors,legs);
        AN3{2} = reshape(AN3{2},[2*d-1,(2*d-1)^2,d,d]);
        tensors = {AN3{3},eye(2*d-1)};
        legs = {[-1,-4,-5,-6],[-2,-3]};
        AN3{3} = ncon(tensors,legs);
        AN3{3} = reshape(AN3{3},[(2*d-1)^2,(2*d-1)^2,d,d]);
        tensors = {AN3{5},eye(2*d-1)};
        legs = {[-1,-6,-4,-5],[-2,-3]};
        AN3{4} = ncon(tensors,legs);
        AN3{4} = reshape(AN3{4},[(2*d-1)^2,2*d-1,d,d]);
    end
    AN = cell([1,5]);
    for nx = 1:d
        for nxp = 1:d
            AN{1}(:,:,nx,nxp) = AN1{1}(:,:,nx,nxp)*kron(AN2{1}(:,:,nx,nxp),AN3{1}(:,:,nx,nxp));
            AN{2}(:,:,nx,nxp) = AN1{1}(:,:,nx,nxp)*kron(AN2{2}(:,:,nx,nxp),AN3{2}(:,:,nx,nxp));
            AN{3}(:,:,nx,nxp) = AN1{1}(:,:,nx,nxp)*kron(AN2{2}(:,:,nx,nxp),AN3{3}(:,:,nx,nxp));
            AN{4}(:,:,nx,nxp) = AN1{1}(:,:,nx,nxp)*kron(AN2{2}(:,:,nx,nxp),AN3{4}(:,:,nx,nxp));
            AN{5}(:,:,nx,nxp) = AN1{1}(:,:,nx,nxp)*kron(AN2{3}(:,:,nx,nxp),AN3{5}(:,:,nx,nxp));
        end
    end
end
if noiserange == 0
    integralforr0 = zeros(n+1);
elseif noiserange == 1
    integralforr0 = integral0*eye(n+1);
elseif noiserange == 2
    integralforr0 = integral0*eye(n+1)+diag(integral1*ones([n,1]),1)+diag(integral1*ones([n,1]),-1);
elseif noiserange == 3
    integralforr0 = integral0*eye(n+1)+diag(integral1*ones([n,1]),1)+diag(integral1*ones([n,1]),-1)+diag(integral2*ones([n-1,1]),2)+diag(integral2*ones([n-1,1]),-2);
end
integralforrp = integralforr0(1:2*floor((n+1)/2),1:n);
if approx ~= 1
    integralforrp = integralforrp.*[-ones([floor((n+1)/2),n]);ones([floor((n+1)/2),n])];
end
f_iter = zeros([2*100,1]);
iter = 0;
while 1
    iter = iter+1;
    A = cell([1,n]);
    for nx = 1:d
        for nxp = 1:d
            A{1}(:,:,nx,nxp) = kron(kron(A0{1}(:,:,nx),conj(A0{1}(:,:,nxp))),AN{1}(:,:,nx,nxp));
            if n ~= 2
                A{2}(:,:,nx,nxp) = kron(kron(A0{2}(:,:,nx),conj(A0{2}(:,:,nxp))),AN{2}(:,:,nx,nxp));
                for x = 3:n-2
                    A{x}(:,:,nx,nxp) = kron(kron(A0{x}(:,:,nx),conj(A0{x}(:,:,nxp))),AN{3}(:,:,nx,nxp));
                end
                A{n-1}(:,:,nx,nxp) = kron(kron(A0{n-1}(:,:,nx),conj(A0{n-1}(:,:,nxp))),AN{4}(:,:,nx,nxp));
            end
            A{n}(:,:,nx,nxp) = kron(kron(A0{n}(:,:,nx),conj(A0{n}(:,:,nxp))),AN{5}(:,:,nx,nxp));
        end
    end
    B = cell([1,n]);
    for nx = 1:d
        for nxp = 1:d
            B{1}(:,:,nx,nxp) = kron([-(1i/sqrt(ti*m))*(nx-nxp)*sum(integralforrp(:,1)),1],A{1}(:,:,nx,nxp));
            for x = 2:n-1
                B{x}(:,:,nx,nxp) = kron([1,0;-(1i/sqrt(ti*m))*(nx-nxp)*sum(integralforrp(:,x)),1],A{x}(:,:,nx,nxp));
            end
            B{n}(:,:,nx,nxp) = kron([1;-(1i/sqrt(ti*m))*(nx-nxp)*sum(integralforrp(:,n))],A{n}(:,:,nx,nxp));
        end
    end
%  2Tr(r'L)-Tr(rLL) maximization
    if figureofmerit == 1
        L1F = cell([1,n-1]);
        L2F = cell([1,n-1]);
        fom1 = zeros([n*100,1]);
        index = 0;
        iter1 = 0;
        while 1
            iter1 = iter1+1;
            tensors = {c(:,1,:,:,n),B{n}};
            legs = {[-1,-3,1,2],[-2,-4,2,1]};
            L1F{n-1} = ncon(tensors,legs);
            tensors = {c(:,1,:,:,n),A{n},c(:,1,:,:,n)};
            legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
            L2F{n-1} = ncon(tensors,legs);
            for x = n-1:-1:2
                tensors = {c(:,:,:,:,x),B{x},L1F{x}};
                legs = {[-1,3,1,2],[-2,4,2,1],[3,4]};
                L1F{x-1} = ncon(tensors,legs);
                tensors = {c(:,:,:,:,x),A{x},c(:,:,:,:,x),L2F{x}};
                legs = {[-1,4,1,2],[-2,5,2,3],[-3,6,3,1],[4,5,6]};
                L2F{x-1} = ncon(tensors,legs);
            end
            tensors = {B{1},L1F{1}};
            legs = {[-5,1,-4,-3],[-2,1,-1]};
            l1 = ncon(tensors,legs);
            l1 = reshape(l1,[bdl*d*d,1]);
            tensors = {A{1},eye(d),L2F{1}};
            legs = {[-9,1,-4,-7],[-8,-3],[-2,1,-6,-1,-5]};
            l2 = ncon(tensors,legs);
            l2 = reshape(l2,[bdl*d*d,bdl*d*d]);
            dl2 = l2+l2.';
            dl1 = 2*l1;
            dl2pinv = pinv2(dl2,tol1);
            dl2pinv = (dl2pinv+dl2pinv.')/2;
            cv = dl2pinv*dl1;
            c(1,:,:,:,1) = reshape(cv,[1,bdl,d,d]);
            if lherm == 1
                c(:,:,:,:,1) = (c(:,:,:,:,1)+conj(permute(c(:,:,:,:,1),[1,2,4,3])))/2;
                cv = reshape(c(1,:,:,:,1),[bdl*d*d,1]);
            end
            index = index+1;
            fom1(index) = real(2*cv.'*l1-cv.'*l2*cv);
            tensors = {c(1,:,:,:,1),B{1}};
            legs = {[-3,-1,1,2],[-4,-2,2,1]};
            l1c = ncon(tensors,legs);
            tensors = {c(1,:,:,:,1),A{1},c(1,:,:,:,1)};
            legs = {[-4,-1,1,2],[-5,-2,2,3],[-6,-3,3,1]};
            l2c = ncon(tensors,legs);
            for x = 2:n-1
                tensors = {l1c,B{x},L1F{x}};
                legs = {[-1,1],[1,2,-4,-3],[-2,2]};
                l1 = ncon(tensors,legs);
                l1 = reshape(l1,[bdl*bdl*d*d,1]);
                tensors = {l2c,A{x},eye(d),L2F{x}};
                legs = {[-1,1,-5],[1,2,-4,-7],[-8,-3],[-2,2,-6]};
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
                tensors = {l1c,c(:,:,:,:,x),B{x}};
                legs = {[3,4],[3,-1,1,2],[4,-2,2,1]};
                l1c = ncon(tensors,legs);
                tensors = {l2c,c(:,:,:,:,x),A{x},c(:,:,:,:,x)};
                legs = {[4,5,6],[4,-1,1,2],[5,-2,2,3],[6,-3,3,1]};
                l2c = ncon(tensors,legs);
            end
            tensors = {l1c,B{n}};
            legs = {[-1,1,-2],[1,-5,-4,-3]};
            l1 = ncon(tensors,legs);
            l1 = reshape(l1,[bdl*d*d,1]);
            tensors = {l2c,A{n},eye(d)};
            legs = {[-1,1,-5,-2,-6,],[1,-9,-4,-7],[-8,-3]};
            l2 = ncon(tensors,legs);
            l2 = reshape(l2,[bdl*d*d,bdl*d*d]);
            dl2 = l2+l2.';
            dl1 = 2*l1;
            dl2pinv = pinv2(dl2,tol1);
            dl2pinv = (dl2pinv+dl2pinv.')/2;
            cv = dl2pinv*dl1;
            c(:,1,:,:,n) = reshape(cv,[bdl,1,d,d]);
            if lherm == 1
                c(:,:,:,:,n) = (c(:,:,:,:,n)+conj(permute(c(:,:,:,:,n),[1,2,4,3])))/2;
                cv = reshape(c(:,1,:,:,n),[bdl*d*d,1]);
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
        tensors = {B{n},B{n}};
        legs = {[-1,-3,1,2],[-2,-4,2,1]};
        l0 = ncon(tensors,legs);
        for x = n-1:-1:2
            tensors = {B{x},B{x},l0};
            legs = {[-1,3,1,2],[-2,4,2,1],[3,4]};
            l0 = ncon(tensors,legs);
        end
        tensors = {B{1},B{1},l0};
        legs = {[-1,3,1,2],[-2,4,2,1],[3,4]};
        l0 = ncon(tensors,legs);
        l0 = abs(l0);
        if lherm == 1
            L1_1F = cell([1,n-1]);
            L1_2F = cell([1,n-1]);
            L2_1F = cell([1,n-1]);
            L2_2F = cell([1,n-1]);
            fom2 = zeros([n*100,1]);
            index = 0;
            iter2 = 0;
            while 1
                iter2 = iter2+1;
                tensors = {c(:,1,:,:,n),B{n},A{n}};
                legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                L1_1F{n-1} = ncon(tensors,legs);
                tensors = {c(:,1,:,:,n),A{n},B{n}};
                legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                L1_2F{n-1} = ncon(tensors,legs);
                tensors = {conj(c(:,1,:,:,n)),A{n},c(:,1,:,:,n),A{n}};
                legs = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                L2_1F{n-1} = ncon(tensors,legs);
                tensors = {c(:,1,:,:,n),c(:,1,:,:,n),A{n},A{n}};
                legs = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                L2_2F{n-1} = ncon(tensors,legs);
                for x = n-1:-1:2
                    tensors = {c(:,:,:,:,x),B{x},A{x},L1_1F{x}};
                    legs = {[-1,4,1,2],[-2,5,2,3],[-3,6,3,1],[4,5,6]};
                    L1_1F{x-1} = ncon(tensors,legs);
                    tensors = {c(:,:,:,:,x),A{x},B{x},L1_2F{x}};
                    legs = {[-1,4,1,2],[-2,5,2,3],[-3,6,3,1],[4,5,6]};
                    L1_2F{x-1} = ncon(tensors,legs);
                    tensors = {conj(c(:,:,:,:,x)),A{x},c(:,:,:,:,x),A{x},L2_1F{x}};
                    legs = {[-1,5,2,1],[-2,6,2,3],[-3,7,3,4],[-4,8,4,1],[5,6,7,8]};
                    L2_1F{x-1} = ncon(tensors,legs);
                    tensors = {c(:,:,:,:,x),c(:,:,:,:,x),A{x},A{x},L2_2F{x}};
                    legs = {[-1,5,1,2],[-2,6,2,3],[-3,7,3,4],[-4,8,4,1],[5,6,7,8]};
                    L2_2F{x-1} = ncon(tensors,legs);
                end
                tensors = {B{1},A{1},L1_1F{1}};
                legs = {[-5,2,-4,1],[-6,3,1,-3],[-2,2,3,-1]};
                l1_1 = ncon(tensors,legs);
                l1_1 = reshape(l1_1,[bdl*d*d,1]);
                tensors = {A{1},B{1},L1_2F{1}};
                legs = {[-5,2,-4,1],[-6,3,1,-3],[-2,2,3,-1]};
                l1_2 = ncon(tensors,legs);
                l1_2 = reshape(l1_2,[bdl*d*d,1]);
                tensors = {A{1},A{1},L2_1F{1}};
                legs = {[-9,1,-4,-7],[-10,2,-8,-3],[-2,1,-6,2,-1,-5]};
                l2_1 = ncon(tensors,legs);
                l2_1 = reshape(l2_1,[bdl*d*d,bdl*d*d]);
                tensors = {eye(d),A{1},A{1},L2_2F{1}};
                legs = {[-4,-7],[-9,2,-8,1],[-10,3,1,-3],[-2,-6,2,3,-1,-5]};
                l2_2 = ncon(tensors,legs);
                l2_2 = reshape(l2_2,[bdl*d*d,bdl*d*d]);
                dl2 = l2_1+l2_1.'+l2_2+l2_2.';
                dl1 = 2*(l1_1+l1_2);
                dl2pinv = pinv2(dl2,tol2);
                dl2pinv = (dl2pinv+dl2pinv.')/2;
                cv = dl2pinv*dl1;
                c(1,:,:,:,1) = reshape(cv,[1,bdl,d,d]);
                c(:,:,:,:,1) = (c(:,:,:,:,1)+conj(permute(c(:,:,:,:,1),[1,2,4,3])))/2;
                cv = reshape(c(1,:,:,:,1),[bdl*d*d,1]);
                index = index+1;
                fom2(index) = abs(1+(-4*(cv.'*l1_1+cv.'*l1_2)+2*(cv.'*l2_1*cv+cv.'*l2_2*cv))/(4*l0));
                tensors = {c(1,:,:,:,1),B{1},A{1}};
                legs = {[-4,-1,1,2],[-5,-2,2,3],[-6,-3,3,1]};
                l1_1c = ncon(tensors,legs);
                tensors = {c(1,:,:,:,1),A{1},B{1}};
                legs = {[-4,-1,1,2],[-5,-2,2,3],[-6,-3,3,1]};
                l1_2c = ncon(tensors,legs);
                tensors = {c(1,:,:,:,1),A{1},c(1,:,:,:,1),A{1}};
                legs = {[-5,-1,1,2],[-6,-2,2,3],[-7,-3,3,4],[-8,-4,4,1]};
                l2_1c = ncon(tensors,legs);
                tensors = {c(1,:,:,:,1),c(1,:,:,:,1),A{1},A{1}};
                legs = {[-5,-1,1,2],[-6,-2,2,3],[-7,-3,3,4],[-8,-4,4,1]};
                l2_2c = ncon(tensors,legs);
                for x = 2:n-1
                    tensors = {l1_1c,B{x},A{x},L1_1F{x}};
                    legs = {[-1,2,3],[2,4,-4,1],[3,5,1,-3],[-2,4,5]};
                    l1_1 = ncon(tensors,legs);
                    l1_1 = reshape(l1_1,[bdl*bdl*d*d,1]);
                    tensors = {l1_2c,A{x},B{x},L1_2F{x}};
                    legs = {[-1,2,3],[2,4,-4,1],[3,5,1,-3],[-2,4,5]};
                    l1_2 = ncon(tensors,legs);
                    l1_2 = reshape(l1_2,[bdl*bdl*d*d,1]);
                    tensors = {l2_1c,A{x},A{x},L2_1F{x}};
                    legs = {[-1,1,-5,2],[1,3,-4,-7],[2,4,-8,-3],[-2,3,-6,4]};
                    l2_1 = ncon(tensors,legs);
                    l2_1 = reshape(l2_1,[bdl*bdl*d*d,bdl*bdl*d*d]);
                    tensors = {l2_2c,eye(d),A{x},A{x},L2_2F{x}};
                    legs = {[-1,-5,2,3],[-4,-7],[2,4,-8,1],[3,5,1,-3],[-2,-6,4,5]};
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
                    tensors = {l1_1c,c(:,:,:,:,x),B{x},A{x}};
                    legs = {[4,5,6],[4,-1,1,2],[5,-2,2,3],[6,-3,3,1]};
                    l1_1c = ncon(tensors,legs);
                    tensors = {l1_2c,c(:,:,:,:,x),A{x},B{x}};
                    legs = {[4,5,6],[4,-1,1,2],[5,-2,2,3],[6,-3,3,1]};
                    l1_2c = ncon(tensors,legs);
                    tensors = {l2_1c,c(:,:,:,:,x),A{x},c(:,:,:,:,x),A{x}};
                    legs = {[5,6,7,8],[5,-1,1,2],[6,-2,2,3],[7,-3,3,4],[8,-4,4,1]};
                    l2_1c = ncon(tensors,legs);
                    tensors = {l2_2c,c(:,:,:,:,x),c(:,:,:,:,x),A{x},A{x}};
                    legs = {[5,6,7,8],[5,-1,1,2],[6,-2,2,3],[7,-3,3,4],[8,-4,4,1]};
                    l2_2c = ncon(tensors,legs);
                end
                tensors = {l1_1c,B{n},A{n}};
                legs = {[-1,2,3,-2],[2,-5,-4,1],[3,-6,1,-3]};
                l1_1 = ncon(tensors,legs);
                l1_1 = reshape(l1_1,[bdl*d*d,1]);
                tensors = {l1_2c,A{n},B{n}};
                legs = {[-1,2,3,-2],[2,-5,-4,1],[3,-6,1,-3]};
                l1_2 = ncon(tensors,legs);
                l1_2 = reshape(l1_2,[bdl*d*d,1]);
                tensors = {l2_1c,A{n},A{n}};
                legs = {[-1,1,-5,2,-2,-6],[1,-9,-4,-7],[2,-10,-8,-3]};
                l2_1 = ncon(tensors,legs);
                l2_1 = reshape(l2_1,[bdl*d*d,bdl*d*d]);
                tensors = {l2_2c,eye(d),A{n},A{n}};
                legs = {[-1,-5,2,3,-2,-6],[-4,-7],[2,-9,-8,1],[3,-10,1,-3]};
                l2_2 = ncon(tensors,legs);
                l2_2 = reshape(l2_2,[bdl*d*d,bdl*d*d]);
                dl2 = l2_1+l2_1.'+l2_2+l2_2.';
                dl1 = 2*(l1_1+l1_2);
                dl2pinv = pinv2(dl2,tol2);
                dl2pinv = (dl2pinv+dl2pinv.')/2;
                cv = dl2pinv*dl1;
                c(:,1,:,:,n) = reshape(cv,[bdl,1,d,d]);
                c(:,:,:,:,n) = (c(:,:,:,:,n)+conj(permute(c(:,:,:,:,n),[1,2,4,3])))/2;
                cv = reshape(c(:,1,:,:,n),[bdl*d*d,1]);
                index = index+1;
                fom2(index) = abs(1+(-4*(cv.'*l1_1+cv.'*l1_2)+2*(cv.'*l2_1*cv+cv.'*l2_2*cv))/(4*l0));
                if fom2(index) < 10^-10 || iter2 >= 20 || (iter2 >= 2 && std(fom2((iter2-2)*n+1:iter2*n),1)/mean(fom2((iter2-2)*n+1:iter2*n)) <= relunc2)
                    break
                end
            end
        else
            L1_1F = cell([1,n-1]);
            L1_2F = cell([1,n-1]);
            L1_3F = cell([1,n-1]);
            L1_4F = cell([1,n-1]);
            L2_1F = cell([1,n-1]);
            L2_2F = cell([1,n-1]);
            L2_3F = cell([1,n-1]);
            fom2 = zeros([n*100,1]);
            index = 0;
            iter2 = 0;
            while 1
                iter2 = iter2+1;
                tensors = {conj(c(:,1,:,:,n)),B{n},A{n}};
                legs = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
                L1_1F{n-1} = ncon(tensors,legs);
                tensors = {conj(c(:,1,:,:,n)),A{n},B{n}};
                legs = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
                L1_2F{n-1} = ncon(tensors,legs);
                tensors = {c(:,1,:,:,n),B{n},A{n}};
                legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                L1_3F{n-1} = ncon(tensors,legs);
                tensors = {c(:,1,:,:,n),A{n},B{n}};
                legs = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                L1_4F{n-1} = ncon(tensors,legs);
                tensors = {conj(c(:,1,:,:,n)),A{n},c(:,1,:,:,n),A{n}};
                legs = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                L2_1F{n-1} = ncon(tensors,legs);
                tensors = {conj(c(:,1,:,:,n)),A{n},A{n},c(:,1,:,:,n)};
                legs = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                L2_2F{n-1} = ncon(tensors,legs);
                tensors = {conj(c(:,1,:,:,n)),c(:,1,:,:,n),A{n},A{n}};
                legs = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                L2_3F{n-1} = ncon(tensors,legs);
                for x = n-1:-1:2
                    tensors = {conj(c(:,:,:,:,x)),B{x},A{x},L1_1F{x}};
                    legs = {[-1,4,2,1],[-2,5,2,3],[-3,6,3,1],[4,5,6]};
                    L1_1F{x-1} = ncon(tensors,legs);
                    tensors = {conj(c(:,:,:,:,x)),A{x},B{x},L1_2F{x}};
                    legs = {[-1,4,2,1],[-2,5,2,3],[-3,6,3,1],[4,5,6]};
                    L1_2F{x-1} = ncon(tensors,legs);
                    tensors = {c(:,:,:,:,x),B{x},A{x},L1_3F{x}};
                    legs = {[-1,4,1,2],[-2,5,2,3],[-3,6,3,1],[4,5,6]};
                    L1_3F{x-1} = ncon(tensors,legs);
                    tensors = {c(:,:,:,:,x),A{x},B{x},L1_4F{x}};
                    legs = {[-1,4,1,2],[-2,5,2,3],[-3,6,3,1],[4,5,6]};
                    L1_4F{x-1} = ncon(tensors,legs);
                    tensors = {conj(c(:,:,:,:,x)),A{x},c(:,:,:,:,x),A{x},L2_1F{x}};
                    legs = {[-1,5,2,1],[-2,6,2,3],[-3,7,3,4],[-4,8,4,1],[5,6,7,8]};
                    L2_1F{x-1} = ncon(tensors,legs);
                    tensors = {conj(c(:,:,:,:,x)),A{x},A{x},c(:,:,:,:,x),L2_2F{x}};
                    legs = {[-1,5,2,1],[-2,6,2,3],[-3,7,3,4],[-4,8,4,1],[5,6,7,8]};
                    L2_2F{x-1} = ncon(tensors,legs);
                    tensors = {conj(c(:,:,:,:,x)),c(:,:,:,:,x),A{x},A{x},L2_3F{x}};
                    legs = {[-1,5,2,1],[-2,6,2,3],[-3,7,3,4],[-4,8,4,1],[5,6,7,8]};
                    L2_3F{x-1} = ncon(tensors,legs);
                end
                tensors = {B{1},A{1},L1_1F{1}};
                legs = {[-5,2,-3,1],[-6,3,1,-4],[-2,2,3,-1]};
                l1_1 = ncon(tensors,legs);
                l1_1 = reshape(l1_1,[bdl*d*d,1]);
                tensors = {A{1},B{1},L1_2F{1}};
                legs = {[-5,2,-3,1],[-6,3,1,-4],[-2,2,3,-1]};
                l1_2 = ncon(tensors,legs);
                l1_2 = reshape(l1_2,[bdl*d*d,1]);
                tensors = {B{1},A{1},L1_3F{1}};
                legs = {[-5,2,-4,1],[-6,3,1,-3],[-2,2,3,-1]};
                l1_3 = ncon(tensors,legs);
                l1_3 = reshape(l1_3,[bdl*d*d,1]);
                tensors = {A{1},B{1},L1_4F{1}};
                legs = {[-5,2,-4,1],[-6,3,1,-3],[-2,2,3,-1]};
                l1_4 = ncon(tensors,legs);
                l1_4 = reshape(l1_4,[bdl*d*d,1]);
                tensors = {A{1},A{1},L2_1F{1}};
                legs = {[-9,1,-3,-7],[-10,2,-8,-4],[-2,1,-6,2,-1,-5]};
                l2_1 = ncon(tensors,legs);
                l2_1 = reshape(l2_1,[bdl*d*d,bdl*d*d]);
                tensors = {A{1},A{1},eye(d),L2_2F{1}};
                legs = {[-9,2,-3,1],[-10,3,1,-7],[-8,-4],[-2,2,3,-6,-1,-5]};
                l2_2 = ncon(tensors,legs);
                l2_2 = reshape(l2_2,[bdl*d*d,bdl*d*d]);
                tensors = {eye(d),A{1},A{1},L2_3F{1}};
                legs = {[-3,-7],[-9,2,-8,1],[-10,3,1,-4],[-2,-6,2,3,-1,-5]};
                l2_3 = ncon(tensors,legs);
                l2_3 = reshape(l2_3,[bdl*d*d,bdl*d*d]);
                dl2 = 2*l2_1+l2_2+l2_3;
                dl2 = (dl2+dl2')/2;
                dl1 = 2*(l1_1+l1_2);
                dl2pinv = pinv2(dl2,tol2);
                dl2pinv = (dl2pinv+dl2pinv')/2;
                cv = dl2pinv*dl1;
                c(1,:,:,:,1) = reshape(cv,[1,bdl,d,d]);
                index = index+1;
                fom2(index) = abs(1+(-2*(cv'*l1_1+cv'*l1_2+cv.'*l1_3+cv.'*l1_4)+2*cv'*l2_1*cv+cv'*l2_2*cv+cv'*l2_3*cv)/(4*l0));
                tensors = {conj(c(1,:,:,:,1)),B{1},A{1}};
                legs = {[-4,-1,2,1],[-5,-2,2,3],[-6,-3,3,1]};
                l1_1c = ncon(tensors,legs);
                tensors = {conj(c(1,:,:,:,1)),A{1},B{1}};
                legs = {[-4,-1,2,1],[-5,-2,2,3],[-6,-3,3,1]};
                l1_2c = ncon(tensors,legs);
                tensors = {c(1,:,:,:,1),B{1},A{1}};
                legs = {[-4,-1,1,2],[-5,-2,2,3],[-6,-3,3,1]};
                l1_3c = ncon(tensors,legs);
                tensors = {c(1,:,:,:,1),A{1},B{1}};
                legs = {[-4,-1,1,2],[-5,-2,2,3],[-6,-3,3,1]};
                l1_4c = ncon(tensors,legs);
                tensors = {conj(c(1,:,:,:,1)),A{1},c(1,:,:,:,1),A{1}};
                legs = {[-5,-1,2,1],[-6,-2,2,3],[-7,-3,3,4],[-8,-4,4,1]};
                l2_1c = ncon(tensors,legs);
                tensors = {conj(c(1,:,:,:,1)),A{1},A{1},c(1,:,:,:,1)};
                legs = {[-5,-1,2,1],[-6,-2,2,3],[-7,-3,3,4],[-8,-4,4,1]};
                l2_2c = ncon(tensors,legs);
                tensors = {conj(c(1,:,:,:,1)),c(1,:,:,:,1),A{1},A{1}};
                legs = {[-5,-1,2,1],[-6,-2,2,3],[-7,-3,3,4],[-8,-4,4,1]};
                l2_3c = ncon(tensors,legs);
                for x = 2:n-1
                    tensors = {l1_1c,B{x},A{x},L1_1F{x}};
                    legs = {[-1,2,3],[2,4,-3,1],[3,5,1,-4],[-2,4,5]};
                    l1_1 = ncon(tensors,legs);
                    l1_1 = reshape(l1_1,[bdl*bdl*d*d,1]);
                    tensors = {l1_2c,A{x},B{x},L1_2F{x}};
                    legs = {[-1,2,3],[2,4,-3,1],[3,5,1,-4],[-2,4,5]};
                    l1_2 = ncon(tensors,legs);
                    l1_2 = reshape(l1_2,[bdl*bdl*d*d,1]);
                    tensors = {l1_3c,B{x},A{x},L1_3F{x}};
                    legs = {[-1,2,3],[2,4,-4,1],[3,5,1,-3],[-2,4,5]};
                    l1_3 = ncon(tensors,legs);
                    l1_3 = reshape(l1_3,[bdl*bdl*d*d,1]);
                    tensors = {l1_4c,A{x},B{x},L1_4F{x}};
                    legs = {[-1,2,3],[2,4,-4,1],[3,5,1,-3],[-2,4,5]};
                    l1_4 = ncon(tensors,legs);
                    l1_4 = reshape(l1_4,[bdl*bdl*d*d,1]);
                    tensors = {l2_1c,A{x},A{x},L2_1F{x}};
                    legs = {[-1,1,-5,2],[1,3,-3,-7],[2,4,-8,-4],[-2,3,-6,4]};
                    l2_1 = ncon(tensors,legs);
                    l2_1 = reshape(l2_1,[bdl*bdl*d*d,bdl*bdl*d*d]);
                    tensors = {l2_2c,A{x},A{x},eye(d),L2_2F{x}};
                    legs = {[-1,2,3,-5],[2,4,-3,1],[3,5,1,-7],[-8,-4],[-2,4,5,-6]};
                    l2_2 = ncon(tensors,legs);
                    l2_2 = reshape(l2_2,[bdl*bdl*d*d,bdl*bdl*d*d]);
                    tensors = {l2_3c,eye(d),A{x},A{x},L2_3F{x}};
                    legs = {[-1,-5,2,3],[-3,-7],[2,4,-8,1],[3,5,1,-4],[-2,-6,4,5]};
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
                    tensors = {l1_1c,conj(c(:,:,:,:,x)),B{x},A{x}};
                    legs = {[4,5,6],[4,-1,2,1],[5,-2,2,3],[6,-3,3,1]};
                    l1_1c = ncon(tensors,legs);
                    tensors = {l1_2c,conj(c(:,:,:,:,x)),A{x},B{x}};
                    legs = {[4,5,6],[4,-1,2,1],[5,-2,2,3],[6,-3,3,1]};
                    l1_2c = ncon(tensors,legs);
                    tensors = {l1_3c,c(:,:,:,:,x),B{x},A{x}};
                    legs = {[4,5,6],[4,-1,1,2],[5,-2,2,3],[6,-3,3,1]};
                    l1_3c = ncon(tensors,legs);
                    tensors = {l1_4c,c(:,:,:,:,x),A{x},B{x}};
                    legs = {[4,5,6],[4,-1,1,2],[5,-2,2,3],[6,-3,3,1]};
                    l1_4c = ncon(tensors,legs);
                    tensors = {l2_1c,conj(c(:,:,:,:,x)),A{x},c(:,:,:,:,x),A{x}};
                    legs = {[5,6,7,8],[5,-1,2,1],[6,-2,2,3],[7,-3,3,4],[8,-4,4,1]};
                    l2_1c = ncon(tensors,legs);
                    tensors = {l2_2c,conj(c(:,:,:,:,x)),A{x},A{x},c(:,:,:,:,x)};
                    legs = {[5,6,7,8],[5,-1,2,1],[6,-2,2,3],[7,-3,3,4],[8,-4,4,1]};
                    l2_2c = ncon(tensors,legs);
                    tensors = {l2_3c,conj(c(:,:,:,:,x)),c(:,:,:,:,x),A{x},A{x}};
                    legs = {[5,6,7,8],[5,-1,2,1],[6,-2,2,3],[7,-3,3,4],[8,-4,4,1]};
                    l2_3c = ncon(tensors,legs);
                end
                tensors = {l1_1c,B{n},A{n}};
                legs = {[-1,2,3,-2],[2,-5,-3,1],[3,-6,1,-4]};
                l1_1 = ncon(tensors,legs);
                l1_1 = reshape(l1_1,[bdl*d*d,1]);
                tensors = {l1_2c,A{n},B{n}};
                legs = {[-1,2,3,-2],[2,-5,-3,1],[3,-6,1,-4]};
                l1_2 = ncon(tensors,legs);
                l1_2 = reshape(l1_2,[bdl*d*d,1]);
                tensors = {l1_3c,B{n},A{n}};
                legs = {[-1,2,3,-2],[2,-5,-4,1],[3,-6,1,-3]};
                l1_3 = ncon(tensors,legs);
                l1_3 = reshape(l1_3,[bdl*d*d,1]);
                tensors = {l1_4c,A{n},B{n}};
                legs = {[-1,2,3,-2],[2,-5,-4,1],[3,-6,1,-3]};
                l1_4 = ncon(tensors,legs);
                l1_4 = reshape(l1_4,[bdl*d*d,1]);
                tensors = {l2_1c,A{n},A{n}};
                legs = {[-1,1,-5,2,-2,-6],[1,-9,-3,-7],[2,-10,-8,-4]};
                l2_1 = ncon(tensors,legs);
                l2_1 = reshape(l2_1,[bdl*d*d,bdl*d*d]);
                tensors = {l2_2c,A{n},A{n},eye(d)};
                legs = {[-1,2,3,-5,-2,-6],[2,-9,-3,1],[3,-10,1,-7],[-8,-4]};
                l2_2 = ncon(tensors,legs);
                l2_2 = reshape(l2_2,[bdl*d*d,bdl*d*d]);
                tensors = {l2_3c,eye(d),A{n},A{n}};
                legs = {[-1,-5,2,3,-2,-6],[-3,-7],[2,-9,-8,1],[3,-10,1,-4]};
                l2_3 = ncon(tensors,legs);
                l2_3 = reshape(l2_3,[bdl*d*d,bdl*d*d]);
                dl2 = 2*l2_1+l2_2+l2_3;
                dl2 = (dl2+dl2')/2;
                dl1 = 2*(l1_1+l1_2);
                dl2pinv = pinv2(dl2,tol2);
                dl2pinv = (dl2pinv+dl2pinv')/2;
                cv = dl2pinv*dl1;
                c(:,1,:,:,n) = reshape(cv,[bdl,1,d,d]);
                index = index+1;
                fom2(index) = abs(1+(-2*(cv'*l1_1+cv'*l1_2+cv.'*l1_3+cv.'*l1_4)+2*cv'*l2_1*cv+cv'*l2_2*cv+cv'*l2_3*cv)/(4*l0));
                if fom2(index) < 10^-10 || iter2 >= 20 || (iter2 >= 2 && std(fom2((iter2-2)*n+1:iter2*n),1)/mean(fom2((iter2-2)*n+1:iter2*n)) <= relunc2)
                    break
                end
            end
        end
        fom2 = fom2(1:index);
        %plot(-log10(fom2))
        tensors = {B{n},c(:,1,:,:,n)};
        legs = {[-1,-3,1,2],[-2,-4,2,1]};
        f = ncon(tensors,legs);
        for x = n-1:-1:2
            tensors = {B{x},c(:,:,:,:,x),f};
            legs = {[-1,3,1,2],[-2,4,2,1],[3,4]};
            f = ncon(tensors,legs);
        end
        tensors = {B{1},c(1,:,:,:,1),f};
        legs = {[-1,3,1,2],[-2,4,2,1],[3,4]};
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
                ch(1,:,nx,nxp,1) = [c(1,:,nx,nxp,1),conj(c(1,:,nxp,nx,1))];
                ch(:,:,nx,nxp,2:n-1) = [c(:,:,nx,nxp,2:n-1),zeros([bdl,bdl,1,1,n-2]);zeros([bdl,bdl,1,1,n-2]),conj(c(:,:,nxp,nx,2:n-1))];
                ch(:,1,nx,nxp,n) = [c(:,1,nx,nxp,n);conj(c(:,1,nxp,nx,n))];
            end
        end
        ch = ch/2^(1/n);
    end
    c2 = zeros([(bdherm*bdl)^2,(bdherm*bdl)^2,d,d,n]);
    for nx = 1:d
        for nxp = 1:d
            for nxpp = 1:d
                c2(1,:,nx,nxp,1) = c2(1,:,nx,nxp,1)+kron(ch(1,:,nx,nxpp,1),ch(1,:,nxpp,nxp,1));
                for x = 2:n-1
                    c2(:,:,nx,nxp,x) = c2(:,:,nx,nxp,x)+kron(ch(:,:,nx,nxpp,x),ch(:,:,nxpp,nxp,x));
                end
                c2(:,1,nx,nxp,n) = c2(:,1,nx,nxp,n)+kron(ch(:,1,nx,nxpp,n),ch(:,1,nxpp,nxp,n));
            end
        end
    end
    cp = zeros([2*bdherm*bdl,2*bdherm*bdl,d,d,n]);
    for nx = 1:d
        for nxp = 1:d
            cp(1,:,nx,nxp,1) = kron([(1i/sqrt(ti*m))*(nx-nxp)*sum(integralforrp(:,1)),1],ch(1,:,nx,nxp,1));
            for x = 2:n-1
                cp(:,:,nx,nxp,x) = kron([1,0;(1i/sqrt(ti*m))*(nx-nxp)*sum(integralforrp(:,x)),1],ch(:,:,nx,nxp,x));
            end
            cp(:,1,nx,nxp,n) = kron([1;(1i/sqrt(ti*m))*(nx-nxp)*sum(integralforrp(:,n))],ch(:,1,nx,nxp,n));
        end
    end
    C2D = cell([1,n]);
    CPD = cell([1,n]);
    for nx = 1:d
        for nxp = 1:d
            C2D{1}(:,:,nx,nxp) = kron(AN{1}(:,:,nx,nxp),c2(1,:,nx,nxp,1));
            CPD{1}(:,:,nx,nxp) = kron(AN{1}(:,:,nx,nxp),cp(1,:,nx,nxp,1));
            if n ~= 2
                C2D{2}(:,:,nx,nxp) = kron(AN{2}(:,:,nx,nxp),c2(:,:,nx,nxp,2));
                CPD{2}(:,:,nx,nxp) = kron(AN{2}(:,:,nx,nxp),cp(:,:,nx,nxp,2));
                for x = 3:n-2
                    C2D{x}(:,:,nx,nxp) = kron(AN{3}(:,:,nx,nxp),c2(:,:,nx,nxp,x));
                    CPD{x}(:,:,nx,nxp) = kron(AN{3}(:,:,nx,nxp),cp(:,:,nx,nxp,x));
                end
                C2D{n-1}(:,:,nx,nxp) = kron(AN{4}(:,:,nx,nxp),c2(:,:,nx,nxp,n-1));
                CPD{n-1}(:,:,nx,nxp) = kron(AN{4}(:,:,nx,nxp),cp(:,:,nx,nxp,n-1));
            end
            C2D{n}(:,:,nx,nxp) = kron(AN{5}(:,:,nx,nxp),c2(:,1,nx,nxp,n));
            CPD{n}(:,:,nx,nxp) = kron(AN{5}(:,:,nx,nxp),cp(:,1,nx,nxp,n));
        end
    end
%  Tr[r0*X] maximization
    L2DF = cell([1,n-1]);
	LPDF = cell([1,n-1]);
    fomd = zeros([n*100,1]);
    indexd = 0;
    iterd = 0;
    while 1
        iterd = iterd+1;
        tensors = {conj(A0{n}),C2D{n},A0{n}};
        legs = {[-1,-4,1],[-2,-5,1,2],[-3,-6,2]};
        L2DF{n-1} = ncon(tensors,legs);
        tensors = {conj(A0{n}),CPD{n},A0{n}};
        legs = {[-1,-4,1],[-2,-5,1,2],[-3,-6,2]};
        LPDF{n-1} = ncon(tensors,legs);
        for x = n-1:-1:2
            tensors = {conj(A0{x}),C2D{x},A0{x},L2DF{x}};
            legs = {[-1,4,1],[-2,5,1,2],[-3,6,2],[4,5,6]};
            L2DF{x-1} = ncon(tensors,legs);
            tensors = {conj(A0{x}),CPD{x},A0{x},LPDF{x}};
            legs = {[-1,4,1],[-2,5,1,2],[-3,6,2],[4,5,6]};
            LPDF{x-1} = ncon(tensors,legs);
        end
        [d1,d2,d3] = size(A0{1});
        tensors = {C2D{1},L2DF{1}};
        legs = {[-7,1,-3,-6],[-2,1,-5,-1,-8,-4]};
        l2d = ncon(tensors,legs);
        l2d = reshape(l2d,[d1*d2*d3,d1*d2*d3]);
        tensors = {CPD{1},LPDF{1}};
        legs = {[-7,1,-3,-6],[-2,1,-5,-1,-8,-4]};
        lpd = ncon(tensors,legs);
        lpd = reshape(lpd,[d1*d2*d3,d1*d2*d3]);
        eiginput = 2*lpd-l2d;
        eiginput = (eiginput+eiginput')/2;
        [a0v,fomdval] = eigs(eiginput,1,'lr');
        a0v = a0v/norm(a0v);
        A0{1} = reshape(a0v,[d1,d2,d3]);
        indexd = indexd+1;
        fomd(indexd) = real(fomdval);
        A0{1} = permute(A0{1},[3,1,2]);
        A0{1} = reshape(A0{1},[d3*d1,d2]);
        [u,s,v] = svd2(A0{1});
        A0{1} = reshape(u,[d3,d1,size(s,1)]);
        A0{1} = permute(A0{1},[2,3,1]);
        tensors = {s*v,A0{2}};
        legs = {[-1,1],[1,-2,-3],};
        A0{2} = ncon(tensors,legs);
        tensors = {conj(A0{1}),C2D{1},A0{1}};
        legs = {[-4,-1,1],[-5,-2,1,2],[-6,-3,2]};
        l2dc = ncon(tensors,legs);
        tensors = {conj(A0{1}),CPD{1},A0{1}};
        legs = {[-4,-1,1],[-5,-2,1,2],[-6,-3,2]};
        lpdc = ncon(tensors,legs);
        for x = 2:n-1
            [d1,d2,d3] = size(A0{x});
            tensors = {l2dc,C2D{x},L2DF{x}};
            legs = {[-1,1,-4],[1,2,-3,-6],[-2,2,-5]};
            l2d = ncon(tensors,legs);
            l2d = reshape(l2d,[d1*d2*d3,d1*d2*d3]);
            tensors = {lpdc,CPD{x},LPDF{x}};
            legs = {[-1,1,-4],[1,2,-3,-6],[-2,2,-5]};
            lpd = ncon(tensors,legs);
            lpd = reshape(lpd,[d1*d2*d3,d1*d2*d3]);
            eiginput = 2*lpd-l2d;
            eiginput = (eiginput+eiginput')/2;
            [a0v,fomdval] = eigs(eiginput,1,'lr');
            a0v = a0v/norm(a0v);
            A0{x} = reshape(a0v,[d1,d2,d3]);
            indexd = indexd+1;
            fomd(indexd) = real(fomdval);
            A0{x} = permute(A0{x},[3,1,2]);
            A0{x} = reshape(A0{x},[d3*d1,d2]);
            [u,s,v] = svd2(A0{x});
            A0{x} = reshape(u,[d3,d1,size(s,1)]);
            A0{x} = permute(A0{x},[2,3,1]);
            tensors = {s*v,A0{x+1}};
            legs = {[-1,1],[1,-2,-3],};
            A0{x+1} = ncon(tensors,legs);
            tensors = {l2dc,conj(A0{x}),C2D{x},A0{x}};
            legs = {[3,4,5],[3,-1,1],[4,-2,1,2],[5,-3,2]};
            l2dc = ncon(tensors,legs);
            tensors = {lpdc,conj(A0{x}),CPD{x},A0{x}};
            legs = {[3,4,5],[3,-1,1],[4,-2,1,2],[5,-3,2]};
            lpdc = ncon(tensors,legs);
        end
        [d1,d2,d3] = size(A0{n});
        tensors = {l2dc,C2D{n}};
        legs = {[-1,1,-4,-2,-8,-5,],[1,-7,-3,-6]};
        l2d = ncon(tensors,legs);
        l2d = reshape(l2d,[d1*d2*d3,d1*d2*d3]);
        tensors = {lpdc,CPD{n}};
        legs = {[-1,1,-4,-2,-8,-5,],[1,-7,-3,-6]};
        lpd = ncon(tensors,legs);
        lpd = reshape(lpd,[d1*d2*d3,d1*d2*d3]);
        eiginput = 2*lpd-l2d;
        eiginput = (eiginput+eiginput')/2;
        [a0v,fomdval] = eigs(eiginput,1,'lr');
        a0v = a0v/norm(a0v);
        A0{n} = reshape(a0v,[d1,d2,d3]);
        indexd = indexd+1;
        fomd(indexd) = real(fomdval);
        for x = n:-1:2
            [d1,d2,d3] = size(A0{x});
            A0{x} = permute(A0{x},[1,3,2]);
            A0{x} = reshape(A0{x},[d1,d3*d2]);
            [u,s,v] = svd2(A0{x});
            A0{x} = reshape(v,[size(s,1),d3,d2]);
            A0{x} = permute(A0{x},[1,3,2]);
            tensors = {A0{x-1},u*s};
            legs = {[-1,1,-3],[1,-2]};
            A0{x-1} = ncon(tensors,legs);
        end
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
if approx ~= 1
    result = result/2;
end
end