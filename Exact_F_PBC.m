function result = Exact_F_PBC(n,d,noiserange,integral0,integral1,integral2,imprecision)
%d = 2;
%noiserange = 2;
%integral0 = 1;
%integral1 = 0.1;
%integral2 = 0.01;
%imprecision = 10^-3;
relunc = 0.1*imprecision;
% Functions
function y = kroncoord(nt)
    coord = zeros([1,n]);
    for xf = 1:n
        coord(xf) = ceil(nt/d^(n-xf));
        nt = nt-d^(n-xf)*(coord(xf)-1);
    end
    y = coord;
end
function y = expintegral(nv,nvp)
    y = exp(-(nv-nvp).'*integralforr*(nv-nvp)/2);
end
function y = h(x)
    if x == 1
        y = 0:d-1;
    else
        y = ones([1,d]);
    end
    for xf = 2:n
        if x == xf
            y = kron(y,0:d-1);
        else
            y = kron(y,ones([1,d]));
        end
    end
    y = diag(y);
end
% Calculations
%  Exact
kroncoordm = zeros([n,d^n]);
for nt = 1:d^n
    kroncoordm(:,nt) = kroncoord(nt);
end
if noiserange == 0
    integralforr = zeros(n);
elseif noiserange == 1
    integralforr = integral0*eye(n);
elseif noiserange == 2
    integralforr = integral0*eye(n)+diag(integral1*ones([n-1,1]),1)+diag(integral1*ones([n-1,1]),-1);
    integralforr(1,n) = integral1;
    integralforr(n,1) = integral1;
elseif noiserange == 3    
    integralforr = integral0*eye(n)+diag(integral1*ones([n-1,1]),1)+diag(integral1*ones([n-1,1]),-1)+diag(integral2*ones([n-2,1]),2)+diag(integral2*ones([n-2,1]),-2);
    integralforr(1,n) = integral1;
    integralforr(n,1) = integral1;
    if n >= 4
        integralforr(1,n-1) = integral2;
        integralforr(2,n) = integral2;
        integralforr(n-1,1) = integral2;
        integralforr(n,2) = integral2;
    end
end
expintegral2 = zeros(d^n);
for nt = 1:d^n
    for ntp = 1:d^n
        expintegral2(nt,ntp) = expintegral(kroncoordm(:,nt),kroncoordm(:,ntp));
    end
end
htot = zeros(d^n);
for x = 1:n
    htot = htot+h(x);
end
psi0 = sqrt(2/(d+1))*sin((1:d)*pi/(d+1));
r0 = psi0'*psi0;
r0timesn = r0;
for x = 1:n-1
    r0timesn = kron(r0timesn,r0);
end
f = zeros([2*100,1]);
iter = 0;
while 1
    rbar = r0timesn.*expintegral2;
    rbarp = -1i*(htot*rbar-rbar*htot);
    rbar = (rbar+rbar')/2;
    [rbareigvec,rbareigval] = eig(rbar);
    lpart1 = zeros(d^n);
    for nt = 1:d^n
        for ntp = 1:d^n
            if abs(rbareigval(nt,nt)+rbareigval(ntp,ntp)) > 10^-10
                lpart1(nt,ntp) = 1/(rbareigval(nt,nt)+rbareigval(ntp,ntp));
            else
                lpart1(nt,ntp) = 0;
            end
        end
    end
    lpart2 = rbareigvec'*rbarp*rbareigvec;
    l = rbareigvec*(2*lpart1.*lpart2)*rbareigvec';
    iter = iter+1;
    f(iter) = real(trace(rbarp*l));
    if iter >= 5 && std(f(iter-4:iter),1)/mean(f(iter-4:iter)) <= relunc
        break
    end
    l2d = l^2.*expintegral2;
    lpd = 1i*(htot*l-l*htot).*expintegral2;
    eiginput = 2*lpd-l2d;
    eiginput = (eiginput+eiginput')/2;
    [psi,fval] = eigs(eiginput,1,'lr');
    iter = iter+1;
    f(iter) = real(fval);
    psi = psi/norm(psi);
    r0timesn = psi*psi';
end
f = f(1:iter);
%plot(f)
result = f(iter);
end