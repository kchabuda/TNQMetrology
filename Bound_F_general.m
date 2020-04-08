function result = Bound_F_general(numof2box,d,integral0,integral1,integral0p)
if numof2box == 1/2
    error('Error occurred')
else
    n = numof2box+1;
end
if integral1 < 0
    phihalf = -integral1/(2*integral0);
else
    phihalf = integral1/integral0;
end
% Parameters
% Functions
function y = kroncoord(nt)
    coord = zeros([1,n]);
    for xf = 1:n
        coord(xf) = ceil(nt/d^(n-xf));
        nt = nt-d^(n-xf)*(coord(xf)-1);
    end
    y = coord;
end
function y = channel(nv,nvp)
    if numof2box == 1
        y = an(1,:,nv(1),nvp(1),1)*an(:,1,nv(2),nvp(2),3);
    elseif numof2box == 2
        y = an(1,:,nv(1),nvp(1),1)*an(:,:,nv(2),nvp(2),2)*an(:,1,nv(3),nvp(3),3);
    elseif numof2box == 3
        y = an(1,:,nv(1),nvp(1),1)*an(:,:,nv(2),nvp(2),2)*an(:,:,nv(3),nvp(3),2)*an(:,1,nv(4),nvp(4),3);
    end
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
kroncoordm = zeros([n,d^n]);
for nt = 1:d^n
    kroncoordm(:,nt) = kroncoord(nt);
end
nxpm = repmat(1:d,[d,1]);
nxm = nxpm.';
nxminusnxp = nxm-nxpm;
bdnoise1 = 1;
an1 = exp(-nxminusnxp.^2*integral0/2);
an1 = permute(an1,[3,4,1,2]);
an1halfleft = exp(-nxminusnxp.^2*integral0p/2);
an1halfleft = permute(an1halfleft,[3,4,1,2]);
an1halfright = exp(-nxminusnxp.^2*(integral0-integral0p)/2);
an1halfright = permute(an1halfright,[3,4,1,2]);
bdnoise2 = 2*d-1;
corr12 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral1);
[u,s,v] = svds(corr12,bdnoise2);
s = diag(s);
us = u*diag(s.^(1/2));
sv = diag(s.^(1/2))*v';
an2 = zeros([bdnoise2,bdnoise2,d,d,3]);
for nx = 1:d
    for nxp = 1:d
        an2(1,:,nx,nxp,1) = us(d*(nx-1)+nxp,:);
        an2(:,:,nx,nxp,2) = sv(:,d*(nx-1)+nxp)*us(d*(nx-1)+nxp,:);
        an2(:,1,nx,nxp,3) = sv(:,d*(nx-1)+nxp);
    end
end
bdnoise = bdnoise1*bdnoise2;
an = zeros([bdnoise,bdnoise,d,d,3]);
for nx = 1:d
    for nxp = 1:d
        an(:,:,nx,nxp,1) = an1halfleft(:,:,nx,nxp)*an2(:,:,nx,nxp,1);
        an(:,:,nx,nxp,2) = an1(:,:,nx,nxp)*an2(:,:,nx,nxp,2);
        an(:,:,nx,nxp,3) = an1halfright(:,:,nx,nxp)*an2(:,:,nx,nxp,3);
    end
end
channel2 = zeros(d^n);
for nt = 1:d^n
    for ntp = 1:d^n
        channel2(nt,ntp) = channel(kroncoordm(:,nt),kroncoordm(:,ntp));
    end
end
htot = zeros(d^n);
for x = 1:n
    if x == 1
        htot = htot+h(x)*phihalf;
    elseif x == n
        htot = htot+h(x)*(1-phihalf);
    else
        htot = htot+h(x);
    end
end
superoperator = diag(channel2(:));
dynamicalmap = reshape(permute(reshape(superoperator,[d^n,d^n,d^n,d^n]),[1,3,2,4]),[d^(2*n),d^(2*n)]);
dynamicalmap = (dynamicalmap+dynamicalmap')/2;
[eigvec,eigval] = eig(dynamicalmap);
eigval = diag(eigval);
if sum(eigval < -10^-10) > 0
    result = Inf;
    return
end
kraus_number_max = sum(eigval > 10^-10);
kraus = cell([1,kraus_number_max]);
dkraus = cell([1,kraus_number_max]);
kraus_number = 0;
for iter = 1:d^(2*n)
    if eigval(iter) > 10^-10
        kraus_number = kraus_number+1;
        kraus{kraus_number} = sqrt(eigval(iter))*reshape(eigvec(:,iter),[d^n,d^n]);
        dkraus{kraus_number} = -1i*htot*kraus{kraus_number};
    end
end
c = ChExt(kraus,dkraus);
f = (1/c^2)/numof2box;
result = f;
end