function [result,qAvar] = Exact_A_OBC_prod(m,d,ti,fullnoise,noiserange,alfa,gamma,beta)
%d = 2;
%fullnoise = 0;
%noiserange = 2;
%alfa = 1; gamma = 2; beta = 0.1;
n = 2*m-1;
if fullnoise == 1
    noiserange = n+1;
end
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
function y = sumintegral(nv,nvp)
    y = sum(integralforrp*(nv-nvp));
end
function y = sigma0(t,alfa,gamma,beta)
    y = alfa*(2*gamma*t+4*exp(-gamma*t)-exp(-2*gamma*t)-3)/(gamma*t)^2+beta/t;
end
% Calculations
%  Exact
kroncoordm = zeros([n,d^n]);
for nt = 1:d^n
    kroncoordm(:,nt) = kroncoord(nt);
end
integral0 = 2*alfa*(gamma*ti+exp(-gamma*ti)-1)/gamma^2+beta*ti;
integralforr0 = integral0*eye(n+1);
for i = 1:noiserange-1
    integral = 2*alfa*(cosh(gamma*ti)-1)*exp(-i*gamma*ti)/gamma^2;
    integralforr0 = integralforr0+diag(integral*ones([n+1-i,1]),i)+diag(integral*ones([n+1-i,1]),-i);
end
integralforr = integralforr0(1:n,1:n);
integralforrp = integralforr0(1:2*floor((n+1)/2),1:n);
integralforrp = integralforrp.*[-ones([floor((n+1)/2),n]);ones([floor((n+1)/2),n])];
expintegral2 = zeros(d^n);
for nt = 1:d^n
    for ntp = 1:d^n
        expintegral2(nt,ntp) = expintegral(kroncoordm(:,nt),kroncoordm(:,ntp));
    end
end
sumintegral2 = zeros(d^n);
for nt = 1:d^n
    for ntp = 1:d^n
        sumintegral2(nt,ntp) = sumintegral(kroncoordm(:,nt),kroncoordm(:,ntp));
    end
end
%{
psi0 = zeros([1,d]);
for i = 1:d
    psi0(i) = sqrt(nchoosek(d-1,i-1));
end
psi0 = psi0/2^((d-1)/2); %prod
%}
psi0 = sqrt(2/(d+1))*sin((1:d)*pi/(d+1)); %sine
%psi0 = [1,zeros([1,d-2]),1]/sqrt(2); %N00N
r0 = psi0'*psi0;
r0timesn = r0;
for x = 1:n-1
    r0timesn = kron(r0timesn,r0);
end
rbar = r0timesn.*expintegral2;
rbarp = -(1i/sqrt(ti*m))*rbar.*sumintegral2;
rbar = (rbar+rbar')/2;
[rbareigvec,rbareigval] = eig(rbar);
if sum(diag(rbareigval) < -10^-10) > 0
    result = NaN;
    return
end
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
result = real(trace(rbarp*l))/2;
qAvar = sigma0(ti*m,alfa,gamma,beta)-result/(ti*m);
end