function result = Exact_Fid_OBC(n,d,A,B,phi)
% Functions
function y = kroncoord(nt)
    coord = zeros([1,n]);
    for xf = 1:n
        coord(xf) = ceil(nt/d^(n-xf));
        nt = nt-d^(n-xf)*(coord(xf)-1);
    end
    y = coord;
end
function y = rbarampo(nv,nvp)
    y = 1;
    for xf = 1:n
        y = y*A{xf}(:,:,nv(xf),nvp(xf));
    end
end
function y = rbarbmpo(nv,nvp)
    y = 1;
    for xf = 1:n
        y = y*B{xf}(:,:,nv(xf),nvp(xf));
    end
end
% Calculations
%  Exact
kroncoordm = zeros([n,d^n]);
for nt = 1:d^n
    kroncoordm(:,nt) = kroncoord(nt);
end
rbar = zeros(d^n);
for nt = 1:d^n
    for ntp = 1:d^n
        rbar(nt,ntp) = rbarampo(kroncoordm(:,nt),kroncoordm(:,ntp));
    end
end
rbar = (rbar+rbar')/2;
rbarphi = zeros(d^n);
for nt = 1:d^n
    for ntp = 1:d^n
        rbarphi(nt,ntp) = rbarbmpo(kroncoordm(:,nt),kroncoordm(:,ntp));
    end
end
rbarphi = (rbarphi+rbarphi')/2;
rbarp = (rbarphi-rbar)/phi;
[rbareigvec,rbareigval] = eig(rbar);
lpart1 = zeros(d^n);
for nt = 1:d^n
    for ntp = 1:d^n
        if abs(rbareigval(nt,nt)+rbareigval(ntp,ntp)) > 10^-5
            lpart1(nt,ntp) = 1/(rbareigval(nt,nt)+rbareigval(ntp,ntp));
        else
            lpart1(nt,ntp) = 0;
        end
    end
end
lpart2 = rbareigvec'*rbarp*rbareigvec;
l = rbareigvec*(2*lpart1.*lpart2)*rbareigvec';
result = real(trace(rbarp*l));
end