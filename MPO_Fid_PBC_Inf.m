function result = MPO_Fid_PBC_Inf(a,b,phi)
% Parameters
figureofmerit = 1;
d = 2;
bdl = 1;
lherm = 1;
imprecision = 10^-(5/2);
tol1 = phi;
tol2 = phi;
d01 = 10^-1;
d02 = 10^-1;
h1 = 10^-10;
h2 = 10^-10;
relunc1 = 0.1*imprecision;
relunc2 = imprecision;
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
function [fom1val,c,cnew] = fom1f(lambda,c,flag,cnew)
    if lambda ~= 1
        if flag ~= 1
            tensorsf = {c,b};
            legsf = {[-1,-3,1,2],[-2,-4,2,1]};
            tmf = ncon(tensorsf,legsf);
            tmf = reshape(tmf,[bdl*bdrhob,bdl*bdrhob]);
            [tmvrf,~] = eigs(tmf,1);
            [tmvlf,~] = eigs(tmf.',1);
            tmvnorm = tmvlf.'*tmvrf;
            l1r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdrhob]);
            l1l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdrhob]);
            tensorsf = {c,a,c};
            legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
            tmf = ncon(tensorsf,legsf);
            tmf = reshape(tmf,[bdl*bdrhoa*bdl,bdl*bdrhoa*bdl]);
            [tmvrf,~] = eigs(tmf,1);
            [tmvlf,~] = eigs(tmf.',1);
            tmvnorm = tmvlf.'*tmvrf;
            l2r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdrhoa,bdl]);
            l2l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdrhoa,bdl]);
            tensorsf = {l1l,b,l1r};
            legsf = {[-1,2],[2,1,-4,-3],[-2,1]};
            l1 = ncon(tensorsf,legsf);
            l1 = reshape(l1,[bdl*bdl*d*d,1]);
            tensorsf = {l2l,a,eye(d),l2r};
            legsf = {[-1,2,-5],[2,1,-4,-7],[-8,-3],[-2,1,-6]};
            l2 = ncon(tensorsf,legsf);
            l2 = reshape(l2,[bdl*bdl*d*d,bdl*bdl*d*d]);
            dl2 = l2+l2.';
            dl1 = 2*l1;
            dl2pinv = pinv2(dl2,tol1);
            dl2pinv = (dl2pinv+dl2pinv.')/2;
            cv = dl2pinv*dl1;
            cnew = reshape(cv,[bdl,bdl,d,d]);
            if lherm == 1
                cnew = (cnew+conj(permute(cnew,[1,2,4,3])))/2;
            end
            tensorsf = {cnew};
            legsf = {[-1,-2,1,1]};
            tmf = ncon(tensorsf,legsf);
            ctrf = eigs(tmf,1);
            ctrf = real(ctrf);
            cnew = d*cnew/ctrf;
            if bdl == 1
                cnew = reshape(cnew,[d,d]);
                cnewmd = mean(diag(cnew));
                cnew = imag(cnew);
                cnew = (cnew+rot90(cnew,2).')/2;
                cnew = cnewmd*eye(d)+1i*cnew;
                cnew = reshape(cnew,[bdl,bdl,d,d]);
            else
                for nxf = 1:d
                    cnew(:,:,nxf,nxf) = zeros(bdl);
                    cnew(1,1,nxf,nxf) = 1;
                end
            end
        end
        c = (cnew*sin(lambda*pi)-c*cos(lambda*pi));
        if lherm == 1
            c = (c+conj(permute(c,[1,2,4,3])))/2;
        end
        tensorsf = {c};
        legsf = {[-1,-2,1,1]};
        tmf = ncon(tensorsf,legsf);
        ctrf = eigs(tmf,1);
        ctrf = real(ctrf);
        c = d*c/ctrf;
        if bdl == 1
            c = reshape(c,[d,d]);
            cmd = mean(diag(c));
            c = imag(c);
            c = (c+rot90(c,2).')/2;
            c = cmd*eye(d)+1i*c;
            c = reshape(c,[bdl,bdl,d,d]);
        else
            for nxf = 1:d
                c(:,:,nxf,nxf) = zeros(bdl);
                c(1,1,nxf,nxf) = 1;
            end
        end
    end
    tensorsf = {c,b};
    legsf = {[-1,-3,1,2],[-2,-4,2,1]};
    tmf = ncon(tensorsf,legsf);
    tmf = reshape(tmf,[bdl*bdrhob,bdl*bdrhob]);
    l1val = eigs(tmf,1);
    tensorsf = {c,a,c};
    legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
    tmf = ncon(tensorsf,legsf);
    tmf = reshape(tmf,[bdl*bdrhoa*bdl,bdl*bdrhoa*bdl]);
    l2val = eigs(tmf,1);
    fom1val = real((2*l1val-l2val-1)/phi^2);
end
function [fom2val,c,cnew] = fom2f(lambda,c,flag,cnew)
    if lherm == 1
        if lambda ~= 1
            if flag ~= 1
                tensorsf = {c,b,a};
                legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdrhob*bdrhoa,bdl*bdrhob*bdrhoa]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l1_1r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdrhob,bdrhoa]);
                l1_1l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdrhob,bdrhoa]);
                tensorsf = {c,a,b};
                legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdrhoa*bdrhob,bdl*bdrhoa*bdrhob]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l1_2r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdrhoa,bdrhob]);
                l1_2l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdrhoa,bdrhob]);
                tensorsf = {c,a,c,a};
                legsf = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdrhoa*bdl*bdrhoa,bdl*bdrhoa*bdl*bdrhoa]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l2_1r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdrhoa,bdl,bdrhoa]);
                l2_1l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdrhoa,bdl,bdrhoa]);
                tensorsf = {c,c,a,a};
                legsf = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdl*bdrhoa*bdrhoa,bdl*bdl*bdrhoa*bdrhoa]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l2_2r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdl,bdrhoa,bdrhoa]);
                l2_2l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdl,bdrhoa,bdrhoa]);
                tensorsf = {l1_1l,b,a,l1_1r};
                legsf = {[-1,4,5],[4,2,-4,1],[5,3,1,-3],[-2,2,3]};
                l1_1 = ncon(tensorsf,legsf);
                l1_1 = reshape(l1_1,[bdl*bdl*d*d,1]);
                tensorsf = {l1_2l,a,b,l1_2r};
                legsf = {[-1,4,5],[4,2,-4,1],[5,3,1,-3],[-2,2,3]};
                l1_2 = ncon(tensorsf,legsf);
                l1_2 = reshape(l1_2,[bdl*bdl*d*d,1]);
                tensorsf = {l2_1l,a,a,l2_1r};
                legsf = {[-1,2,-5,4],[2,1,-4,-7],[4,3,-8,-3],[-2,1,-6,3]};
                l2_1 = ncon(tensorsf,legsf);
                l2_1 = reshape(l2_1,[bdl*bdl*d*d,bdl*bdl*d*d]);
                tensorsf = {l2_2l,eye(d),a,a,l2_2r};
                legsf = {[-1,-5,4,5],[-4,-7],[4,2,-8,1],[5,3,1,-3],[-2,-6,2,3]};
                l2_2 = ncon(tensorsf,legsf);
                l2_2 = reshape(l2_2,[bdl*bdl*d*d,bdl*bdl*d*d]);
                dl2 = l2_1+l2_1.'+l2_2+l2_2.';
                dl1 = 2*(l1_1+l1_2);
                dl2pinv = pinv2(dl2,tol2);
                dl2pinv = (dl2pinv+dl2pinv.')/2;
                cv = dl2pinv*dl1;
                cnew = reshape(cv,[bdl,bdl,d,d]);
                cnew = (cnew+conj(permute(cnew,[1,2,4,3])))/2;
                tensorsf = {cnew};
                legsf = {[-1,-2,1,1]};
                tmf = ncon(tensorsf,legsf);
                ctrf = eigs(tmf,1);
                ctrf = real(ctrf);
                cnew = d*cnew/ctrf;
                if bdl == 1
                    cnew = reshape(cnew,[d,d]);
                    cnewmd = mean(diag(cnew));
                    cnew = imag(cnew);
                    cnew = (cnew+rot90(cnew,2).')/2;
                    cnew = cnewmd*eye(d)+1i*cnew;
                    cnew = reshape(cnew,[bdl,bdl,d,d]);
                else
                    for nxf = 1:d
                        cnew(:,:,nxf,nxf) = zeros(bdl);
                        cnew(1,1,nxf,nxf) = 1;
                    end
                end
            end
            c = (cnew*sin(lambda*pi)-c*cos(lambda*pi));
            c = (c+conj(permute(c,[1,2,4,3])))/2;
            tensorsf = {c};
            legsf = {[-1,-2,1,1]};
            tmf = ncon(tensorsf,legsf);
            ctrf = eigs(tmf,1);
            ctrf = real(ctrf);
            c = d*c/ctrf;
            if bdl == 1
                c = reshape(c,[d,d]);
                cmd = mean(diag(c));
                c = imag(c);
                c = (c+rot90(c,2).')/2;
                c = cmd*eye(d)+1i*c;
                c = reshape(c,[bdl,bdl,d,d]);
            else
                for nxf = 1:d
                    c(:,:,nxf,nxf) = zeros(bdl);
                    c(1,1,nxf,nxf) = 1;
                end
            end
        end
        tensorsf = {c,b,a};
        legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdrhob*bdrhoa,bdl*bdrhob*bdrhoa]);
        l1_1val = eigs(tmf,1);
        tensorsf = {c,a,b};
        legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdrhoa*bdrhob,bdl*bdrhoa*bdrhob]);
        l1_2val = eigs(tmf,1);
        tensorsf = {c,a,c,a};
        legsf = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdrhoa*bdl*bdrhoa,bdl*bdrhoa*bdl*bdrhoa]);
        l2_1val = eigs(tmf,1);
        tensorsf = {c,c,a,a};
        legsf = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdl*bdrhoa*bdrhoa,bdl*bdl*bdrhoa*bdrhoa]);
        l2_2val = eigs(tmf,1);
        fom2val = abs(1+(-4*(l1_1val+l1_2val)+2*(l2_1val+l2_2val))/(4*l0val));
    else
        if lambda ~= 1
            if flag ~= 1
                tensorsf = {conj(c),b,a};
                legsf = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdrhob*bdrhoa,bdl*bdrhob*bdrhoa]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l1_1r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdrhob,bdrhoa]);
                l1_1l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdrhob,bdrhoa]);
                tensorsf = {conj(c),a,b};
                legsf = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdrhoa*bdrhob,bdl*bdrhoa*bdrhob]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l1_2r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdrhoa,bdrhob]);
                l1_2l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdrhoa,bdrhob]);
                tensorsf = {conj(c),a,c,a};
                legsf = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdrhoa*bdl*bdrhoa,bdl*bdrhoa*bdl*bdrhoa]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l2_1r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdrhoa,bdl,bdrhoa]);
                l2_1l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdrhoa,bdl,bdrhoa]);
                tensorsf = {conj(c),a,a,c};
                legsf = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdrhoa*bdrhoa*bdl,bdl*bdrhoa*bdrhoa*bdl]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l2_2r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdrhoa,bdrhoa,bdl]);
                l2_2l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdrhoa,bdrhoa,bdl]);
                tensorsf = {conj(c),c,a,a};
                legsf = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdl*bdrhoa*bdrhoa,bdl*bdl*bdrhoa*bdrhoa]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l2_3r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdl,bdrhoa,bdrhoa]);
                l2_3l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdl,bdrhoa,bdrhoa]);
                tensorsf = {l1_1l,b,a,l1_1r};
                legsf = {[-1,4,5],[4,2,-3,1],[5,3,1,-4],[-2,2,3]};
                l1_1 = ncon(tensorsf,legsf);
                l1_1 = reshape(l1_1,[bdl*bdl*d*d,1]);
                tensorsf = {l1_2l,a,b,l1_2r};
                legsf = {[-1,4,5],[4,2,-3,1],[5,3,1,-4],[-2,2,3]};
                l1_2 = ncon(tensorsf,legsf);
                l1_2 = reshape(l1_2,[bdl*bdl*d*d,1]);
                tensorsf = {l2_1l,a,a,l2_1r};
                legsf = {[-1,2,-5,4],[2,1,-3,-7],[4,3,-8,-4],[-2,1,-6,3]};
                l2_1 = ncon(tensorsf,legsf);
                l2_1 = reshape(l2_1,[bdl*bdl*d*d,bdl*bdl*d*d]);
                tensorsf = {l2_2l,a,a,eye(d),l2_2r};
                legsf = {[-1,4,5,-5],[4,2,-3,1],[5,3,1,-7],[-8,-4],[-2,2,3,-6]};
                l2_2 = ncon(tensorsf,legsf);
                l2_2 = reshape(l2_2,[bdl*bdl*d*d,bdl*bdl*d*d]);
                tensorsf = {l2_3l,eye(d),a,a,l2_3r};
                legsf = {[-1,-5,4,5],[-3,-7],[4,2,-8,1],[5,3,1,-4],[-2,-6,2,3]};
                l2_3 = ncon(tensorsf,legsf);
                l2_3 = reshape(l2_3,[bdl*bdl*d*d,bdl*bdl*d*d]);
                dl2 = 2*l2_1+l2_2+l2_3;
                dl2 = (dl2+dl2')/2;
                dl1 = 2*(l1_1+l1_2);
                dl2pinv = pinv2(dl2,tol2);
                dl2pinv = (dl2pinv+dl2pinv')/2;
                cv = dl2pinv*dl1;
                cnew = reshape(cv,[bdl,bdl,d,d]);
                tensorsf = {cnew};
                legsf = {[-1,-2,1,1]};
                tmf = ncon(tensorsf,legsf);
                ctrf = eigs(tmf,1);
                ctrf = real(ctrf);
                cnew = d*cnew/ctrf;
                if bdl == 1
                    cnew = reshape(cnew,[d,d]);
                    cnewmd = mean(diag(cnew));
                    cnew = imag(cnew);
                    cnew = (cnew+rot90(cnew,2).')/2;
                    cnew = cnewmd*eye(d)+1i*cnew;
                    cnew = reshape(cnew,[bdl,bdl,d,d]);
                else
                    for nxf = 1:d
                        cnew(:,:,nxf,nxf) = zeros(bdl);
                        cnew(1,1,nxf,nxf) = 1;
                    end
                end
            end
            c = (cnew*sin(lambda*pi)-c*cos(lambda*pi));
            tensorsf = {c};
            legsf = {[-1,-2,1,1]};
            tmf = ncon(tensorsf,legsf);
            ctrf = eigs(tmf,1);
            ctrf = real(ctrf);
            c = d*c/ctrf;
            if bdl == 1
                c = reshape(c,[d,d]);
                cmd = mean(diag(c));
                c = imag(c);
                c = (c+rot90(c,2).')/2;
                c = cmd*eye(d)+1i*c;
                c = reshape(c,[bdl,bdl,d,d]);
            else
                for nxf = 1:d
                    c(:,:,nxf,nxf) = zeros(bdl);
                    c(1,1,nxf,nxf) = 1;
                end
            end
        end
        tensorsf = {conj(c),b,a};
        legsf = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdrhob*bdrhoa,bdl*bdrhob*bdrhoa]);
        l1_1val = eigs(tmf,1);
        tensorsf = {conj(c),a,b};
        legsf = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdrhoa*bdrhob,bdl*bdrhoa*bdrhob]);
        l1_2val = eigs(tmf,1);
        tensorsf = {c,b,a};
        legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdrhob*bdrhoa,bdl*bdrhob*bdrhoa]);
        l1_3val = eigs(tmf,1);
        tensorsf = {c,a,b};
        legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdrhoa*bdrhob,bdl*bdrhoa*bdrhob]);
        l1_4val = eigs(tmf,1);
        tensorsf = {conj(c),a,c,a};
        legsf = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdrhoa*bdl*bdrhoa,bdl*bdrhoa*bdl*bdrhoa]);
        l2_1val = eigs(tmf,1);
        tensorsf = {conj(c),a,a,c};
        legsf = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdrhoa*bdrhoa*bdl,bdl*bdrhoa*bdrhoa*bdl]);
        l2_2val = eigs(tmf,1);
        tensorsf = {conj(c),c,a,a};
        legsf = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdl*bdrhoa*bdrhoa,bdl*bdl*bdrhoa*bdrhoa]);
        l2_3val = eigs(tmf,1);
        fom2val = abs(1+(-2*(l1_1val+l1_2val+l1_3val+l1_4val)+2*l2_1val+l2_2val+l2_3val)/(4*l0val));
    end
end
function [fom1,c] = fom1opt(c)
    fom1 = zeros([100,1]);
    fom1_1 = fom1f(1,c,0);
    [fom1_05,c_05] = fom1f(1/2,c,0);
    if fom1_05 > fom1_1
        c = c_05;
        fom1(1) = fom1_05;
    else
        fom1(1) = fom1_1;
    end
    fom1_1ph = fom1f(1+h1,c,0);
    if fom1_1ph > fom1(1)
        d1 = d01;
    else
        d1 = -d01;
    end
    testf = 0;
    flag = 0;
    cnew = [];
    iter1 = 1;
    while 1
        iter1 = iter1+1;
        [fom1val,c_1,c_2] = fom1f(1+d1,c,flag,cnew);
        if fom1val > fom1(iter1-1)
            c = c_1;
            fom1(iter1) = fom1val;
            testf = testf+1;
            flag = 0;
        else
            d1 = d1/2;
            fom1(iter1) = fom1(iter1-1);
            testf = 0;
            flag = 1;
            cnew = c_2;
        end
        if abs(d1) < h1 || (testf >= 4 && all(fom1(iter1-3:iter1) > 0) && std(fom1(iter1-3:iter1),1)/mean(fom1(iter1-3:iter1)) <= relunc1)
            break
        end
    end
    fom1 = fom1(1:iter1);
    %plot(fom1)
    fom1 = fom1(iter1);
end
function [fom2,c] = fom2opt(c)
    fom2 = zeros([100,1]);
    fom2_1 = fom2f(1,c,0);
    [fom2_05,c_05] = fom2f(1/2,c,0);
    if fom2_05 < fom2_1
        c = c_05;
        fom2(1) = fom2_05;
    else
        fom2(1) = fom2_1;
    end
    fom2_1ph = fom2f(1+h2,c,0);
    if fom2_1ph < fom2(1)
        d2 = d02;
    else
        d2 = -d02;
    end
    testf = 0;
    flag = 0;
    cnew = [];
    iter2 = 1;
    while 1
        iter2 = iter2+1;
        [fom2val,c_1,c_2] = fom2f(1+d2,c,flag,cnew);
        if fom2val < fom2(iter2-1)
            c = c_1;
            fom2(iter2) = fom2val;
            testf = testf+1;
            flag = 0;
        else
            d2 = d2/2;
            fom2(iter2) = fom2(iter2-1);
            testf = 0;
            flag = 1;
            cnew = c_2;
        end
        if abs(d2) < h2 || iter2 >= 40 || (testf >= 4 && std(fom2(iter2-3:iter2),1)/mean(fom2(iter2-3:iter2)) <= relunc2)
            break
        end
    end
    fom2 = fom2(1:iter2);
    %plot(-log10(fom2))
    fom2 = fom2(iter2);
end
% Calculations
%  MPO
bdrhoa = size(a,1);
bdrhob = size(b,1);
c = triu(ones(d)-eye(d));
c = 1i*phi*exp(-1)*c;
c = c+c';
c = eye(d)+c;
c = reshape(c,[bdl,bdl,d,d]);
%  [2Tr(rpL1)-Tr(rL1L1)-1]/phi^2 meanimization
if figureofmerit == 1
    [fom1,~] = fom1opt(c);
    f = fom1;
%  ||2*rp-rL1-L1r2||^2/||2*rp||^2 minimization
elseif figureofmerit == 2
    tensors = {a,a};
    legs = {[-1,-3,1,2],[-2,-4,2,1]};
    tm = ncon(tensors,legs);
    tm = reshape(tm,[bdrhoa*bdrhoa,bdrhoa*bdrhoa]);
    anorm = eigs(tm,1);
    anorm = abs(anorm)^(1/2);
    a = a/anorm;
    b = b/anorm;
    tensors = {b,b};
    legs = {[-1,-3,1,2],[-2,-4,2,1]};
    tm = ncon(tensors,legs);
    tm = reshape(tm,[bdrhoa*bdrhoa,bdrhoa*bdrhoa]);
    l0val = eigs(tm,1);
    l0val = abs(l0val);
    [fom2,c] = fom2opt(c);
    a = a*anorm;
    b = b*anorm;
    f = fom1f(1,c);
end
result = f;
fids = 1/4*result;
fid = 1-1/8*result*phi^2;
end