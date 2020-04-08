function [result,a0,c] = MPO_A_PBC_Inf(figureofmerit,d,bdpsi,bdl,ti,noiserange,integral0,integral1,integral2,lherm,phi,imprecision,a0,c)
% Parameters
tol1 = imprecision*phi^2;
tol2 = imprecision*phi^2;
told = imprecision*phi^2;
d01 = 10^-1;
d02 = 10^-1;
d0d = 10^-1;
h1 = 10^-10;
h2 = 10^-10;
hd = 10^-10;
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
function [fom1val,c,cnew] = fom1f(lambda,c,flag,cnew)
    if lambda ~= 1
        if flag ~= 1
            tensorsf = {c,b};
            legsf = {[-1,-3,1,2],[-2,-4,2,1]};
            tmf = ncon(tensorsf,legsf);
            tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise]);
            [tmvrf,~] = eigs(tmf,1);
            [tmvlf,~] = eigs(tmf.',1);
            tmvnorm = tmvlf.'*tmvrf;
            l1r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise]);
            l1l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise]);
            tensorsf = {c,a,c};
            legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
            tmf = ncon(tensorsf,legsf);
            tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdl,bdl*bdpsi^2*bdnoise*bdl]);
            [tmvrf,~] = eigs(tmf,1);
            [tmvlf,~] = eigs(tmf.',1);
            tmvnorm = tmvlf.'*tmvrf;
            l2r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdl]);
            l2l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdl]);
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
    tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise]);
    l1val = eigs(tmf,1);
    tensorsf = {c,a,c};
    legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
    tmf = ncon(tensorsf,legsf);
    tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdl,bdl*bdpsi^2*bdnoise*bdl]);
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
                tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l1_1r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise]);
                l1_1l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise]);
                tensorsf = {c,a,b};
                legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l1_2r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise]);
                l1_2l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise]);
                tensorsf = {c,a,c,a};
                legsf = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdl*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdl*bdpsi^2*bdnoise]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l2_1r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise]);
                l2_1l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise]);
                tensorsf = {c,c,a,a};
                legsf = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l2_2r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise]);
                l2_2l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise]);
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
        tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
        l1_1val = eigs(tmf,1);
        tensorsf = {c,a,b};
        legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
        l1_2val = eigs(tmf,1);
        tensorsf = {c,a,c,a};
        legsf = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdl*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdl*bdpsi^2*bdnoise]);
        l2_1val = eigs(tmf,1);
        tensorsf = {c,c,a,a};
        legsf = {[-1,-5,1,2],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
        l2_2val = eigs(tmf,1);
        fom2val = abs(1+(-4*(l1_1val+l1_2val)+2*(l2_1val+l2_2val))/(4*l0val));
    else
        if lambda ~= 1
            if flag ~= 1
                tensorsf = {conj(c),b,a};
                legsf = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l1_1r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise]);
                l1_1l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise]);
                tensorsf = {conj(c),a,b};
                legsf = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l1_2r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise]);
                l1_2l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise]);
                tensorsf = {conj(c),a,c,a};
                legsf = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdl*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdl*bdpsi^2*bdnoise]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l2_1r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise]);
                l2_1l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdl,bdpsi^2*bdnoise]);
                tensorsf = {conj(c),a,a,c};
                legsf = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise*bdl,bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise*bdl]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l2_2r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise,bdl]);
                l2_2l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise,bdl]);
                tensorsf = {conj(c),c,a,a};
                legsf = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
                tmf = ncon(tensorsf,legsf);
                tmf = reshape(tmf,[bdl*bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
                [tmvrf,~] = eigs(tmf,1);
                [tmvlf,~] = eigs(tmf.',1);
                tmvnorm = tmvlf.'*tmvrf;
                l2_3r = reshape(tmvrf/tmvnorm^(1/2),[bdl,bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise]);
                l2_3l = reshape(tmvlf.'/tmvnorm^(1/2),[bdl,bdl,bdpsi^2*bdnoise,bdpsi^2*bdnoise]);
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
        tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
        l1_1val = eigs(tmf,1);
        tensorsf = {conj(c),a,b};
        legsf = {[-1,-4,2,1],[-2,-5,2,3],[-3,-6,3,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
        l1_2val = eigs(tmf,1);
        tensorsf = {c,b,a};
        legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
        l1_3val = eigs(tmf,1);
        tensorsf = {c,a,b};
        legsf = {[-1,-4,1,2],[-2,-5,2,3],[-3,-6,3,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
        l1_4val = eigs(tmf,1);
        tensorsf = {conj(c),a,c,a};
        legsf = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdl*bdpsi^2*bdnoise,bdl*bdpsi^2*bdnoise*bdl*bdpsi^2*bdnoise]);
        l2_1val = eigs(tmf,1);
        tensorsf = {conj(c),a,a,c};
        legsf = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise*bdl,bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise*bdl]);
        l2_2val = eigs(tmf,1);
        tensorsf = {conj(c),c,a,a};
        legsf = {[-1,-5,2,1],[-2,-6,2,3],[-3,-7,3,4],[-4,-8,4,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdl*bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdl*bdl*bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
        l2_3val = eigs(tmf,1);
        fom2val = abs(1+(-2*(l1_1val+l1_2val+l1_3val+l1_4val)+2*l2_1val+l2_2val+l2_3val)/(4*l0val));
    end
end
function [fomdval,a0,a0new] = fomdf(lambda,a0,flag,a0new)
    if lambda ~= 1
        if flag ~= 1
            tensorsf = {conj(a0),cpd,a0};
            legsf = {[-1,-4,1],[-2,-5,1,2],[-3,-6,2]};
            tmf = ncon(tensorsf,legsf);
            tmf = reshape(tmf,[bdpsi*bdherm*bdl*bdnoise*bdpsi,bdpsi*bdherm*bdl*bdnoise*bdpsi]);
            [tmvrf,~] = eigs(tmf,1);
            [tmvlf,~] = eigs(tmf.',1);
            tmvnorm = tmvlf.'*tmvrf;
            lpdr = reshape(tmvrf/tmvnorm^(1/2),[bdpsi,bdherm*bdl*bdnoise,bdpsi]);
            lpdl = reshape(tmvlf.'/tmvnorm^(1/2),[bdpsi,bdherm*bdl*bdnoise,bdpsi]);
            tensorsf = {conj(a0),c2d,a0};
            legsf = {[-1,-4,1],[-2,-5,1,2],[-3,-6,2]};
            tmf = ncon(tensorsf,legsf);
            tmf = reshape(tmf,[bdpsi*(bdherm*bdl)^2*bdnoise*bdpsi,bdpsi*(bdherm*bdl)^2*bdnoise*bdpsi]);
            [tmvrf,~] = eigs(tmf,1);
            [tmvlf,~] = eigs(tmf.',1);
            tmvnorm = tmvlf.'*tmvrf;
            l2dr = reshape(tmvrf/tmvnorm^(1/2),[bdpsi,(bdherm*bdl)^2*bdnoise,bdpsi]);
            l2dl = reshape(tmvlf.'/tmvnorm^(1/2),[bdpsi,(bdherm*bdl)^2*bdnoise,bdpsi]);
            tensorsf = {conj(a0),a0};
            legsf = {[-1,-3,1],[-2,-4,1]};
            tmf = ncon(tensorsf,legsf);
            tmf = reshape(tmf,[bdpsi*bdpsi,bdpsi*bdpsi]);
            [tmvrf,~] = eigs(tmf,1);
            [tmvlf,~] = eigs(tmf.',1);
            tmvnorm = tmvlf.'*tmvrf;
            psinormr = reshape(tmvrf/tmvnorm^(1/2),[bdpsi,bdpsi]);
            psinorml = reshape(tmvlf.'/tmvnorm^(1/2),[bdpsi,bdpsi]);
            tensorsf = {lpdl,cpd,lpdr};
            legsf = {[-1,2,-4],[2,1,-3,-6],[-2,1,-5]};
            lpd = ncon(tensorsf,legsf);
            lpd = reshape(lpd,[bdpsi*bdpsi*d,bdpsi*bdpsi*d]);
            tensorsf = {l2dl,c2d,l2dr};
            legsf = {[-1,2,-4],[2,1,-3,-6],[-2,1,-5]};
            l2d = ncon(tensorsf,legsf);
            l2d = reshape(l2d,[bdpsi*bdpsi*d,bdpsi*bdpsi*d]);
            psinormrpinv = pinv2(psinormr,told);
            psinormrpinv = (psinormrpinv+psinormrpinv')/2;
            psinormlpinv = pinv2(psinorml,told);
            psinormlpinv = (psinormlpinv+psinormlpinv')/2;
            tensorsf = {psinormlpinv,eye(d),psinormrpinv};
            legsf = {[-1,-4],[-3,-6],[-2,-5]};
            psinormpinv = ncon(tensorsf,legsf);
            psinormpinv = reshape(psinormpinv,[bdpsi*bdpsi*d,bdpsi*bdpsi*d]);
            eiginput = 2*lpd-l2d;
            eiginput = (eiginput+eiginput')/2;
            eiginput = psinormpinv*eiginput;
            [a0v,~] = eigs(eiginput,1,'lr');
            a0new = reshape(a0v,[bdpsi,bdpsi,d]);
            a0new = (a0new+conj(permute(a0new,[2,1,3])))/2;
            a0new = (a0new+permute(flip(a0new,3),[2,1,3]))/2;
            a0new = (a0new+permute(rot90(a0new,2),[2,1,3]))/2;
            tensorsf = {conj(a0new),a0new};
            legsf = {[-1,-3,1],[-2,-4,1]};
            tmf = ncon(tensorsf,legsf);
            tmf = reshape(tmf,[bdpsi*bdpsi,bdpsi*bdpsi]);
            a0normf = eigs(tmf,1);
            a0normf = abs(a0normf)^(1/2);
            a0new = a0new/a0normf;
        end
        a0 = (a0new*sin(lambda*pi)-a0*cos(lambda*pi));
        a0 = (a0+conj(permute(a0,[2,1,3])))/2;
        a0 = (a0+permute(flip(a0,3),[2,1,3]))/2;
        a0 = (a0+permute(rot90(a0,2),[2,1,3]))/2;
        tensorsf = {conj(a0),a0};
        legsf = {[-1,-3,1],[-2,-4,1]};
        tmf = ncon(tensorsf,legsf);
        tmf = reshape(tmf,[bdpsi*bdpsi,bdpsi*bdpsi]);
        a0normf = eigs(tmf,1);
        a0normf = abs(a0normf)^(1/2);
        a0 = a0/a0normf;
    end
    tensorsf = {conj(a0),cpd,a0};
    legsf = {[-1,-4,1],[-2,-5,1,2],[-3,-6,2]};
    tmf = ncon(tensorsf,legsf);
    tmf = reshape(tmf,[bdpsi*bdherm*bdl*bdnoise*bdpsi,bdpsi*bdherm*bdl*bdnoise*bdpsi]);
    lpdval = eigs(tmf,1);
    tensorsf = {conj(a0),c2d,a0};
    legsf = {[-1,-4,1],[-2,-5,1,2],[-3,-6,2]};
    tmf = ncon(tensorsf,legsf);
    tmf = reshape(tmf,[bdpsi*(bdherm*bdl)^2*bdnoise*bdpsi,bdpsi*(bdherm*bdl)^2*bdnoise*bdpsi]);
    l2dval = eigs(tmf,1);
    fomdval = real((2*lpdval-l2dval-1)/phi^2);
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
function [fomd,a0] = fomdopt(a0)
    fomd = zeros([100,1]);
    fomd_1 = fomdf(1,a0,0);
    [fomd_05,a0_05] = fomdf(1/2,a0,0);
    if fomd_05 > fomd_1
        a0 = a0_05;
        fomd(1) = fomd_05;
    else
        fomd(1) = fomd_1;
    end
    fomd_1ph = fomdf(1+hd,a0,0);
    if fomd_1ph > fomd(1)
        dd = d0d;
    else
        dd = -d0d;
    end
    testf = 0;
    flag = 0;
    a0new = [];
    iterd = 1;
    while 1
        iterd = iterd+1;
        [fomdval,a0_1,a0_2] = fomdf(1+dd,a0,flag,a0new);
        if fomdval > fomd(iterd-1)
            a0 = a0_1;
            fomd(iterd) = fomdval;
            testf = testf+1;
            flag = 0;
        else
            dd = dd/2;
            fomd(iterd) = fomd(iterd-1);
            testf = 0;
            flag = 1;
            a0new = a0_2;
        end
        if abs(dd) < hd || (testf >= 4 && all(fomd(iterd-3:iterd) > 0) && std(fomd(iterd-3:iterd),1)/mean(fomd(iterd-3:iterd)) <= reluncd)
            break
        end
    end
    fomd = fomd(1:iterd);
    %plot(fomd)
    fomd = fomd(iterd);
end
% Calculations
%  MPO
if noiserange == 0
    bdnoise = 1;
	an = ones([1,1,d,d]);
    integralp = 0;
elseif noiserange == 1
    nxpm = repmat(1:d,[d,1]);
    nxm = nxpm.';
    nxminusnxp = nxm-nxpm;
    bdnoise = 1;
    an = exp(-nxminusnxp.^2*integral0/2);
    an = permute(an,[3,4,1,2]);
    integralp = integral0;
elseif noiserange == 2
    nxpm = repmat(1:d,[d,1]);
    nxm = nxpm.';
    nxminusnxp = nxm-nxpm;
    bdnoise1 = 1;
    an1 = exp(-nxminusnxp.^2*integral0/2);
    an1 = permute(an1,[3,4,1,2]);
    bdnoise2 = 2*d-1;
    corr12 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral1);
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
    integralp = integral0+2*integral1;
elseif noiserange == 3
    nxpm = repmat(1:d,[d,1]);
    nxm = nxpm.';
    nxminusnxp = nxm-nxpm;
    bdnoise1 = 1;
    an1 = exp(-nxminusnxp.^2*integral0/2);
    an1 = permute(an1,[3,4,1,2]);
    bdnoise2 = 2*d-1;
    corr12 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral1);
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
    bdnoise3ini = 2*d-1;
    bdnoise3 = bdnoise3ini^2;
    corr13 = exp(-nxminusnxp(:)*nxminusnxp(:).'*integral2);
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
    bdnoise = bdnoise1*bdnoise2*bdnoise3;
    an = zeros([bdnoise,bdnoise,d,d]);
    for nx = 1:d
        for nxp = 1:d
            an(:,:,nx,nxp) = an1(:,:,nx,nxp)*kron(an2(:,:,nx,nxp),an3(:,:,nx,nxp));
        end
    end
    integralp = integral0+2*integral1+2*integral2;
end
f_iter = zeros([2*100,1]);
iter = 0;
while 1
    iter = iter+1;
    a = zeros([bdpsi^2*bdnoise,bdpsi^2*bdnoise,d,d]);
    b = zeros([bdpsi^2*bdnoise,bdpsi^2*bdnoise,d,d]);
    for nx = 1:d
        for nxp = 1:d
            a(:,:,nx,nxp) = kron(kron(a0(:,:,nx),conj(a0(:,:,nxp))),an(:,:,nx,nxp));
            b(:,:,nx,nxp) = (1-(1i/sqrt(ti))*(nx-nxp)*integralp*phi)*a(:,:,nx,nxp);
        end
    end
%  [2Tr(rpL1)-Tr(rL1L1)-1]/phi^2 meanimization
    if figureofmerit == 1
        [fom1,c] = fom1opt(c);
        f = fom1;
%  ||2*rp-rL1-L1r2||^2/||2*rp||^2 minimization
    elseif figureofmerit == 2
        tensors = {a,a};
        legs = {[-1,-3,1,2],[-2,-4,2,1]};
        tm = ncon(tensors,legs);
        tm = reshape(tm,[bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
        anorm = eigs(tm,1);
        anorm = abs(anorm)^(1/2);
        a = a/anorm;
        b = b/anorm;
        tensors = {b,b};
        legs = {[-1,-3,1,2],[-2,-4,2,1]};
        tm = ncon(tensors,legs);
        tm = reshape(tm,[bdpsi^2*bdnoise*bdpsi^2*bdnoise,bdpsi^2*bdnoise*bdpsi^2*bdnoise]);
        l0val = eigs(tm,1);
        l0val = abs(l0val);
        [fom2,c] = fom2opt(c);
        a = a*anorm;
        b = b*anorm;
        f = fom1f(1,c);
    end
    f_iter(2*iter-1) = f;
	if (iter >= 5 && std(f_iter(2*iter-8:2*iter-1),1)/mean(f_iter(2*iter-8:2*iter-1)) <= relunc) || iter >= 100
        break
	end
    if figureofmerit == 2 && iter >= 40
        break
    end
%  MPO (dual)
    if lherm == 1
        bdherm = 1;
        ch = c;
    else
        bdherm = 2;
        ch = zeros([bdherm*bdl,bdherm*bdl,d,d]);
        for nx = 1:d
            for nxp = 1:d
                ch(:,:,nx,nxp) = [c(:,:,nx,nxp),zeros(bdl);zeros(bdl),conj(c(:,:,nxp,nx))];
            end
        end
        ch = ch/2;
    end
    c2 = zeros([(bdherm*bdl)^2,(bdherm*bdl)^2,d,d]);
    for nx = 1:d
        for nxp = 1:d
            for nxpp = 1:d
                c2(:,:,nx,nxp) = c2(:,:,nx,nxp)+kron(ch(:,:,nx,nxpp),ch(:,:,nxpp,nxp));
            end
        end
    end
    cd = zeros([bdherm*bdl*bdnoise,bdherm*bdl*bdnoise,d,d]);
    cpd = zeros([bdherm*bdl*bdnoise,bdherm*bdl*bdnoise,d,d]);
    c2d = zeros([(bdherm*bdl)^2*bdnoise,(bdherm*bdl)^2*bdnoise,d,d]);
    for nx = 1:d
        for nxp = 1:d
            cd(:,:,nx,nxp) = kron(an(:,:,nx,nxp),ch(:,:,nx,nxp));
            cpd(:,:,nx,nxp) = (1+1i*(nx-nxp)*phi)*cd(:,:,nx,nxp);
            c2d(:,:,nx,nxp) = kron(an(:,:,nx,nxp),c2(:,:,nx,nxp));
        end
    end
%  Tr[r0*X] meanimization
    [fomd,a0] = fomdopt(a0);
    f_iter(2*iter) = fomd;
end
f_iter = f_iter(1:2*iter-1);
%plot(f_iter)
result = f_iter(2*iter-1);
end