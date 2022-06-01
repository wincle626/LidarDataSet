function preprocess() 
    %% Data import: CBCS for real LiDAR data 
    % Outline:
    %           1. Select and load appropriate LiDAR scene (5 available)
    %                   2 options: 25% pattern density and 50% pattern density at 64 patterns total                
    %           2. Reformat data into A, YQ, YI 
    % To adjust sample numbers for processing
    % e.g. m<n with m=8
    %           A(:,1:8,:), YQ(:,1:8), YI(:,1:8)
    %           samples can be chosen randomly and/or from different start points
    %           e.g XI = admm_lasso(A, YI, lambda, rho, alpha); for m=64 > n=6
    %               -> XI = admm_lasso(A(:,4:12,:), YI(:,4:12), lambda, rho, alpha); for m=8 < n=16
    %%
    close all
    clear all
    clc
    %%
    ratio = '50perc'; % options: '25perc' and '50perc' - ratio of non-zeros/zeros in pattern
    scene_select = 3;
    %
    switch scene_select
        case 1
            load(['neopropeneBCS_64_' ratio '.mat']);
        case 2
            load(['concblocksBCS_64_' ratio '.mat']);
        case 3
            load(['floatingBCS_64_' ratio '.mat']);
        case 4
            load(['sandblocksBCS_64_' ratio '.mat']);
        case 5
            load(['pipesBCS_64_' ratio '.mat']);
    end
    
    
    %% pinv data
    A = cbcs_data.A; YQ = cbcs_data.YQ; YI = cbcs_data.YI;
    samplingnum = 2304;
    sparsenum = 48;
    photonnum = 16;
    para.ROW = photonnum;
    para.COL = sparsenum;
    para.DIM = samplingnum;
    para.RHO = 1.2;
    para.ALPHA = 1.4;
    para.A_ROW = photonnum;
    para.A_COL = sparsenum;
    para.Y_ROW = samplingnum;
    para.Y_COL = sparsenum;
    para.X_ROW = samplingnum;
    para.X_COL = photonnum;
    para.dim_cb = [4,4];
    para.dim = [192,192];
    para.climit = 300;
    para.dmax = 39;
    para.dmin = 27;
    % algorithm parameters
    para.rho = 1.2;
    para.alpha = 1.4;
    para.MAX_ITER = 5; % gradient descent maximum iteration number
    rndind1 = 1:sparsenum;
    rndind2 = 1:samplingnum;
    dctm = dctmtx(photonnum);
    dctminv = dctm';
    A = A(rndind2,rndind1,:);
    yq = YQ(rndind2,rndind1);
    yi = YI(rndind2,rndind1);
    At = zeros(samplingnum,photonnum,sparsenum);
    L = zeros(samplingnum,sparsenum,sparsenum);
    U = zeros(samplingnum,sparsenum,sparsenum);
    Linv = zeros(samplingnum,sparsenum,sparsenum);
    Uinv = zeros(samplingnum,sparsenum,sparsenum);
    lambda_yq_dct = zeros(1,samplingnum);
    lambda_yi_dct = zeros(1,samplingnum);
    for i=1:samplingnum
        A(i,:,:) = reshape(A(i,:,:),sparsenum,photonnum);
        At(i,:,:) = reshape(A(i,:,:),sparsenum,photonnum).';
        [l,u] = lu(reshape(A(i,:,:),sparsenum,photonnum)*reshape(A(i,:,:),sparsenum,photonnum).');
        L(i,:,:) = l; U(i,:,:) = u;
        Linv(i,:,:) = pinv(l); Uinv(i,:,:) = pinv(u);
        lambda_max = norm( reshape(A(i,:,:),sparsenum,photonnum).'*yq(i,:).', 'inf' );
        lambda_yq_dct(i) = 0.1*lambda_max;
        lambda_max = norm( reshape(A(i,:,:),sparsenum,photonnum).'*yi(i,:).', 'inf' );
        lambda_yi_dct(i) = 0.1*lambda_max;
    end        
    
    %% test pinv
    Ap = zeros(samplingnum,photonnum,sparsenum);
    for i=1:samplingnum
        tmp = reshape(A(i,:,:),sparsenum,photonnum);
        Ap(i,:,:) = pinv(tmp'*tmp)*tmp';
    end
    [xd,xq,xi] = pseudoinv1(para, Ap, yi, yq);  
    tempx = 1;
    tempy = 1;
    for i=1:samplingnum
        XD_temp = reshape(xd(i, 1:para.dim_cb(1)*para.dim_cb(2)), [para.dim_cb(1) para.dim_cb(2)]);
        XD_reshaped(tempx:tempx+para.dim_cb(1)-1, tempy:tempy+para.dim_cb(2)-1) = XD_temp;
        tempx = tempx+para.dim_cb(1);
        if tempx == para.dim(1)+1
            tempx = 1;
            tempy = tempy + para.dim_cb(2);
        end
    end
    XD = XD_reshaped;
    figure
    imagesc(XD);
    title('Depth (PINV)')
    title(colorbar,'Distance, cm')
    caxis manual
    caxis([para.dmin para.dmax])
    pbaspect([para.dim(1)/para.dim(2) 1 1])
    
    %% admm data
    A = cbcs_data.A; YQ = cbcs_data.YQ; YI = cbcs_data.YI;
    samplingnum = 2304;
    sparsenum = 8;
    photonnum = 16;
    para.ROW = photonnum;
    para.COL = sparsenum;
    para.DIM = samplingnum;
    para.RHO = 1.2;
    para.ALPHA = 1.4;
    para.A_ROW = photonnum;
    para.A_COL = sparsenum;
    para.Y_ROW = samplingnum;
    para.Y_COL = sparsenum;
    para.X_ROW = samplingnum;
    para.X_COL = photonnum;
    para.dim_cb = [4,4];
    para.dim = [192,192];
    para.climit = 300;
    para.dmax = 39;
    para.dmin = 27;
    % algorithm parameters
    para.rho = 1.2;
    para.alpha = 1.4;
    para.MAX_ITER = 5; % gradient descent maximum iteration number
    rndind1 = 1:sparsenum;
    rndind2 = 1:samplingnum;
    dctm = dctmtx(photonnum);
    dctminv = dctm';
    A = A(rndind2,rndind1,:);
    yq = YQ(rndind2,rndind1);
    yi = YI(rndind2,rndind1);
    At = zeros(samplingnum,photonnum,sparsenum);
    L = zeros(samplingnum,sparsenum,sparsenum);
    U = zeros(samplingnum,sparsenum,sparsenum);
    Linv = zeros(samplingnum,sparsenum,sparsenum);
    Uinv = zeros(samplingnum,sparsenum,sparsenum);
    lambda_yq_dct = zeros(1,samplingnum);
    lambda_yi_dct = zeros(1,samplingnum);
    for i=1:samplingnum
        A(i,:,:) = reshape(A(i,:,:),sparsenum,photonnum)*dctm;
        At(i,:,:) = reshape(A(i,:,:),sparsenum,photonnum).';
        [l,u] = lu(reshape(A(i,:,:),sparsenum,photonnum)*reshape(A(i,:,:),sparsenum,photonnum).');
        L(i,:,:) = l; U(i,:,:) = u;
        Linv(i,:,:) = pinv(l); Uinv(i,:,:) = pinv(u);
        lambda_max = norm( reshape(A(i,:,:),sparsenum,photonnum).'*yq(i,:).', 'inf' );
        lambda_yq_dct(i) = 0.1*lambda_max;
        lambda_max = norm( reshape(A(i,:,:),sparsenum,photonnum).'*yi(i,:).', 'inf' );
        lambda_yi_dct(i) = 0.1*lambda_max;
    end 
    
    %% test admm
    [xd,xq,xi] = admmcs1(para, A, At, yi, yq, lambda_yi_dct, lambda_yq_dct, ...
                    dctminv, Linv, Uinv);                
    tempx = 1;
    tempy = 1;
    for i=1:samplingnum
        XD_temp = reshape(xd(i, 1:para.dim_cb(1)*para.dim_cb(2)), [para.dim_cb(1) para.dim_cb(2)]);
        XD_reshaped(tempx:tempx+para.dim_cb(1)-1, tempy:tempy+para.dim_cb(2)-1) = XD_temp;
        tempx = tempx+para.dim_cb(1);
        if tempx == para.dim(1)+1
            tempx = 1;
            tempy = tempy + para.dim_cb(2);
        end
    end
    XD = XD_reshaped;
    figure
    imagesc(XD);
    title('Depth (ADMM)')
    title(colorbar,'Distance, cm')
    caxis manual
    caxis([para.dmin para.dmax])
    pbaspect([para.dim(1)/para.dim(2) 1 1])
    
    %% test pgd
    [xd,xq,xi] = pgdcs1(para, A, yi, yq, dctminv); 
    tempx = 1;
    tempy = 1;
    for i=1:samplingnum
        XD_temp = reshape(xd(i, 1:para.dim_cb(1)*para.dim_cb(2)), [para.dim_cb(1) para.dim_cb(2)]);
        XD_reshaped(tempx:tempx+para.dim_cb(1)-1, tempy:tempy+para.dim_cb(2)-1) = XD_temp;
        tempx = tempx+para.dim_cb(1);
        if tempx == para.dim(1)+1
            tempx = 1;
            tempy = tempy + para.dim_cb(2);
        end
    end
    XD = XD_reshaped;
    figure
    imagesc(XD);
    title('Depth (PGD)')
    title(colorbar,'Distance, cm')
    caxis manual
    caxis([para.dmin para.dmax])
    pbaspect([para.dim(1)/para.dim(2) 1 1])
    
end


function [xd,xq,xi] = pseudoinv1(para, A,yi,yq)
    xd = zeros(para.X_ROW,para.X_COL);
    xq = zeros(para.X_ROW,para.X_COL);
    xi = zeros(para.X_ROW,para.X_COL);
    for cbiter = 1 : para.X_ROW
        for  row = 1 : para.A_ROW
            xi = 0;
            xq = 0;
            for col = 1 : para.A_COL
                xi = xi + A(cbiter,row,col)*yi(cbiter,col);
                xq = xq + A(cbiter,row,col)*yq(cbiter,col);
            end
            if xi~=0
                temp = xq/xi;
                temp = min(temp, 300);
                temp = max(temp, 0);
                xd(cbiter,row) = temp;
            else
                xd(cbiter,row) = 0;
            end
        end
    end
end


function [xd,xq,xi] = admmcs1(para, A, At, yi, yq, ...
                             lambda_yi_dct, lambda_yq_dct, ...
                             dctinv, Linv, Uinv)  
    % para.rho = 1.2; para.alpha = 1.4; % EUSIPCO 2020
    xd = zeros(para.X_ROW,para.X_COL);
    for i=1:para.X_ROW
        [xi, ~] = admm(para, reshape(A(i,:,:),8,16), reshape(At(i,:,:),16,8), yi(i,:)', lambda_yi_dct(i), ...
                              reshape(Linv(i,:,:),8,8), reshape(Uinv(i,:,:),8,8), para.rho, para.alpha, 0);
        [xq, ~] = admm(para, reshape(A(i,:,:),8,16), reshape(At(i,:,:),16,8), yq(i,:)', lambda_yq_dct(i),...
                              reshape(Linv(i,:,:),8,8), reshape(Uinv(i,:,:),8,8), para.rho, para.alpha, 0);
        xi(xi==0) = 1;
        xq(xi==0) = 0;
        xd_dct = (dctinv * xq) ./ (dctinv * xi);
%         xd_dct = xq ./ xi;
    %     xd_dct = xq ./ xi;
    %     xd_dct(xd_dct>300) = 300;
    %     xd_dct(xd_dct<0) = 0;
        xd(i,:) = xd_dct';
    end
end


function [xd,xq,xi] = pgdcs1(para, A, yi, yq, dctinv)
    for i=1:2304
        A2 = reshape(A(i,:,:),8,16);
        AtA = A2' * A2;
        L = max(eig(AtA));
        lambad = 1.8/L; % EUSIPCO 2020 fp64
        Atyi = A2' * yi(i,:)';
        Atyq = A2' * yq(i,:)';
    %     gamma_max = norm(Atyi,'inf');
    %     gamma = 0.1*gamma_max;
    %     xi = pgd(AtA, Atyi, 1, 0.01, gamma);
        xi = pgd(para, AtA, Atyi, lambad, 1, 0); % EUSIPCO 2020
    %     gamma_max = norm(Atyq,'inf');
    %     gamma = 0.1*gamma_max;
    %     xq = pgd(AtA, Atyq, 1, 0.01, gamma);
        xq = pgd(para, AtA, Atyq, lambad, 1, 0); % EUSIPCO 2020
        xi(xi==0) = 1;
        xq(xi==0) = 0;
        xd_dct = (dctinv * xq) ./ (dctinv * xi);
%         xd_dct = xq ./ xi;
        xd_dct(xd_dct>300) = 300;
        xd_dct(xd_dct<0) = 0;
        xd(i,:) = xd_dct';
    end

end