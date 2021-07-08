function p=calculatePosteriorLikelihoodCMGPDS(model,Testmeans,Testcovars,Yts,index,scale,origBias,origScale)
if ~exist('scale','var')
    scale = ones(1,length(index));
end
Yts(:,index) = Yts(:,index).*repmat(scale,[size(Yts,1),1]);
jitter = 1e-5;
switch model.type
    case 'varmgplvm'
        convolveKernTmp=model.kern.comp{1,1};
        Ainv=model.P1'*model.P1;
        alpha=Ainv*model.Psi1'*model.m;
        Z=model.X_u;
        beta=model.beta;
        D=convolveKernTmp.outputDimension;
        Kuu=model.K_uu;
        
        Nf=50;
        Ny=size(Testmeans,1);
        
        % Ny=10;
        for i=1:Ny
            Xmu=Testmeans(i,:)';
            Xvar=diag(Testcovars(i,:));
            
            Xtests=gsamp(Xmu,Xvar,Nf);
%             Xtests = Xmu;
            for k=1:Nf % number of sample F 
                Xtest=Xtests(k,:);
                Kfu=convolveKernComputeKfu(convolveKernTmp,Xtest,Z);
                Kff=convolveKernComputeKff(convolveKernTmp,Xtest);
                Fmu=Kfu*alpha;
                Fvar=abs(Kfu*Ainv*Kfu'/beta+Kff-Kfu*Kuu^-1*Kfu');
%                 Fvar=Fvar.*(model.scale'*model.scale);
                if det(Fvar)==0
%                 Fvar=diag(diag(Fvar));
                    Fvar = Fvar + jitter*eye(D);
                end
                if det(Fvar)<0
                    Fvar=diag(diag(Fvar));
                   
                end
%                 Fvar=diag(diag(Fvar));
                
                Ftest=gsamp(Fmu',Fvar,1);
                Ymu=Ftest;
                Yvar=eye(D)/beta+diag(diag(Fvar));
                % Rescale the mean
                Ymu = Ymu.*model.scale;

                % Add the bias back in
                Ymu = Ymu + model.bias;
                Yvar(index,index) = Yvar(index,index).*diag(model.scale(index).^2);
                %for angle space or scaled space
                Ymu(index) = Ymu(index).*scale;
                
                Yvar(index,index) = Yvar(index,index).*diag(scale).^2;
                pdf4y(i,k)=mvnpdf(Yts(i,index),Ymu(index),Yvar(index,index));
%                 pdf4y(i,k)=-D/2*log(2*pi)-0.5*log(det(Yvar(index,index)))-0.5*(Yts(i,index)-Ymu(index))*Yvar(index,index)^-1*(Yts(i,index)-Ymu(index))';
%                 diff(i,k) = sum((Yts(i,index)-Ymu(index)).^2);
            end
        end  
        p0=mean(pdf4y,2);   
       
    case 'vargplvm'
        kernTmp_v=model.kern_v.comp{1,1};
        Ainv_v=model.P1_v'*model.P1_v;
        alpha_v=Ainv_v*model.Psi0'*model.m;%M*D
        
        Kinvk_u = cell(model.d,model.J);
        vars_u = zeros(model.d,model.J);
        
        kernTmp_u=model.kern_u.comp{1,1};
        Ainv_u = cell(model.d,model.J);
        alpha_u = zeros(model.k,model.d,model.J);%M*D*J
        for j = 1:model.J
            for d = 1:model.d
                Ainv_u{d,j} = model.P1_u{d,j}' * model.P1_u{d,j};
                alpha_u(:,d,j) = Ainv_u{d,j}*model.Psi1'*model.m(:,d); % size: M*1
                Kinvk_u{d,j} = (model.invK_uu - (1/model.beta)*Ainv_u{d,j});
                vars_u(d,j) = Psi3_star - sum(sum(Kinvk_u{d,j}.*Psi5_star));
            end
        end
        
        
        Z=model.X_u;
        Zh=model.X_v;
        
        beta=model.beta;
        D=model.d;
        Kuu=model.K_uu;
        Kvv=model.K_vv;
        W = model.W;
        
        Nf=50;
        Ny=size(Testmeans,1);  
        vard = vardistCreate(zeros(1,model.q), model.q, 'gaussian');%1*9 varidist
        Kinvk_v = (model.invK_vv - (1/model.beta)*Ainv_v);%model.invK_vv - inv(model.K_vv+model.beta*model.Psi4)
       
        
%         varsigmah = zeros(size(vardistX.means,1),model.d);
%         varsigmag = zeros(size(vardistX.means,1),model.d,model.J);
%         varsigma = zeros(size(vardistX.means,1),model.d);
        
        for i=1:Ny
            
            varsigmah = zeros(1,model.d);
            varsigmag = zeros(1,model.d,model.J);
            varsigma = zeros(1,model.d);

            Xmu=Testmeans(i,:)';
            Xvar=diag(Testcovars(i,:));
            Xtests=gsamp(Xmu,Xvar,Nf);
            
            vard.means = Testmeans(i,:);
            vard.covars = Testcovars(i,:);
                  % compute psi0 term
                  Psi0_star = kernVardistPsi1Compute(model.kern_v, vard, Zh);%1*M
                  Psi1_star = kernVardistPsi1Compute(model.kern_u, vard, Z);%1*M
                  
                  Psi2_star = kernVardistPsi0Compute(model.kern_v, vard);%1*1
                  Psi3_star = kernVardistPsi0Compute(model.kern_u, vard);%1*1
                  % compute psi2 term
                  Psi4_star = kernVardistPsi2Compute(model.kern_v, vard, Zh);%M*M
                  Psi5_star = kernVardistPsi2Compute(model.kern_u, vard, Z);%M*M

                  vars_v = Psi2_star - sum(sum(Kinvk_v.*Psi4_star));%M*M

            for k=1:Nf % number of sample F
                Xtest=Xtests(k,:);
%                 Kgu=rbfard2KernCompute(kernTmp_u,Xtest,Z);%1*M
%                 Kgg=rbfard2KernCompute(kernTmp_u,Xtest);%1*1
%                 
%                 Khv=rbfard2KernCompute(kernTmp_v,Xtest,Z);%1*M
%                 Khh=rbfard2KernCompute(kernTmp_u,Xtest);%1*1

                
                Hmu = Psi0_star*alpha_v;%1*D

                Gmu = zeros(1,D,model.J);
                
                Gvar = zeros(D,D,model.J);
                Gtest = zeros(1,D,model.J);
                
                for d = 1:model.d
                    varsigmah(1,d) = alpha_v(:,d)'*(Psi4_star*alpha_v(:,d)) - Hmu(1,d)^2;%1*1
                    for j = 1:model.J
                        Gmu(1,d,j) = Psi1_star*alpha_u(1,d,j);%1*1
                        varsigmag(1,d,j) = alpha_u(:,d,j)'*(Psi5_star*alpha_u(:,d,j)) - Gmu(1,d,j)^2;                        
                    end
                end

                 
                
%                 tempVarg = sum(varsigmag,3);
                varsigmah(1,:) = varsigmah(1,:) + repmat(vars_v,1,model.d);
%                 varsigma(n,:) = varsigma(n,:) + varsigmah(n,:) + repmat(vars_v,1,model.d) + tempVarg(n,:) + sum(vars_u,2)'; 


%                 Hvar=eye(D)*abs(Khv*Ainv_v*Khv'/beta+Khh-Khv*Kvv^-1*Khv');%D*D
                Hvar = diag(varsigmah(1,:));
                Htest=gsamp(Hmu',Hvar,1);%1*D
                
                Ymu=Htest;%1*D
                
                
                for j = 1:model.J
                    varsigmag(1,:,j) = varsigmag(1,:,j) + vars_u(:,j);
                    Gvar(:,:,j) = diag(varsigmag(1,:,j));
                    Gtest(:,:,j) = gsamp(Gmu(:,:,j)',Gvar(:,:,j),1);%1*D
                    Ymu = Ymu + W(:,j)'.*Gtest(:,:,j);%1*D
                end
                
               
               %
               
                  %
                  
                
                Ymu = Ymu./repmat(origScale, size(Ymu, 1), 1);
                Ymu = Ymu + repmat(origBias, size(Ymu, 1), 1);
                   %for angle space or scaled space
                Yvar(index,index) = Yvar(index,index).*diag(model.scale(index).^2);
                Ymu(index) = Ymu(index).*scale;
                
%                 Yvar(index,index) = Yvar(index,index).*diag(scale).^2;
%                 pdf4y(i,k)=mvnpdf(Yts(i,index),Ymu(index),Yvar(index,index));
                pdf4y(i,k)=-D/2*log(2*pi)-0.5*log(det(Yvar(index,index)))-0.5*(Yts(i,index)-Ymu(index))*Yvar(index,index)^-1*(Yts(i,index)-Ymu(index))';
                diff(i,k) = sum((Yts(i,index)-Ymu(index)).^2);
            end
        end
        p0=mean(pdf4y,2);   
end
        
%         p=mean(log(p0));
        p=mean(p0);
        diff = mean(mean(diff));
end
