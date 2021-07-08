function p=calculatePosteriorLikelihood(model,Testmeans,Testcovars,Yts,index,scale)
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
        kernTmp=model.kern.comp{1,1};
        Ainv=model.P1'*model.P1;
        alpha=Ainv*model.Psi1'*model.m;
        Z=model.X_u;
        beta=model.beta;
        D=model.d;
        Kuu=model.K_uu;
        Nf=50;
        Ny=size(Testmeans,1);    
        for i=1:Ny
            Xmu=Testmeans(i,:)';
            Xvar=diag(Testcovars(i,:));
            Xtests=gsamp(Xmu,Xvar,Nf);
            
            for k=1:Nf % number of sample F
                Xtest=Xtests(k,:);
                Kfu=rbfard2KernCompute(kernTmp,Xtest,Z);
                Kff=rbfard2KernCompute(kernTmp,Xtest);

                Fmu=Kfu*alpha;
                Fvar=eye(D)*abs(Kfu*Ainv*Kfu'/beta+Kff-Kfu*Kuu^-1*Kfu');
%                 Fvar=Fvar.*(model.scale'*model.scale);
                Ftest=gsamp(Fmu',Fvar,1);
                
                Ymu=Ftest;
                Yvar=eye(D)/beta+diag(diag(Fvar));
                % Rescale the mean
                Ymu = Ymu.*model.scale;

                % Add the bias back in
                Ymu = Ymu + model.bias;
                   %for angle space or scaled space
                Yvar(index,index) = Yvar(index,index).*diag(model.scale(index).^2);
                Ymu(index) = Ymu(index).*scale;
                
                Yvar(index,index) = Yvar(index,index).*diag(scale).^2;
                pdf4y(i,k)=mvnpdf(Yts(i,index),Ymu(index),Yvar(index,index));
%                 pdf4y(i,k)=-D/2*log(2*pi)-0.5*log(det(Yvar(index,index)))-0.5*(Yts(i,index)-Ymu(index))*Yvar(index,index)^-1*(Yts(i,index)-Ymu(index))';
%                 diff(i,k) = sum((Yts(i,index)-Ymu(index)).^2);
            end
        end
        p0=mean(pdf4y,2);   
end
        
        p=mean(log(p0));
%         diff = mean(mean(diff));
end
