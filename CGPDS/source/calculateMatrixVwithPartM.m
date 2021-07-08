function tempMatrix = calculateMatrixVwithPartM( model,samd,samr )
% calculate matrixs for the computation of elbo and gradient
    if model.Rtr == 1
        lend = length(samd);
        lenr = length(samr);
        
        tempMatrix.TrYYdr = trace(model.partm*model.partm');

        tempMatrix.C_v = model.invLm_v * model.Psi4 * model.invLmT_v;
        if isempty(find(isnan(tempMatrix.C_v)))&&isempty(find(isinf(tempMatrix.C_v)))
        else
            n = size(tempMatrix.C_v,1);
            X = diag(rand(n,1));
            U = orth(rand(n,n));
            tempMatrix.C_v = U'*X*U;
        end
        
        tempMatrix.TrC_v = sum(diag(tempMatrix.C_v)); % Tr(C)
        % Matrix At replaces the matrix A of the old implementation; At is more stable
        % since it has a much smaller condition number than A=sigma^2 K_uu + Psi2
        tempMatrix.At_v = (1/model.beta) * eye(size(tempMatrix.C_v,1)) + tempMatrix.C_v; % At = beta^{-1} I + C
        tempMatrix.Lat_v = jitChol(tempMatrix.At_v)';
        tempMatrix.invLat_v = tempMatrix.Lat_v\eye(size(tempMatrix.Lat_v,1));  
        tempMatrix.invLatT_v = tempMatrix.invLat_v';
        tempMatrix.logDetAt_v = 2*(sum(log(diag(tempMatrix.Lat_v)))); % log |At|

        tempMatrix.P1_v = tempMatrix.invLat_v * model.invLm_v; % M x M

        % First multiply the two last factors; so, the large N is only involved
        % once in the calculations (P1: MxM, Psi1':MxN, Y: NxD)
        tempMatrix.P_v = tempMatrix.P1_v * (model.Psi0' * model.partm);%M*D

        tempMatrix.P1TP1_v = (tempMatrix.P1_v' * tempMatrix.P1_v);
        % Needed for both, the bound's and the derivs. calculations.
        tempMatrix.TrPP_v = sum(sum(tempMatrix.P_v .* tempMatrix.P_v));%1*1


        %%% Precomputations for the derivatives (of the likelihood term) of the bound %%%

        %model.B = model.invLmT * model.invLatT * model.P; %next line is better
        tempMatrix.B_v = tempMatrix.P1_v' * tempMatrix.P_v;%M*D

        Tb_v = (1/model.beta) * lend* (tempMatrix.P1_v' * tempMatrix.P1_v);
            Tb_v = Tb_v + (tempMatrix.B_v * tempMatrix.B_v');
        tempMatrix.T1_v = lend* model.invK_vv - Tb_v;
        
    elseif model.Rtr>1
        
    end
   

end

