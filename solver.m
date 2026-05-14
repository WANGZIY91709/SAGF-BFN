function [Z,Wa,S,norm_error,objAll]=solver_xrv5(X,lambda1,lambda2,k)
num_view=size(X,2);
num_X=size(X{1,1},2);
% error X-XZ=E，<Yv,X-XZ-E>
Yv=cell(1,num_view);
E=cell(1,num_view);
% XTZ,Z=Q,<Bv,Z-Q>
Bv=cell(1,num_view);
Q=cell(1,num_view);
% <W,Z-J>,Z=J
Wv=cell(1,num_view);
J=cell(1,num_view);

for v=1:num_view
    % Lagrange multiplier
    Yv{v} = zeros(num_X,num_X);
    Bv{v} = zeros(num_X,num_X);
    Wv{v} = zeros(num_X,num_X);
    E{v} = zeros(num_X,num_X);
    Q{v} = zeros(num_X,num_X);
    J{v} = zeros(num_X,num_X);
end
Isconverg = zeros(1,num_view);
epson = 1e-5;
iter = 0;
rho=0.001;max_rho = 10e12; pho_rho =1.3;
%  Z
Z=cell(1,num_view);
distX=cell(1,num_view);
distX1=cell(1,num_view);
idx=cell(1,num_view);
for v=1:num_view
    distX{v} = L2_distance_1(X{v},X{v});
    [distX1{v}, idx{v}] = sort(distX{v},2);
    Z{1,v} = zeros(num_X);
    for i = 1:num_X
        di = distX1{v}(i,2:k+2);
        id = idx{v}(i,2:k+2);
        Z{1,v}(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
    end
    Z{1,v}=(Z{1,v}+Z{1,v}');
end
%  S
S=zeros(num_X,num_X);
%  Wa
Wa=ones(num_view,1)/num_view;
%% 
while(sum(Isconverg) == 0)
    %% A
    temp_Z=Z;
    temp_Yv=Yv;
    temp_E=E;
    temp_Bv=Bv;
    temp_Q=Q;
    temp_Wv=Wv;
    temp_J=J;
    %% 
    for v=1:num_view
        E_B = X{v} - temp_E{v}+ temp_Yv{v}/rho;
        Q_B = pinv(X{v})*(temp_Q{v}-temp_Bv{v}/rho+temp_J{v}-temp_Wv{v}/rho+Wa(v)*S);
        for j = 1:num_X
            cols = [1:j-1, j+1:num_X];
            X_tilde = X{v}(:, cols); 
            eb_j = E_B(:, j);
            qb = Q_B(:,j);
            XtX = X_tilde' * X_tilde + (2+Wa(v))*eye(num_X-1); 
            temp_1 = pinv(XtX);
            z_j = temp_1* (X_tilde' * (eb_j+qb));
            Z{v}(cols, j) = z_j;
            Z{v}(j, j) = 0;
        end
        Z{v} = normalize_columns(Z{v});
    end
    %% E(v)
    for v=1:num_view
        Dv=X{v}-X{v}*temp_Z{v}+temp_Yv{v}/rho;
        E{1,v}=zeros(num_X);
        for i=1:num_X
            threshold = 1 / rho; 
            norm_dvi = norm(Dv(:,i), 2);
            if norm_dvi > threshold
                E{v}(:,i) = (1 - threshold / norm_dvi) * Dv(:,i);
            else
                E{v}(:,i) = 0;
            end
        end
    end
    %% Q(v)
    beta1=cell(1,num_view);
    F1=cell(1,num_view);
    M1=cell(1,num_view);
    for v=1:num_view
        distx = L2_distance_1(X{v}',X{v}');
        distx=distx-diag(diag(distx));
        beta1{v}=distx;
        F1{v}=temp_Z{v}+temp_Bv{v}/rho;
        M1{v}=F1{v}-beta1{v}/rho;
        Q{1,v} = zeros(num_X);
        for i=1:num_X
            threshold = lambda1/rho;
            norm_mi=norm(M1{v}(:,i),2);
            if norm_mi > threshold
                Q{v}(:,i) = (1 - threshold / norm_mi) * M1{v}(:,i);
            else
                Q{v}(:,i) = 0;
            end
        end
    end
    %% J
    Z_tensor = cat(3, temp_Z{:,:});
    W_tensor = cat(3, temp_Wv{:,:});
    temp_Zt = Z_tensor(:);
    temp_Wt = W_tensor(:);
    sX = [num_X, num_X, num_view];
    %twist-version
    [j, objV] = Gshrink(temp_Zt + 1/rho*temp_Wt,(num_X*lambda2)/rho,sX,0,3);
    J_tensor = reshape(j, sX);
    for v=1:num_view
        J{v} = J_tensor(:,:,v);
    end
    %% B，
    if iter>0
        for v=1:num_view
            v_fnorm=norm(S-Z{v},'fro');
            Wa(v)=1/(v_fnorm+eps);
        end
        Wa=Wa/sum(Wa);
    end
    %% C，
    Z_ALL=cell(1,num_view);
    Z_AVE=zeros(num_X);
    for v=1:num_view
        Z_ALL{v}=Wa(v)*Z{v};
        Z_AVE=Z_AVE+Z_ALL{v};
    end
    S = zeros(num_X);
    for i=1:num_X
        temp_i=1:num_X;
        S(i,temp_i~=i) = update_L2(Z_AVE(i,temp_i~=i));
    end
    %% 
    Z_tensor_upd=cat(3, Z{:,:});
    temp_Zt_upd = Z_tensor_upd(:);
    temp_Wt = temp_Wt + rho*(temp_Zt_upd - j);
    W_tensor = reshape(temp_Wt , sX);
    for v=1:num_view
        % Lagrange multiplier
        Yv{v} = Yv{v}+rho*(X{v}-X{v}*Z{v}-E{v});
        Bv{v} = Bv{v}+rho*(Z{v}-Q{v});
        Wv{v} = W_tensor(:,:,v);
    end
    %% 
    norm_error.XXZ_E(iter+1)=0;
    norm_error.Z_Q(iter+1)=0;
    norm_error.Z_J(iter+1)=0;
    objAll(iter+1)=0;
    Isconverg = ones(1,num_view);
    for v=1:num_view
        % XXZ_E
        norm_XXZ_E=0;
        norm_XXZ_E=abs(norm(X{v}-X{v}*Z{v}-E{v},'fro'));
        norm_error.XXZ_E(iter+1) = norm_error.XXZ_E(iter+1)+norm_XXZ_E;
        % Z_Q
        norm_Z_Q=0;
        norm_Z_Q=abs(norm(Z{v}-Q{v},'fro'));
        norm_error.Z_Q(iter+1) = norm_error.Z_Q(iter+1)+norm_Z_Q;
        % Z_J
        norm_Z_J=0;
        norm_Z_J=lambda2*abs(norm(Z{v}-J{v},'fro'));
        norm_error.Z_J(iter+1) = norm_error.Z_J(iter+1)+norm_Z_J;
        if (norm_XXZ_E>epson||norm_Z_Q>epson||norm_Z_J>epson)
            Isconverg(v) = 0;
        end
        Isconverg(v) = 0;
    end
    if (iter>30)
        Isconverg  = ones(1,num_view);
    end
    iter = iter + 1;
    rho = min(rho*pho_rho, max_rho);
end
end