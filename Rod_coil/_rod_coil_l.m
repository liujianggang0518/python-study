close all
clear all
clc

Beta = 2;
Eita = 4
chiN = 8;
f = 0.5;
Lambda = 30;
D = 4.0;
Nx = 64;
z = linspace(0, D-D/Nx, Nx)';
dt = 0.02;
s = 0:dt:1;
Ns = length(s);
Nf_1 = find(s == f);
Nf_2 = find(s == 1 - f);
K = 2 * pi *[0:Nx/2-1, -Nx/2:-1]/D;
Ntheta = 8;
Nphi = 17;
dphi = 2 * pi / Nphi;
phi = 0: dphi: 2 * pi - dphi;
[costheta, wtheta]=lgwt(Ntheta,-1,1);
theta = acos(costheta);
mu_positive = -2 * cos(2 * pi * z / D);
mu_negative = 2 * cos(2 * pi * z / D);
M_z = [-0.4, 0, 0;
       0, -0.4, 0;
       0, 0, 0.8];
M = zeros(3, 3, Nx);

for i=1:Nx
    M(:,:,i) = M_z
end
%准备工作
coeff_1 = exp(-K'.^2 * dt / 2);
W_M = zeros(3, 3, Ntheta, Nphi);
coeff_3 = zeros(Nx, Ntheta);
coeff_4 = zeros(Nx, Ntheta);
for j=1:Ntheta
    coeff_3(:,j) = exp(-Beta * cos(theta(j)) * K' * 1i * dt / 2);
    coeff_4(:,j) = exp(Beta * cos(theta(j)) * K' * 1i * dt / 2);
    for n=1:Nphi
        W_M(:,:,j,n) = [cos(phi(n)) * sin(theta(j)), sin(phi(n)) * sin(theta(j)), cos(theta(j))]' * [cos(phi(n)) * sin(theta(j)), sin(phi(n)) * sin(theta(j)), cos(theta(j))] - 1 / 3 * eye(3,3);
    end
end
l = 0:Ntheta;
Nl = length(l);
YPlm = prepare(Nl, Ntheta, Nphi, Ntheta);
l_matrix = zeros(Ntheta+1, Ntheta+1);
for i=1:Ntheta+1
    l_matrix(i, i) = exp(-l(i) * (l(i)+1) / 2 / Lambda * dt);
end
step = 0;
while(1)
    q_A = ones(Nx, Nf_1);
    q_B = zeros(Ntheta, Nphi, Nx, Ns - Nf_1 + 1);
    q_A_plus = zeros(Nx, Ns - Nf_2 + 1);
    q_B_plus = zeros(Ntheta, Nphi, Nx, Nf_2);
    w_a = mu_positive - mu_negative;
    w_b = mu_positive + mu_negative;
    coeff_2 = exp(-w_a .* dt);
    %q_A求解：
    for i=1:Nf_1 -1
        q_A(:, i+1) = fft(q_A(:, i)).*coeff_1;
        q_A(:, i+1) = ifft(q_A(:, i+1)).*coeff_2;
        q_A(:, i+1) = ifft(fft(q_A(:, i+1)).*coeff_1);
    end
    %q_B求解:
    %赋初值
    for r=1:Nx
        q_B(:,:,r,1) = q_A(r, Nf_1) / 4 / pi;
    end
    %第一个算子半步:固定r,u
    for i=1:Ns - Nf_1
        for r=1:Nx
            for j=1:Ntheta
                for n=1:Nphi
                    q_B(j, n, r, i+1) = q_B(j, n, r, i) * exp((-w_b(r) + sum(sum(M(:,:,r).* W_M(:,:,j,n))))*dt/2); 
                end
            end
        end
        %第二个算子半步：固定theta,phi
        for n=1:Nphi
            for j=1:Ntheta
                q_B(j, n, :, i+1) = reshape(fft(q_B(j, n, :, i+1)), Nx, 1).*coeff_3(:, j);
                q_B(j, n, Nx/2+1, i+1)=real(q_B(j, n, Nx/2+1, i+1));
                q_B(j, n, :, i+1) = ifft(q_B(j, n, :, i+1));
            end
        end
        %第三个算子利用球谐函数展开求解
        for r=1:Nx
            q_B(:, :, r, i+1) = spherical2real_c2r(Ntheta + 1, Nphi, Ntheta, Ntheta, l_matrix*real2spherical_r2c(Ntheta + 1, Nphi, Ntheta, Ntheta, q_B(:,:, r, i+1),YPlm),YPlm);
        end
        %第二个算子一步：固定theta,phi
        for n=1:Nphi
            for j=1:Ntheta
                q_B(j, n, :, i+1) = reshape(fft(q_B(j, n, :, i+1)), Nx, 1).*coeff_3(:,j);
                q_B(j, n, Nx/2+1, i+1)=real(q_B(j, n, Nx/2+1, i+1));
                q_B(j, n, :, i+1) = ifft(q_B(j, n, :, i+1));
            end
        end
       
       %第一个算子一步:固定r,u
       for j=1:Ntheta
           for n=1:Nphi
               for r=1:Nx
                   q_B(j, n, r, i+1) = q_B(j, n, r, i) * exp((-w_b(r) + sum(sum(M(:,:,r).* W_M(:,:,j,n))))*dt/2);
               end
           end
       end
    end    
    %q_B_plus求解：
    %赋初值
    q_B_plus(:, :, :, 1) = 1 / 4 / pi;
    %第一个算子半步:固定r,u
    for i=1:Nf_2 - 1
        for j=1:Ntheta
            for n=1:Nphi
                for r=1:Nx
                    q_B_plus(j, n, r, i+1) =  q_B_plus(j, n, r, i) * exp((-w_b(r) + sum(sum(M(:,:,r).* W_M(:,:,j,n))))*dt/2);
                end
            end
        end
        %第二个算子半步：固定theta,phi
        for n=1:Nphi
            for j=1:Ntheta
                q_B_plus(j, n, :, i+1) = reshape(fft(q_B_plus(j, n, :, i+1)), Nx, 1).*coeff_4(:,j);
                q_B_plus(j, n, Nx/2+1, i+1)=real(q_B_plus(j, n, Nx/2+1, i+1));
                q_B_plus(j, n, :, i+1) = ifft(q_B_plus(j, n, :, i+1));
            end
        end
        %第三个算子利用球谐函数展开求解
        for r=1:Nx
           q_B_plus(:, :, r, i+1) = spherical2real_c2r(Ntheta + 1, Nphi, Ntheta, Ntheta, l_matrix*real2spherical_r2c(Ntheta + 1, Nphi, Ntheta, Ntheta, q_B_plus(:,:, r, i+1),YPlm),YPlm); 
        end
        %第二个算子一步：固定theta,phi
        for n=1:Nphi
            for j=1:Ntheta
                q_B_plus(j, n, :, i+1) = reshape(fft(q_B_plus(j, n, :, i+1)), Nx, 1).*coeff_4(:,j);
                q_B_plus(j, n, Nx/2+1, i+1)=real(q_B_plus(j, n, Nx/2+1, i+1));
                q_B_plus(j, n, :, i+1) = ifft(q_B_plus(j, n, :, i+1));
            end
        end
        %第一个算子一步:固定r,u

        for j=1:Ntheta
            for n=1:Nphi
                for r=1:Nx
                    q_B_plus(j, n, r, i+1) =  q_B_plus(j, n, r, i) * exp((-w_b(r) + sum(sum(M(:,:,r).* W_M(:,:,j,n))))*dt/2);
                end
            end
        end
    end
    %q_A_plus求解：
    for r=1:Nx
        q_A_plus(r, 1)=sum(sum(q_B_plus(:,:,r,Nf_2),2)*dphi.*wtheta)
    end
    for i=1:Ns - Nf_2
        q_A_plus(:, i+1) = fft(q_A_plus(:, i)).*coeff_1;
        q_A_plus(:, i+1) = ifft(q_A_plus(:, i+1)).*coeff_2;
        q_A_plus(:, i+1) = ifft(fft(q_A_plus(:, i+1)).*coeff_1);
    end
    %计算Q
    QQ = zeros(Nx,1);
    q_B_plus_tixing1 = q_B_plus(:,:,r,Nf_2).*4*pi.*q_B(:,:,:,1);
    q_B_sita1=permute(sum(q_B_plus_tixing1, 2),[1,3,2])*dphi;
    for r=1:Nx
        QQ(r)=sum(q_B_sita1(:,r).*wtheta);
    end
    QQ = real(fft(QQ)/Nx);
    Q = QQ(1)
    fprintf('%.6e \n', Q);
    
    %对q_A_plus,q_B_plus进行翻转:
    q_A_plus = q_A_plus(:,end:-1:1);
    q_B_plus = q_B_plus(:,:,:,end:-1:1);
    %求解faiA,faiB:
    %faiA
    qqplus_A = q_A .* q_A_plus;
    sum_A = 0;
    for i=1:Nf_1
        sum_A = sum_A + qqplus_A(:, i);
    end
    faiA = 1 / Q * (-5 / 8 * (qqplus_A(:,1) + qqplus_A(:,Nf_1)) + 1 / 6 * (qqplus_A(:,2) + qqplus_A(:,Nf_1 -1)) - 1 / 24 * (qqplus_A(:,3) + qqplus_A(:,Nf_1 - 2)) + sum_A) * dt;
    qqplus_B = q_B .* q_B_plus;
    qq_B = zeros(Nx, Nf_2);
    for i=1:Nf_2
        for r=1:Nx
            qq_B(r, i) = sum(sum(qqplus_B(:,:,r,i),2) * dphi .* wtheta);
        end
    end
    sum_B = 0;
    for n=1:Nf_2
        sum_B = sum_B + qq_B(:, n);
    end
    faiB = 4 * pi / Q * (-5 / 8 * (qq_B(:,1) + qq_B(:,Nf_2)) + 1 / 6 * (qq_B(:,2) + qq_B(:,Nf_2 -1)) - 1 / 24 * (qq_B(:,3) + qq_B(:,Nf_2 - 2)) + sum_B) * dt;
    %求解S(r):
    q_u_qplus = zeros(3, 3, Ntheta, Nphi, Nx, Nf_2);
    for i=1:Nf_2
        for r=1:Nx
            for j=1:Ntheta
                for n=1:Nphi
                    q_u_qplus(:,:,j,n,r,i) = q_B(j,n,r,i) * W_M(:,:,j,n) * q_B_plus(j,n,r,i);
                end
            end
        end
    end
    qq_s1 = permute(sum(q_u_qplus,4) * dphi, [1, 2, 3, 5, 6, 4]);
    sum_s1 = zeros(3, 3, Ntheta, Nx, Nf_1);
    for i=1:Nf_2
        for r=1:Nx
            for j=1:Ntheta
                sum_s1(:,:,j,r,i) = qq_s1(:,:,j,r,i) * wtheta(j);
            end
        end
    end
    qq_s2 = permute(sum(sum_s1,3), [1, 2, 4, 5, 3]);
    sum_s2 = 0
    for n=1:Nf_2
        sum_s2 = sum_s2 + qq_s2(:,:,:,n)
    end
    S = 4 * pi / Q * (-5 / 8 * (qq_s2(:,:,:,1) + qq_s2(:,:,:,Nf_2)) + 1 / 6 * (qq_s2(:,:,:,2) + qq_s2(:,:,:,Nf_2 -1)) - 1 / 24 * (qq_s2(:,:,:,3) + qq_s2(:,:,:,Nf_2 - 2)) + sum_s2) * dt; 
    %迭代条件
    update_S = zeros(3,3);
    SS = S - M / Eita / chiN;
    for i=1:3
        for j=1:3
            update_S(i,j) = norm(permute(SS(i,j,:),[3, 1, 2]))^2;
        end
    end
    step = step + 1;
    E_M = sqrt(sum(update_S(:)) / 9);
    fprintf('迭代步数:%d \t 误差1：%.15f \t 误差2：%.15f \t 误差3：%.15f \n', step, norm(faiA + faiB - 1), norm( faiA - faiB - 2 / chiN * mu_negative ), E_M);
    if (norm(faiA + faiB - 1)>10^-5 || norm(faiA - faiB - 2 / chiN * mu_negative)>10^-5 ||E_M>10^-5)
        mu_poitive = mu_positive + (faiA + faiB - 1) * 1;
        mu_negative =  mu_negative + 1 * (faiA - faiB - 2 * mu_negative./ chiN);
        M = M + 1*(S - M / Eita / chiN);
        
        fprintf('\n ');
        
    else
        break;
    end
    
end

    
    
    
    
