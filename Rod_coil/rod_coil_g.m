clc;clear;
%参数的定义
beta=2; eita=4; kN=8; lamda=30; D=4;
%把时间进行剖分:Ns份,步长为dt:0.02
tao=0.02;s=0:tao:1;Ns=length(s);
%表示A部分所占的体积分数
f=0.5;Nf=find(s==f);%正向
N1_f=find(s==1-f);%逆向
%把空间[0,D]进行剖分，剖分Nx份.
Nx=64;x=linspace(0,D-D/Nx,Nx)';
k=2*pi*[0:Nx/2-1 -Nx/2:-1]/D;
%方向u(sita,phi)
Nphi=17;dt1=2*pi/Nphi;
phi=0:dt1:2*pi-dt1;%phi方向是均分的
Nsita=8;
[cosita,wsita]=lgwt(Nsita,-1,1);%sita方向是不等分，调用gauss函数包[cossita,权重]
sita=acos(cosita);
%初始场函数
u_plus = -2*cos(2*pi*x/D);  %u+
u_dis = 2*cos(2*pi*x/D);  %u-
%M是与空间剖分x有关的
matrix=[-0.4,0,0;
         0,-0.4,0;
         0,0,0.8;];
% matrix=[-2  0  0 ;
%          0  3  2;
%          0  2 -1;];
 M=zeros(3,3,Nx);
 for i=1:Nx
%  M(:,:,i)=2*matrix*(u_dis(i)+2);
   M(:,:,i)=matrix;
 end
%准备工作
coff1=exp(-k'.^2*tao/2);
jzb=zeros(3,3,Nsita,Nphi);
coff3=zeros(Nx,Nsita);
coff4=zeros(Nx,Nsita);
for k1=1:Nsita %sita
    coff3(:,k1)=exp(-beta*cos(sita(k1))*k'*1i*tao/2);
    coff4(:,k1)=exp(beta*cos(sita(k1))*k'*1i*tao/2);
    for n=1:Nphi  %phi
        jzb(:,:,k1,n)=[cos(phi(n))*sin(sita(k1)),sin(phi(n))*sin(sita(k1)),cos(sita(k1))]'*[cos(phi(n))*sin(sita(k1)),sin(phi(n))*sin(sita(k1)),cos(sita(k1))]-1/3*eye(3,3);
    end
end
l=0:Nsita;Nl=length(l);
YPlm=prepare(Nl,Nsita,Nphi,Nsita);
l_matrix=zeros(Nsita+1,Nsita+1);
for i=1:Nsita+1
l_matrix(i,i)=exp(-l(i)*(l(i)+1)/2/lamda*tao);
end
step=0;
 while(1)
 qA=ones(Nx,Nf);%正向传播qA
 qA_plus=zeros(Nx,Ns-N1_f+1);%逆向传播qA+
 qB=zeros(Nsita,Nphi,Nx,Ns-Nf+1);%正向传播qB,它与(sita,phi,r,s)有关
 qB_plus=zeros(Nsita,Nphi,Nx,N1_f);%逆向传播qA+
 wa=u_plus-u_dis;
 wb=u_plus+u_dis;
 coff2=exp(-wa.*tao);
 %正向传播子，从柔性链qA，到刚性链qB
 %算子qA
 for i=1:Nf-1
     qA(:,i+1)=fft(qA(:,i)).*coff1;
     qA(:,i+1)=ifft(qA(:,i+1)).*coff2;
     qA(:,i+1)=ifft(fft(qA(:,i+1)).*coff1);
 end
 %给qB赋初值
 for r=1:Nx
 qB(:,:,r,1)=qA(r,Nf)/4/pi;
 end
 for i=1:Ns-Nf %时间s,这是一个大循环:
% 固定r,u;第一阶段,实空间
      for r=1:Nx %空间r   
           for k1=1:Nsita %sita
              for n=1:Nphi  %phi
                  qB(k1,n,r,i+1)=qB(k1,n,r,i)*exp((-wb(r)+sum(sum(M(:,:,r).*jzb(:,:,k1,n))))*tao/2);                           
              end
           end
       end
% % 第二阶段复空间,固定u
    for n=1:Nphi  %phi
        for k1=1:Nsita %sita
         qB(k1,n,:,i+1)=reshape(fft(qB(k1,n,:,i+1)),Nx,1).*coff3(:,k1);
         qB(k1,n,Nx/2+1,i+1)=real(qB(k1,n,Nx/2+1,i+1));
          qB(k1,n,:,i+1)=ifft(qB(k1,n,:,i+1));
        end
    end
%     qB(1,1,:,2)
% % 第三阶段求调和
%      YPlm=prepare(Nl,Nsita,Nphi,Nsita); 
     for r=1:Nx
         qB(:,:,r,i+1)=spherical2real_c2r(Nsita+1,Nphi,Nsita,Nsita, l_matrix*real2spherical_r2c(Nsita+1,Nphi,Nsita,Nsita,qB(:,:,r,i+1),YPlm),YPlm);
     end

% %第四阶段复空间,固定u     
    for  n=1:Nphi  %phi
         for k1=1:Nsita %sita
         qB(k1,n,:,i+1)=reshape(fft(qB(k1,n,:,i+1)),Nx,1).*coff3(:,k1);
         qB(k1,n,Nx/2+1,i+1)=real(qB(k1,n,Nx/2+1,i+1));
          qB(k1,n,:,i+1)=ifft(qB(k1,n,:,i+1)); 
         end
     end
% %第五阶段实空间,固定r,u;   
       for k1=1:Nsita %sita
           for n=1:Nphi  %phi
               for r=1:Nx %空间r
                   qB(k1,n,r,i+1)=qB(k1,n,r,i+1)*exp((-wb(r)+sum(sum(M(:,:,r).*jzb(:,:,k1,n))))*tao/2);                           
               end
           end
       end  
 end
     %%%%%%%%%%%%%%%%%%%%%逆向传播子%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %给qB赋初值
 qB_plus(:,:,:,1)=1/4/pi;
 for i=1:N1_f-1 %时间s,这是一个大循环:
%固定r,u;第一阶段,实空间
  for k1=1:Nsita %sita
       for n=1:Nphi  %phi 
             for r=1:Nx %空间r
                qB_plus(k1,n,r,i+1)=qB_plus(k1,n,r,i)*exp((-wb(r)+sum(sum(M(:,:,r).*jzb(:,:,k1,n))))*tao/2);                           
             end
        end
   end
%第二阶段复空间,固定u
     for n=1:Nphi  %phi
         for k1=1:Nsita %sita
         qB_plus(k1,n,:,i+1)=reshape(fft(qB_plus(k1,n,:,i+1)),Nx,1).*coff4(:,k1);
         qB_plus(k1,n,Nx/2+1,i+1)=real(qB_plus(k1,n,Nx/2+1,i+1));
          qB_plus(k1,n,:,i+1)=ifft(qB_plus(k1,n,:,i+1));
         end
     end
%第三阶段求调和
%      YPlm=prepare(Nl,Nsita,Nphi,Nsita); 
     for r=1:Nx
        qB_plus(:,:,r,i+1)=spherical2real_c2r(Nsita+1,Nphi,Nsita,Nsita, l_matrix*real2spherical_r2c(Nsita+1,Nphi,Nsita,Nsita,qB_plus(:,:,r,i+1),YPlm),YPlm);
     end
%第四阶段复空间,固定u     
     for n=1:Nphi  %phi
         for k1=1:Nsita %sita
         qB_plus(k1,n,:,i+1)=reshape(fft(qB_plus(k1,n,:,i+1)),Nx,1).*coff4(:,k1);
         qB_plus(k1,n,Nx/2+1,i+1)=real(qB_plus(k1,n,Nx/2+1,i+1));
          qB_plus(k1,n,:,i+1)=ifft(qB_plus(k1,n,:,i+1));
         end
     end
%第五阶段实空间,固定r,u;   
     for k1=1:Nsita %sita
         for n=1:Nphi  %phi
             for r=1:Nx %空间r
                 qB_plus(k1,n,r,i+1)=qB_plus(k1,n,r,i+1)*exp((-wb(r)+sum(sum(M(:,:,r).*jzb(:,:,k1,n))))*tao/2);                           
             end
         end
     end 
 end
%算子qA_plus
for r=1:Nx
    qA_plus(r,1)=sum(sum(qB_plus(:,:,r,N1_f),2)*dt1.*wsita);
end
for i=1:Ns-N1_f
     qA_plus(:,i+1)=fft(qA_plus(:,i)).*coff1;
     qA_plus(:,i+1)=ifft(qA_plus(:,i+1)).*coff2;
     qA_plus(:,i+1)=ifft(fft(qA_plus(:,i+1)).*coff1);
end
%%%%%计算Q%%%%%
QQ=zeros(Nx,1);
%  for i=N1_f:-1:1
%用qA来算
% QQ=qA_plus(:,i).*qA(:,Nf-i+1);
% 用qB来算
% qB_plus_tixing1=qB_plus(:,:,r,N1_f)*qA(r,Nf);
  qB_plus_tixing1=qB_plus(:,:,:,N1_f).*4*pi.*qB(:,:,:,1);
% 验证   
%  qB_plus_tixing1=qB_plus(:,:,:,i).*4*pi.*qB(:,:,:,N1_f+1-i);
 qB_sita1=permute(sum(qB_plus_tixing1,2),[1,3,2])*dt1;
  for r=1:Nx
      QQ(r)=sum(qB_sita1(:,r).*wsita);
  end
    QQ=real(fft(QQ)/Nx);
    Q=QQ(1);
        fprintf('%.6e \n',Q);
%  end

%把qA_plus和qB_plus变成正向存储
 qA_plus=fliplr(qA_plus);
 qB_plus=qB_plus(:,:,:,end:-1:1);
    
 %计算faia,faib
    qA_muti=qA.*qA_plus;qB_muti=qB.*qB_plus;
    sum1=0;sum2=0;
    for m=1:Nf
        sum1=sum1+qA_muti(:,m);
    end
    qB_muti_sum=zeros(Nx,N1_f);
    qB_sita2=zeros(Nsita,Nx,N1_f);
    for i=1:N1_f
        for r=1:Nx
        qB_muti_sum(r,i)=sum(sum(qB_muti(:,:,r,i),2)*dt1.*wsita);
        end
    end  
    for n=1:N1_f
        sum2=sum2+qB_muti_sum(:,n);
    end
    faia=tao*(-5/8*(qA_muti(:,1)+qA_muti(:,Nf))+1/6*(qA_muti(:,2)+qA_muti(:,Nf-1))-1/24*(qA_muti(:,3)+qA_muti(:,Nf-2))+sum1);
    faia=faia/Q;
    
    faib=tao*(-5/8*(qB_muti_sum(:,1)+qB_muti_sum(:,N1_f))+1/6*(qB_muti_sum(:,2)+qB_muti_sum(:,N1_f-1))-1/24*(qB_muti_sum(:,3)+qB_muti_sum(:,N1_f-2))+sum2);
    faib=4*pi*faib/Q;
%     fprintf('%d \n' ,faia);
%     fprintf('\n\n\n ' );
%     fprintf('%d \n' ,faib);
%     fprintf('\n\n\n ' );
%      fprintf('%d \n' ,faia+faib);
%      fprintf('\n\n\n ' );
    
%     计算S
  mid_result=zeros(3,3,Nsita,Nphi,Nx,N1_f);
    for i=1:N1_f
        for r=1:Nx
            for k1=1:Nsita
                for n=1:Nphi
                    mid_result(:,:,k1,n,r,i)=qB(k1,n,r,i)*jzb(:,:,k1,n)*qB_plus(k1,n,r,i);
                end
            end
        end 
    end
    mid_1=permute(sum(mid_result,4)*dt1,[1,2,3,5,6,4]);
    sum3=zeros(3,3,Nsita,Nx,N1_f);

    for i=1:N1_f
        for r=1:Nx
            for k1=1:Nsita
                sum3(:,:,k1,r,i)=mid_1(:,:,k1,r,i)*wsita(k1);
            end
        end
    end
    mid_2=permute(sum(sum3,3),[1,2,4,5,3]);
    sum4=0;
    for n=1:N1_f
        sum4=sum4+mid_2(:,:,:,n);
    end
    S=tao*(-5/8*(mid_2(:,:,:,1)+mid_2(:,:,:,N1_f))+1/6*(mid_2(:,:,:,2)+mid_2(:,:,:,N1_f-1))-1/24*(mid_2(:,:,:,3)+mid_2(:,:,:,N1_f-2))+sum4);
    S=4*pi*S/Q;
  
    %终止条件
   upde_S=zeros(3,3);
   ss=S-M/eita/kN;
   for i=1:3
       for j=1:3
           upde_S(i,j)=norm(permute(ss(i,j,:),[3,1,2]))^2;
       end
   end
   step=step+1;
   error_M=sqrt(sum(upde_S(:))/9);
   fprintf('迭代步数:%d \t 误差1：%.15f \t 误差2 : %.15f \t 误差3:%.15f \n', step,norm(faia+faib-1),norm(faia-faib-2/kN*u_dis),error_M);
    if(norm(faia+faib-1)>10^-5||norm(faia-faib-2/kN*u_dis)>10^-5||error_M>10^-5 )
       u_plus=u_plus+(faia+faib-1)*1;
       u_dis=u_dis+1*(faia-faib-2*u_dis./kN); 
       M=M+1*(S-M/eita/kN);
       
       
       fprintf('\n ' );
    else
           break;
     end
 end
 
 u_dis_1=ifft(u_dis.^2);u_plus_1=ifft(u_plus); M_dot=zeros(Nx,1);
 for r=1:Nx
M_dot(r)=sum(sum(M(:,:,r).*M(:,:,r)));
 end
M_dot=ifft(M_dot);M_dot=M_dot(1);
H=u_dis_1(1)/kN-u_plus_1(1)-log(Q)+M_dot/2/eita/kN;