%求除e^(im*phi)外的Y_lm(sita,phi),即h(sita)
function  YPlm=prepare(N1,BNC,N3C,N2C)

% legendre=zeros(N1,N1);
% spherical=zeros(N1,N1,N1);
% legendre(1,1)=1;
% legendre(2,1)=0;
% legendre(2,2)=1;
% for n=2:N1-1
%     for j=1:n+1
%         if j>1
%             legendre(n+1,j)=(2*(n-1)+1)*1.0/(n)*legendre(n,j-1)-(n-1)*1.0/n*legendre(n-1,j);%定义勒让德函数
%         else
%             legendre(n+1,j)=-(n-1)*1.0/n*legendre(n-1,j);
%         end
%     end
% end
% for l=1:BNC
%     for i=1:N1
%         spherical(l,1,i)=legendre(l,i);%定义连带勒让德函数（m=0）
%     end
% end
% for l=1:BNC
%     for m=2:l
%         for i=1:N1-1
%             spherical(l,m,i)=i*spherical(l,m-1,i+1);%定义除（1-x^2）^(m/2)的连带勒让德函数(m≠0)
%         end
%     end
% end
gauss=get_gaussian(N2C);
cossita=zeros(N2C,1);sinsita=zeros(N2C,1);%cosphi=zeros(N3C,1);sinphi=zeros(N3C,1);
for i=1:N2C
    if i<=(N2C+1)/2
        cossita(i)=gauss(i,1);
    else
        cossita(i)=-gauss(N2C+1-i,1);
    end
    sinsita(i)=sqrt(1-cossita(i)*cossita(i));
end
%for j=1:N3C
%    cosphi(j)=cos(2*pi/(N3C)*(j-1));
%    sinphi(j)=sin(2*pi/(N3C)*(j-1));
%end
% jj = 0:1:N3C-1;
% cosphi = cos(2*pi.*jj/N3C);
% sinphi = sin(2*pi.*jj/N3C);

%Y_star=zeros((BNC+1)*(BNC+2)/2,N2C,N3C);
YPlm=zeros((BNC+1)*(BNC+2)/2,N2C);

for l=1:BNC+1
    for m=1:l
        for i=1:N2C
            YPlm(((l-1)*l/2+m),i)=spherical_harmonic(N1,BNC,l,m,cossita(i));%求除e^(im*phi)外的Y_lm(sita,phi),即h(sita)
             %for j=1:N3C
              %  Y_star((l-1)*l/2+m,i,j,1)=spherical_harmonic(N1,BNC,l,m,cossita(i))*cos(m*2*pi/(N3C)*(j-1));
 		       %Y_star((l-1)*l/2+m,i,j,2)=-spherical_harmonic(N1,BNC,l,m,cossita(i))*sin(m*2*pi/(N3C)*(j-1));%定义Y_lm(sita,phi)
%              end
        end
    end
end
end


