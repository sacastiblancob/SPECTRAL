%Creating matrices for solve the system


Globals1D

%
%  Create mass matrix
%
M = W;
if deal==1
    Md = Wd;
end
% M=zeros(Np,Np);
% for i=1:Np
%   %M(i,i)=wg(i)*dx1/2.0;
%   M(i,i) = w(i);
% end

%
%  Create Advection matrix
%
A = Dr.*w';
if deal==1
    Ad = Drd.*wd';
end
% A=zeros(Np,Np);
% for i=1:Np
%   for j=1:Np
%     A(i,j)=Dr(i,j)*w(i);
%   end
% end

%
%  Create Stiffness matrix
%
%S=zeros(Np,Np);
%S = Dr'*(A);
%S = Dr'*Dr;
ep = 1/(dt*nu);
S = Dr'*W*Dr;
% S=zeros(Np,Np);
% for i=1:Np
%   for j=1:Np
%     for k=1:Np
%       S(i,j)=S(i,j)+(Dr(k,i)*Dr(k,j)*w(k)*J(1,1));
%     end
%   end
% end
%S = Dr;

%
%  Global Assembly
%

AG=zeros(Ns,Ns);
MG=zeros(Ns,Ns);
if deal==1
    AGd = zeros(Nsd,Nsd);
    MGd = zeros(Nsd,Nsd);
end
SG=zeros(Ns,Ns);
ap=0;
for k=1:K
  for i=1:Np
    for j=1:Np
      SG(Np*(k-1)+i-ap,Np*(k-1)+j-ap) = SG(Np*(k-1)+i-ap,Np*(k-1)+j-ap) + ...
          S(i,j)*(1/J(i,k)); % + ep*M(i,j)*J(i,k);
        
      AG(Np*(k-1)+i-ap,Np*(k-1)+j-ap) = AG(Np*(k-1)+i-ap,Np*(k-1)+j-ap) + A(i,j); %*(J(i,k));

      MG(Np*(k-1)+i-ap,Np*(k-1)+j-ap) = MG(Np*(k-1)+i-ap,Np*(k-1)+j-ap) + M(i,j)*J(i,k);
    end
  end
  ap=ap+1;
end

if deal==1
    ap = 0;
    for k=1:K
      for i=1:Npd
        for j=1:Npd
          AGd(Npd*(k-1)+i-ap,Npd*(k-1)+j-ap) = AGd(Npd*(k-1)+i-ap,Npd*(k-1)+j-ap) + Ad(i,j); %*(J(i,k));
          
          MGd(Npd*(k-1)+i-ap,Npd*(k-1)+j-ap) = MGd(Npd*(k-1)+i-ap,Npd*(k-1)+j-ap) + Md(i,j)*Jd(i,k);
        end
      end
      ap=ap+1;
    end
end

%
%  Inversion of matrix MG
%
IMG = MG;
for i=1:Ns
    IMG(i,i) = 1/MG(i,i);
end

if deal==1
    IMGd = MGd;
    for i=1:Ns
        IMGd(i,i) = 1/MGd(i,i);
    end
end
%
%  For solve in time
%
AG = IMG*AG;
if deal == 1
    AGd = IMGd*AGd;
end
%SG = IMG*SG;
%GM = AG
% AG(1,:) = 0.0;
% AG(Ns,:) = 0.0;
% AG(1,1) = 1;
% AG(Ns,Ns) = 1;
% SG(1,:) = 0.0;
% SG(Ns,:) = 0.0;
% SG(1,1) = 1;
% SG(Ns,Ns) = 1;





