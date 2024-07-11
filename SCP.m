%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program can calculat eigenvalu befor and after pertobation     %
%   pertobation in this program is space charge effect in HQFC Matrix   %
%   I  use Morita and toprek paper                                      %
%                       5/24/2013                                       %
%                       N S Mirian                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
syms N k KK s DKx DKy

%eq.16 of Mortia and Iwashita STAB paper
R=[cos(k*s/(2*N)) 0 sin(k*s/(2*N)) 0
0 cos(k*s/(2*N)) 0 sin(k*s/(2*N))
-sin(k*s/(2*N)) 0 cos(k*s/(2*N)) 0
0 -sin(k*s/(2*N)) 0 cos(k*s/(2*N))];

%eq.17 of Morita and Iwashita STAB paper
M=[1 s/N 0 0
-KK*s/N 1 0 0
0 0 1 s/N
0 0 KK*s/N 1];

C1=R*M;
%--------------------------------------
%calcualte 4 eigenvecotrs (v) and 
%4 eignvalues (d) of C1 matrix
%eq. 18 of Morita dn Iwashita STAB paper
[v,d]=eig(C1);
%first eigenvalue to power of N
d1N=d(1,1)^N;
%take limit when N->00
lambda1=limit(d1N, N, inf);
disp('------------------lambda01------------------------------------ ')
pretty(lambda1)
%simple(d(1));

d2N=d(2,2)^N;
lambda2=limit(d2N, N, inf);
disp('------------------lambda02------------------------------------- ')
pretty(lambda2)
%simple(d(2));

d3N=d(3,3)^N;
lambda3=limit(d3N, N, inf);
disp('------------------lambda03------------------------------------- ')
pretty(lambda3)
%simple(d(3));

d4N=d(4,4)^N;
lambda4=limit(d4N, N, inf);
disp('------------------lambda04------------------------------------- ')
pretty(lambda4)
%------------------------------------
% now we operate perturbation-----------
%------Toprek paper ---------------------- 
% --- Dkx and Dky are space chrge effect in eq(4a) & eq(4b)
deltaM=[0 0 0 0 
    DKx*s/N 0 0 0
    0 0 0 0
    0 0 DKy*s/N 0];
disp('-----------[delta_K]=Rotitian Matrix * space-chrge Matrix -- ')
deltaC1=R*deltaM

VV1=[v(1,1); v(2,1); v(3,1); v(4,1)];
VV2=[v(1,2); v(2,2); v(3,2); v(4,2)];
VV3=[v(1,3); v(2,3); v(3,3); v(3,4)];
VV4=[v(1,4); v(2,4); v(3,4); v(4,4)];
disp('-----------------lambda(s/N) are: ')
d1=d(1,1)+VV1.'*deltaC1*VV1,
d2=d(2,2)+VV2.'*deltaC1*VV2,
d3=d(3,3)+VV3.'*deltaC1*VV3,
d4=d(4,4)+VV4.'*deltaC1*VV4,


 %----------------
    %first eigenvalue to power of N
d1N=d1^N;
%take limit when N->00
lambda1=limit(d1N, N, inf);
L1=simple(lambda1);
disp('-----------------lambda1= ---------------------- ')
pretty(L1)
%-----------------
d2N=d2^N;
lambda2=limit(d2N, N, inf);
L2=simple(lambda2);
disp('-----------------lambda2= ---------------------- ')
pretty(L2)
%-----------------
d3N=d3^N;
lambda3=limit(d3N, N, inf);
%pretty(lambda3)
L3=simple(lambda3);
disp('-----------------lambda3= ---------------------- ')
pretty(L3)
%-----------------
d4N=d4^N;
lambda4=limit(d4N, N, inf);
L4=simple(lambda4);
disp('-----------------lambda4= ---------------------- ')
pretty(L4)
%-----------------
