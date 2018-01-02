%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% 8/16/2017                                    %
% Author: Jaman Mohebujjaman                %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all


 
nModes = 8;
endTime = 0.332;
snapIndex = 1000;
 

My_tol_DDROM = 1e-4;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
load snapshotData35Kdt002SV_Re100
%load Gsnap_SV35K_r8_d16_N16_166
load Gsnap_SV35K_r8_d16_N16_166


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate A and B tilde

%Import offline matrices
load ROMtestSV35K_N16_166  Snapshots MassROM StiffROM TriLinROM2 NLlift NLdrag vdmass vdstiff vlmass vlstiff GlobalV PhiR MassMatrix T dt nu BalanceTable nodeco GradDivMatrix elnode
tic 

endTimestep = 166;

N = nModes;
MassROM = MassROM(1:N,1:N);
StiffROM = StiffROM(1:N,1:N);
%GradDivROM = GradDivROM(1:N,1:N);
TriLinROM2 = TriLinROM2(1:N,1:N,1:N);
NLlift = NLlift(1:N,1:N);
NLdrag = NLdrag(1:N,1:N);
vdmass = vdmass(1:N);
vdstiff = vdstiff(1:N);
vlmass = vlmass(1:N);
vlstiff = vlstiff(1:N);

r=nModes;

% L2 project the initial condition (held in GlobalV) into ROM basis - coeff
% vector is put into "velInit"

for timestep = 1:endTimestep
    vvv = Snapshots(:,1000+timestep); 
    RHS(:,timestep) = zeros(r,1); 

    for i = 1:r
        RHS(i,timestep) = vvv' * (MassMatrix * PhiR(:,i) );
    end
    
    A = MassROM;
    RHS(:,timestep) = RHS(:,timestep) - A(:,1)*1;
    
    A(1,:) = 0;
    A(:,1) = 0;
    A(1,1) = 1;
    RHS(1,timestep) = 1;
    velInit(:,timestep) = A  \ RHS(:,timestep);
end


    

MySystemMatrixForAB = [];

%loop over a single period

for timestep = 1:endTimestep
   
    MyLocalmatForAB = zeros(r,r^2+r^3);
    
    for row = 1:r
        MyLocalmatForAB(row,(row-1)*r+(1:r)) = velInit(1:r,timestep)'; 
    end
    
    vec = zeros(r^2,1);
    
    for i = 1:r
        for j = 1:r
            k = (i-1)*r+j;
            vec(k,1) = velInit(i,timestep)'* velInit(j,timestep);
        end
    end
    
    for row = 1:r
        MyLocalmatForAB(row,row*r^2+(1:r^2)) = vec;
    end
    
    MySystemMatrixForAB = [MySystemMatrixForAB;MyLocalmatForAB];

end

[myU,mySigma,myV] = svd(MySystemMatrixForAB);
s = diag(mySigma);


%Number of singular values bigger than tol

mycountforAB = 0;
for i =1:length(s)
    if s(i)>= My_tol_DDROM;
        mycountforAB = mycountforAB +1;
    end
end

%Truncated SVD algorithm

ABtildeU = myU(:,1:mycountforAB);

ABtildeS = mySigma(1:mycountforAB,1:mycountforAB)

ABtildeV = myV(:,1:mycountforAB);

tildeD = ABtildeU * ABtildeS * ABtildeV';

 
F = [];
ExtractRightHandSide = Gsnap;
for i = 1:endTimestep
    F = [F;ExtractRightHandSide(1:r,i)];
end

X = ABtildeV*inv(ABtildeS)*ABtildeU'*F;

%check the l2 norm of the residual
norm(F-MySystemMatrixForAB*X)


XA = X(1:r^2,1);
XB = X(r^2+1:end,1);


ABtildeA = reshape(XA,[r,r])'; 

ABtildeB = zeros(r,r,r);

for i = 1:r
    for j = 1:r
            ABtildeB(i,:,j) = XB((1:r)+(j-1)*r+(i-1)*r^2);
    end
end
toc

%save the matrices

save ABtilde_N16_r8_d16_166

