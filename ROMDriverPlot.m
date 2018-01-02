%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% 7/30/2017                                    %
% Author: Jaman Mohebujjaman                %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all


 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
load snapshotData35Kdt002SV_Re100
%load ABtilde_N32_r8_d9_166.mat
%load ABtilde_N32_r8_d9.mat
%load Gsnap_SV35K_r8_d32_N32

%load ABtilde_N32_r8_d8_166
%load ABtilde_N16_r8_d16
%load ABtilde_N32_r8_d9_1500.mat
%load ABtilde_N32_r8_d32_166.mat
%load ABtilde_N32_r8_d24_166.mat
%load ABtilde_N32_r8_d16_1500.mat
%load ABtilde_N32_r6_d18_166.mat
%load ABtilde_N32_r4_d12_166.mat
%load Gsnap_SV35K_r8_d16_N16
%load Gsnap_SV35K_r8_d9_N32_1500

load ABtilde_N16_r8_d16_166
%load Gsnap_SV35K_r8_d16_N16_2
load Atilde_N16_r8_d16

%load ABtilde_N9_r8_d9_1500



%load ABtilde_N9_r8_d9_1500
%load Atildeonly
%load Gsnap_SV35K_r8_d8_N32_166
%load Gsnap_SV35K_r8_d24_N32_n166.mat
%load Gsnap_SV35K_r6_d18_N32.mat
%load Gsnap_SV35K_r4_d12_N32
%load Gsnap_SV35K_r8_d32_N32_1500.mat
% load ABtilde_N32_r8_d32_1500
% load Gsnap_SV35K_r8_d32_N32_1500

% load ABtilde_N32_r8_d16_1500
% load Gsnap_SV35K_r8_d16_N32_1500
%load Gsnap_SV35K_r8_d32_N32.mat
%load ABtilde_N32_r8_d8_166
%load ABtilde_N32_r8_d16
%load Gsnap_SV35K_r8_d8_N8
load DNSProjectionMatrix_r8
%load DNStable10
%load ABtilde_N32_r8_d32

nModes = 8;
endTime = 3.0;
snapIndex = 1000;
 
%endTimestep = 166;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate Error
%%
% Atilde only
load ROMtestSV35K_N16_166  Snapshots MassROM StiffROM  TriLinROM2 NLlift NLdrag vdmass vdstiff vlmass vlstiff GlobalV PhiR MassMatrix T dt nu BalanceTable nodeco GradDivMatrix elnode
 

GlobalV = Snapshots(:,snapIndex);
%GlobalVPrev = Snapshots(:,snapIndex-1);


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


% L2 project the initial condition (held in GlobalV) into ROM basis - coeff
% vector is put into "velInit"
vvv = GlobalV;
RHS = zeros(N,1);
 for i=1:N
     RHS(i) = vvv' * (MassMatrix * PhiR(:,i) );
 end
A = MassROM;
RHS = RHS - A(:,1)*1;
A(1,:)=0;
A(:,1)=0;
A(1,1)=1;
RHS(1)=1;
velInit = A  \ RHS;

velPrevPrev=velInit;


dataTableAtilde=[];
dataTableDNS = [];

% simulation parameters

T=endTime;
numTimeSteps =  round(T/dt);


solns = zeros(N,numTimeSteps+1);
solns(:,1) = velInit; 

% MyL2norm = zeros(numTimeSteps, 1);
%initial condition

 ts  = 1;
% 
 b = 1.0/dt * MassROM * velPrevPrev;
% % build matrix
 NLmat = 0*MassROM;
 for k=1:N
    % NLmat = NLmat + velPrevPrev(k)*(TriLinROM2(:,:,k));
     NLmat = NLmat + velPrevPrev(k)*(TriLinROM2(:,:,k));
 end
  % A = 1.0/dt * MassROM + nu*StiffROM + NLmat +myAnew;
 A = 1.0/dt * MassROM + nu*StiffROM + NLmat +Atilde;     
% 
 A(1,:) = 0;
 A(1,1) = 1;
 b(1) = 1;
%     
% % solve the linear system
 velSoln = A \ b;
 
% compute L2 norm of error

% diff = velSoln-DNSProjectionMatrix(:,1001);

%  MyL2norm (1) = sqrt(diff'* (MassROM * diff));
 %MyL2norm (1) = norm(diff);

% First put ROM solution back in the FE basis
%     u_ROM = 0*GlobalV;
%     for j=1:r
%       u_ROM = u_ROM + velSoln(j)*PhiR(:,j);
%     end
%     
%     diff = u_ROM - Snapshots(:,1001); 
%     MyL2norm (1) = sqrt(diff' * (MassMatrix * diff));
%     
%lift = 1/dt (u-uprev,v_l) + nu(grad u,grad v_l) + b(u,u,v_l)
vtDNS = 1/dt * (DNSProjectionMatrix(:,1000+ts)-DNSProjectionMatrix(:,1000));
vt = 1/dt*(velSoln - velPrevPrev);
lift = -20*( vt'*vlmass' + nu*velSoln'*vlstiff' + velSoln' * NLlift * velSoln );

liftDNS = -20*( vtDNS'*vlmass' + nu*DNSProjectionMatrix(:,1000+ts)'*vlstiff' + DNSProjectionMatrix(:,1000+ts)' * NLlift * DNSProjectionMatrix(:,1000+ts) );

drag = -20*( vt'*vdmass' + nu*velSoln'*vdstiff' + velSoln' * NLdrag * velSoln );


energy = 1/2 * sqrt((velSoln' * (MassROM * velSoln) ));
energyDNS = 1.0/2.0 * sqrt((DNSProjectionMatrix(:,1000+ts)' *(MassROM * DNSProjectionMatrix(:,1000+ts)) ));
dragDNS = -20*( vtDNS'*vdmass' + nu*DNSProjectionMatrix(:,1000+ts)'*vdstiff' + DNSProjectionMatrix(:,1000+ts)' * NLdrag * DNSProjectionMatrix(:,1000+ts) );
%divL2 = sqrt((velSoln' * (GradDivROM * velSoln) ));
display([ num2str(ts*dt) '  ' num2str(lift) '   '  num2str(drag) '   ' num2str(energy) ])
dataTableAtilde = [dataTableAtilde; ts*dt, lift, drag,   energy];
dataTableDNS = [dataTableDNS; ts*dt, liftDNS, dragDNS, energyDNS];
%      
%     
 solns(:,ts+1) = velSoln;
 velPrev = velSoln;



%BDF2 with linear scheme
   
for ts=2:numTimeSteps
    %RHS
    
    
    
    b = 2/dt * MassROM * velPrev -0.5/dt*MassROM*velPrevPrev ;
    % build matrix
    NLmat = 0*MassROM;
    for k=1:N
       % NLmat = NLmat + (2*velPrev(k)-velPrevPrev(k))*(TriLinROM2(:,:,k) );
       NLmat = NLmat + (2*velPrev(k)-velPrevPrev(k))*(TriLinROM2(:,:,k));
    end
    A = 1.5/dt * MassROM + nu*StiffROM + NLmat +Atilde;
   % A = 1.5/dt * MassROM + nu*StiffROM + NLmat +myAnew;

    A(1,:)=0;
    A(1,1)=1;
    b(1)=1;
    
    % solve the linear system
    velSoln = A \ b;

    
    
    % compute L2 norm error
    
%     diff = velSoln-DNSProjectionMatrix(:,1000+ts);
    %MyL2norm (ts) = norm(diff);
%     MyL2norm (ts) = sqrt(diff'* (MassROM * diff));
%     
    % First put ROM solution back in the FE basis
%     u_ROM = 0*GlobalV;
%     for j=1:N
%       u_ROM = u_ROM + velSoln(j)*PhiR(:,j);
%     end
%     
%     diff = u_ROM - Snapshots(:,ts+1000); 
%     MyL2norm (ts) = sqrt(diff' * (MassMatrix * diff));
%     
    
    %     
    %lift = 1/dt (u-uprev,v_l) + nu(grad u,grad v_l) + b(u,u,v_l)
    vtDNS = 1/dt * (DNSProjectionMatrix(:,1000+ts)-DNSProjectionMatrix(:,999+ts));
     vt = 1/dt*(velSoln - velPrev);
     lift = -20*( vt'*vlmass' + nu*velSoln'*vlstiff' + velSoln' * NLlift * velSoln );
     drag = -20*( vt'*vdmass' + nu*velSoln'*vdstiff' + velSoln' * NLdrag * velSoln );


     energy = 1/2 * sqrt((velSoln' * (MassROM * velSoln) ));
     liftDNS = -20*( vtDNS'*vlmass' + nu*DNSProjectionMatrix(:,1000+ts)'*vlstiff' + DNSProjectionMatrix(:,1000+ts)' * NLlift * DNSProjectionMatrix(:,1000+ts) );
     energyDNS = 1.0/2.0 * sqrt((DNSProjectionMatrix(:,1000+ts)' *(MassROM * DNSProjectionMatrix(:,1000+ts)) ));
     dragDNS = -20*( vtDNS'*vdmass' + nu*DNSProjectionMatrix(:,1000+ts)'*vdstiff' + DNSProjectionMatrix(:,1000+ts)' * NLdrag * DNSProjectionMatrix(:,1000+ts) );
     %divL2 = sqrt((velSoln' * (GradDivROM * velSoln) ));
     display([ num2str(ts*dt) '  ' num2str(lift) '   '  num2str(drag) '     ' num2str(energy) ])
     dataTableAtilde = [dataTableAtilde; ts*dt, lift, drag,   energy];
     dataTableDNS = [dataTableDNS; ts*dt,liftDNS, dragDNS, energyDNS];
        
    solns(:,ts+1)=velSoln;
    velPrevPrev=velPrev;
    velPrev = velSoln;
    
    
    

end
display('AB tilde');
% aveErrorDD = sum(MyL2norm)/length(MyL2norm)


load ROMtestSV35K_N16_166  Snapshots MassROM StiffROM  TriLinROM2 NLlift NLdrag vdmass vdstiff vlmass vlstiff GlobalV PhiR MassMatrix T dt nu BalanceTable nodeco GradDivMatrix elnode
 



GlobalV = Snapshots(:,snapIndex);
%GlobalVPrev = Snapshots(:,snapIndex-1);


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




% L2 project the initial condition (held in GlobalV) into ROM basis - coeff
% vector is put into "velInit"
vvv = GlobalV;
RHS = zeros(N,1);
 for i=1:N
     RHS(i) = vvv' * (MassMatrix * PhiR(:,i) );
 end
A = MassROM;
RHS = RHS - A(:,1)*1;
A(1,:)=0;
A(:,1)=0;
A(1,1)=1;
RHS(1)=1;
velInit = A  \ RHS;

velPrevPrev=velInit;


dataTableAB=[];
%dataTableDNS = [];

% simulation parameters

T=endTime;
numTimeSteps =  round(T/dt);


solns = zeros(N,numTimeSteps+1);
solns(:,1) = velInit; 

% MyL2norm = zeros(numTimeSteps, 1);
%initial condition

 ts  = 1;
% 
 b = 1.0/dt * MassROM * velPrevPrev;
% % build matrix
 NLmat = 0*MassROM;
 for k=1:N
    % NLmat = NLmat + velPrevPrev(k)*(TriLinROM2(:,:,k));
     NLmat = NLmat + velPrevPrev(k)*(TriLinROM2(:,:,k) + ABtildeB(:,:,k));
 end
  % A = 1.0/dt * MassROM + nu*StiffROM + NLmat +myAnew;
 A = 1.0/dt * MassROM + nu*StiffROM + NLmat +ABtildeA;     
% 
 A(1,:) = 0;
 A(1,1) = 1;
 b(1) = 1;
%     
% % solve the linear system
 velSoln = A \ b;
 
% compute L2 norm of error

% diff = velSoln-DNSProjectionMatrix(:,1001);

%  MyL2norm (1) = sqrt(diff'* (MassROM * diff));
 %MyL2norm (1) = norm(diff);

% First put ROM solution back in the FE basis
%     u_ROM = 0*GlobalV;
%     for j=1:r
%       u_ROM = u_ROM + velSoln(j)*PhiR(:,j);
%     end
%     
%     diff = u_ROM - Snapshots(:,1001); 
%     MyL2norm (1) = sqrt(diff' * (MassMatrix * diff));
%     
%lift = 1/dt (u-uprev,v_l) + nu(grad u,grad v_l) + b(u,u,v_l)
vt = 1/dt*(velSoln - velPrevPrev);
lift = -20*( vt'*vlmass' + nu*velSoln'*vlstiff' + velSoln' * NLlift * velSoln );
drag = -20*( vt'*vdmass' + nu*velSoln'*vdstiff' + velSoln' * NLdrag * velSoln );


energy = 1/2 * sqrt((velSoln' * (MassROM * velSoln) ));
%energyDNS = 1.0/2.0 * sqrt((DNSProjectionMatrix(:,ts)' *(MassROM * DNSProjectionMatrix(:,ts)) ));
%divL2 = sqrt((velSoln' * (GradDivROM * velSoln) ));
display([ num2str(ts*dt) '  ' num2str(lift) '   '  num2str(drag) '   ' num2str(energy) ])
dataTableAB = [dataTableAB; ts*dt, lift, drag,   energy];
%dataTableDNS = [dataTableDNS; ts*dt, energyDNS];
%      
%     
 solns(:,ts+1) = velSoln;
 velPrev = velSoln;



%BDF2 with linear scheme
   
for ts=2:numTimeSteps
    %RHS
    
    
    
    b = 2/dt * MassROM * velPrev -0.5/dt*MassROM*velPrevPrev ;
    % build matrix
    NLmat = 0*MassROM;
    for k=1:N
       % NLmat = NLmat + (2*velPrev(k)-velPrevPrev(k))*(TriLinROM2(:,:,k) );
       NLmat = NLmat + (2*velPrev(k)-velPrevPrev(k))*(TriLinROM2(:,:,k) +ABtildeB(:,:,k));
    end
    A = 1.5/dt * MassROM + nu*StiffROM + NLmat +ABtildeA;
   % A = 1.5/dt * MassROM + nu*StiffROM + NLmat +myAnew;

    A(1,:)=0;
    A(1,1)=1;
    b(1)=1;
    
    % solve the linear system
    velSoln = A \ b;

    
    
    % compute L2 norm error
    
%     diff = velSoln-DNSProjectionMatrix(:,1000+ts);
    %MyL2norm (ts) = norm(diff);
%     MyL2norm (ts) = sqrt(diff'* (MassROM * diff));
%     
    % First put ROM solution back in the FE basis
%     u_ROM = 0*GlobalV;
%     for j=1:N
%       u_ROM = u_ROM + velSoln(j)*PhiR(:,j);
%     end
%     
%     diff = u_ROM - Snapshots(:,ts+1000); 
%     MyL2norm (ts) = sqrt(diff' * (MassMatrix * diff));
%     
    
    %     
    %lift = 1/dt (u-uprev,v_l) + nu(grad u,grad v_l) + b(u,u,v_l)
     vt = 1/dt*(velSoln - velPrev);
     lift = -20*( vt'*vlmass' + nu*velSoln'*vlstiff' + velSoln' * NLlift * velSoln );
     drag = -20*( vt'*vdmass' + nu*velSoln'*vdstiff' + velSoln' * NLdrag * velSoln );


     energy = 1/2 * sqrt((velSoln' * (MassROM * velSoln) ));
     %energyDNS = 1.0/2.0 * sqrt((DNSProjectionMatrix(:,ts)' *(MassROM * DNSProjectionMatrix(:,ts)) ));
     %divL2 = sqrt((velSoln' * (GradDivROM * velSoln) ));
     display([ num2str(ts*dt) '  ' num2str(lift) '   '  num2str(drag) '     ' num2str(energy) ])
     dataTableAB = [dataTableAB; ts*dt, lift, drag,   energy];
     %dataTableDNS = [dataTableDNS; ts*dt, energyDNS];
        
    solns(:,ts+1)=velSoln;
    velPrevPrev=velPrev;
    velPrev = velSoln;
    
    
    

end
display('GROM');
load ROMtestSV35K_N16_166  Snapshots MassROM StiffROM  TriLinROM2 NLlift NLdrag vdmass vdstiff vlmass vlstiff GlobalV PhiR MassMatrix T dt nu BalanceTable nodeco GradDivMatrix elnode
 



GlobalV = Snapshots(:,snapIndex);
%GlobalVPrev = Snapshots(:,snapIndex-1);


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




% L2 project the initial condition (held in GlobalV) into ROM basis - coeff
% vector is put into "velInit"
vvv = GlobalV;
RHS = zeros(N,1);
 for i=1:N
     RHS(i) = vvv' * (MassMatrix * PhiR(:,i) );
 end
A = MassROM;
RHS = RHS - A(:,1)*1;
A(1,:)=0;
A(:,1)=0;
A(1,1)=1;
RHS(1)=1;
velInit = A  \ RHS;

velPrevPrev=velInit;


dataTableGROM=[];
%dataTableDNS = [];

% simulation parameters

T=endTime;
numTimeSteps =  round(T/dt);


solns = zeros(N,numTimeSteps+1);
solns(:,1) = velInit; 

% MyL2norm = zeros(numTimeSteps, 1);
%initial condition

 ts  = 1;
% 
 b = 1.0/dt * MassROM * velPrevPrev;
% % build matrix
 NLmat = 0*MassROM;
 for k=1:N
    % NLmat = NLmat + velPrevPrev(k)*(TriLinROM2(:,:,k));
     NLmat = NLmat + velPrevPrev(k)*(TriLinROM2(:,:,k));
 end
  % A = 1.0/dt * MassROM + nu*StiffROM + NLmat +myAnew;
 A = 1.0/dt * MassROM + nu*StiffROM + NLmat;     
% 
 A(1,:) = 0;
 A(1,1) = 1;
 b(1) = 1;
%     
% % solve the linear system
 velSoln = A \ b;
 
% compute L2 norm of error

% diff = velSoln-DNSProjectionMatrix(:,1001);

%  MyL2norm (1) = sqrt(diff'* (MassROM * diff));
 %MyL2norm (1) = norm(diff);

% First put ROM solution back in the FE basis
%     u_ROM = 0*GlobalV;
%     for j=1:r
%       u_ROM = u_ROM + velSoln(j)*PhiR(:,j);
%     end
%     
%     diff = u_ROM - Snapshots(:,1001); 
%     MyL2norm (1) = sqrt(diff' * (MassMatrix * diff));
%     
%lift = 1/dt (u-uprev,v_l) + nu(grad u,grad v_l) + b(u,u,v_l)
vt = 1/dt*(velSoln - velPrevPrev);
lift = -20*( vt'*vlmass' + nu*velSoln'*vlstiff' + velSoln' * NLlift * velSoln );
drag = -20*( vt'*vdmass' + nu*velSoln'*vdstiff' + velSoln' * NLdrag * velSoln );


energy = 1/2 * sqrt((velSoln' * (MassROM * velSoln) ));
%energyDNS = 1.0/2.0 * sqrt((DNSProjectionMatrix(:,ts)' *(MassROM * DNSProjectionMatrix(:,ts)) ));
%divL2 = sqrt((velSoln' * (GradDivROM * velSoln) ));
display([ num2str(ts*dt) '  ' num2str(lift) '   '  num2str(drag) '   ' num2str(energy) ])
dataTableGROM = [dataTableGROM; ts*dt, lift, drag,   energy];
%dataTableDNS = [dataTableDNS; ts*dt, energyDNS];
%      
%     
 solns(:,ts+1) = velSoln;
 velPrev = velSoln;



%BDF2 with linear scheme
   
for ts=2:numTimeSteps
    %RHS
    
    
    
    b = 2/dt * MassROM * velPrev -0.5/dt*MassROM*velPrevPrev ;
    % build matrix
    NLmat = 0*MassROM;
    for k=1:N
       % NLmat = NLmat + (2*velPrev(k)-velPrevPrev(k))*(TriLinROM2(:,:,k) );
       NLmat = NLmat + (2*velPrev(k)-velPrevPrev(k))*(TriLinROM2(:,:,k));
    end
    A = 1.5/dt * MassROM + nu*StiffROM + NLmat;
   % A = 1.5/dt * MassROM + nu*StiffROM + NLmat +myAnew;

    A(1,:)=0;
    A(1,1)=1;
    b(1)=1;
    
    % solve the linear system
    velSoln = A \ b;

    
    
    % compute L2 norm error
    
%     diff = velSoln-DNSProjectionMatrix(:,1000+ts);
    %MyL2norm (ts) = norm(diff);
%     MyL2norm (ts) = sqrt(diff'* (MassROM * diff));
%     
    % First put ROM solution back in the FE basis
%     u_ROM = 0*GlobalV;
%     for j=1:N
%       u_ROM = u_ROM + velSoln(j)*PhiR(:,j);
%     end
%     
%     diff = u_ROM - Snapshots(:,ts+1000); 
%     MyL2norm (ts) = sqrt(diff' * (MassMatrix * diff));
%     
    
    %     
    %lift = 1/dt (u-uprev,v_l) + nu(grad u,grad v_l) + b(u,u,v_l)
     vt = 1/dt*(velSoln - velPrev);
     lift = -20*( vt'*vlmass' + nu*velSoln'*vlstiff' + velSoln' * NLlift * velSoln );
     drag = -20*( vt'*vdmass' + nu*velSoln'*vdstiff' + velSoln' * NLdrag * velSoln );


     energy = 1/2 * sqrt((velSoln' * (MassROM * velSoln) ));
     %energyDNS = 1.0/2.0 * sqrt((DNSProjectionMatrix(:,ts)' *(MassROM * DNSProjectionMatrix(:,ts)) ));
     %divL2 = sqrt((velSoln' * (GradDivROM * velSoln) ));
     display([ num2str(ts*dt) '  ' num2str(lift) '   '  num2str(drag) '     ' num2str(energy) ])
     dataTableGROM = [dataTableGROM; ts*dt, lift, drag,   energy];
     %dataTableDNS = [dataTableDNS; ts*dt, energyDNS];
        
    solns(:,ts+1)=velSoln;
    velPrevPrev=velPrev;
    velPrev = velSoln;
    
    
    

end
% aveErrorDD = sum(MyL2norm)/length(MyL2norm)

 

% aveErrorG
% aveErrorDD
% figure
% plot(dataTable1(:,1),dataTable1(:,3),'g', dataTable1(:,1),dataTable2(:,3),'r',dataTable1(:,1),dataTable3(:,3),'k', DNStable10(1:numTimeSteps,1)-7, DNStable10(1:numTimeSteps,4),'b-.','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Drag','FontSize',20)
% title(['N=' num2str(N)],'FontSize',20)
% I = legend('G-ROM','iDD-ROM','DD-ROM', 'DNS')
% set(I,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight
% 
% 
% 
% figure
% plot(dataTable1(:,1),dataTable1(:,2),'g',dataTable1(:,1),dataTable2(:,2),'r',dataTable1(:,1),dataTable3(:,2),'k',DNStable10(1:numTimeSteps,1)-7, DNStable10(1:numTimeSteps,5),'b-.','LineWidth',2)
% %plot(dataTable1(:,1),dataTable1(:,2)- BalanceTable1(end-1499:end,5),'k-',dataTable2(:,1),dataTable2(:,2)- BalanceTable2(end-1499:end,5),'r-.','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Lift','FontSize',20)
% title(['N=' num2str(N)],'FontSize',20)
% J = legend('G-ROM','iDD-ROM','DD-ROM', 'DNS')
% set(J,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight

%,



figure
plot(dataTableAB(1:1500,1), dataTableDNS(1:1500,2),'b-.',dataTableGROM(:,1),dataTableGROM(:,2),'g',dataTableAB(:,1),dataTableAB(:,2),'r',dataTableAtilde(:,1),dataTableAtilde(:,2),'k','LineWidth',2)
xlabel('t','FontSize',20)
ylabel('Lift','FontSize',20)
title(['N=' num2str(N)],'FontSize',20)
K = legend('DNS','G-ROM','DDF-ROM-quadratic','DDF-ROM-linear')
set(K,'Interpreter','Latex');
set(gca,'FontSize',20)
axis tight

figure
plot(dataTableAB(1:1500,1), dataTableDNS(1:1500,3),'b-.',dataTableGROM(:,1),dataTableGROM(:,3),'g',dataTableAB(:,1),dataTableAB(:,3),'r',dataTableAtilde(:,1),dataTableAtilde(:,3),'k','LineWidth',2)
xlabel('t','FontSize',20)
ylabel('Drag','FontSize',20)
title(['N=' num2str(N)],'FontSize',20)
K = legend('DNS','G-ROM','DDF-ROM-quadratic','DDF-ROM-linear')
set(K,'Interpreter','Latex');
set(gca,'FontSize',20)
axis tight

figure
plot(dataTableAB(1:1500,1), dataTableDNS(1:1500,4),'b-.',dataTableGROM(:,1),dataTableGROM(:,4),'g',dataTableAB(:,1),dataTableAB(:,4),'r',dataTableAtilde(:,1),dataTableAtilde(:,4),'k','LineWidth',2)
xlabel('t','FontSize',20)
ylabel('Energy','FontSize',20)
title(['N=' num2str(N)],'FontSize',20)
K = legend('DNS','G-ROM','DDF-ROM-quadratic','DDF-ROM-linear')
set(K,'Interpreter','Latex');
set(gca,'FontSize',20)
axis tight

% figure
% plot(dataTable1(:,1),dataTable1(:,4),'g', dataTable1(:,1), dataTableDNS(:,2),'b',dataTable1(:,1),dataTable2(:,4),'r',a,bmm,'k','LineWidth',2)
% %plot(dataTable1(:,1), dataTableDNS(:,2),'b',dataTable1(:,1),dataTable2(:,4),'r',dataTable1(:,1),dataTable3(:,4),'k','LineWidth',2)
% %plot(dataTable1(:,1),dataTable1(:,5)- BalanceTable1(end-1499:end,2),'k-',dataTable2(:,1),dataTable2(:,5)- BalanceTable2(end-1499:end,2),'r-.','LineWidth',2)
% xlabel('t','FontSize',20)
% ylabel('Energy','FontSize',20)
% %title(['N=' num2str(N)],'FontSize',20)
% K = legend('G-ROM','DNS','DDF-ROM + Constraints','DDF-ROM')
% set(K,'Interpreter','Latex');
% set(gca,'FontSize',20)
% axis tight

