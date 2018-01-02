%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Author: Dr. Leo G. Rebholz. Edited by Jaman Mohebujjaman    %
%  Date: 7/30/2017                                            %
%                                                              %
%  CREATER ROMDRIVER FOR 35K DOF N=16                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the snapshot data from DNS

load snapshotData35Kdt002SV_Re100.mat


% get the ROM basis - cutoff evals that are 1e-10 times smaller than
% biggest eval
cutoff = 1e-16;


%%
%get ROM basis

U = Snapshots(:,1001:end);
S = MassMatrix;

averaging=1;

if averaging==1
    Uave = sum(U')/size(U,2);
    Uave = Uave';
    for i=1:size(U,2);
         U(:,i) = U(:,i) - Uave;
    end
end


M = U' * S * U;

[V,D] = eig(M);

% find the cutoff
d = abs(diag(D));
e = d > d(1)*cutoff;
N = sum(e)

% do the svd
AlphaR = V(:,1:N);
D = diag(D);
LambdaR = D(1:N);





if averaging==1

    PhiR = zeros(NVU,N+1);
    PhiR(:,1) = Uave;
    for i=2:N+1
        x = U*AlphaR(:,i-1);
        PhiR(:,i) = x;
        denom = sqrt( x' * MassMatrix * x);
        PhiR(:,i) = PhiR(:,i) / denom;
    end
else
    PhiR = zeros(NVU,N);
    for i=1:N
        x = U*AlphaR(:,i);
        PhiR(:,i) = x;
        denom = sqrt( x' * MassMatrix * x);
        PhiR(:,i) = PhiR(:,i) / denom;
    end
    
end
%%

% % get the ROM matrices
%%
N=16;
% getROMmatrices;

% Generate ROM mass matrix

% MassROM = zeros(N,N);
% StiffROM = zeros(N,N);
% GradDivROM = zeros(N,N);

clear temp

% % create mass, graddiv and stiffness matrices
% for i=1:N
%     for j=1:N
%         for itrg=1:size(elnode,1)
%             % Description of triangle.
%             triag_no = itrg ;
%             cotri(1:3,1) = nodeco(elnode(triag_no, 1:3), 1) ;
%             cotri(1:3,2) = nodeco(elnode(triag_no, 1:3), 2) ;
%             Jmat = [(cotri(2,1) - cotri(1,1)), (cotri(3,1) - cotri(1,1)) ; ...
%                 (cotri(2,2) - cotri(1,2)) , (cotri(3,2) - cotri(1,2)) ] ;
%             detJ = abs(Jmat(1,1)*Jmat(2,2) - Jmat(1,2)*Jmat(2,1));
%             JInv = inv(Jmat) ;
% 
%             [quad_pts, quad_wghts] = feval('quad_75') ;
%                         
%             nqpts = size(quad_pts,1) ;    
%             quad_wghts = detJ * quad_wghts ;
%             
%                         
%             basvals = basisValues(:,:,itrg);
%             Gradtrue1 = gradBasisValues(:,:,:,itrg);
%             Gradx(:,:) = Gradtrue1(:,1,:) ;
%             Grady(:,:) = Gradtrue1(:,2,:) ;
% 
%             Vstart = [2*(elnode(triag_no,1:3) - 1) + 1 , 2*(elnode(triag_no,4:6) + nVert - 1) + 1 ] ;
%             Vel1i = PhiR([Vstart], i) ;
%             Vel1j = PhiR([Vstart], j) ;
% 
%             
%             Vstart = Vstart + 1 ;
%             Vel2i = PhiR([Vstart], i) ;   
%             Vel2j = PhiR([Vstart], j) ;   
% 
%             vi1 = Vel1i.' * basvals ;
%             vi2 = Vel2i.' * basvals ;
%             
%             vj1 = Vel1j.' * basvals ;
%             vj2 = Vel2j.' * basvals ;
% 
%             
%             MassROM(i,j)=MassROM(i,j) + quad_wghts * (  vi1.*vj1 + vi2.*vj2 ).';
%             
% 
%             localmati(1,1,:) = Vel1i.' * Gradx ;
%             localmati(1,2,:) = Vel2i.' * Gradx ;
%             localmati(2,1,:) = Vel1i.' * Grady ;
%             localmati(2,2,:) = Vel2i.' * Grady ;
% 
%             localmatj(1,1,:) = Vel1j.' * Gradx ;
%             localmatj(1,2,:) = Vel2j.' * Gradx ;
%             localmatj(2,1,:) = Vel1j.' * Grady ;
%             localmatj(2,2,:) = Vel2j.' * Grady ;
%             
%                 temp(1:nqpts) = localmati(1,1,:).*localmatj(1,1,:)  + ...
%                                 localmati(1,2,:).*localmatj(1,2,:)  + ...
%                                 localmati(2,1,:).*localmatj(2,1,:)  + ...
%                                 localmati(2,2,:).*localmatj(2,2,:); 
%                 StiffROM(i,j) = StiffROM(i,j) + quad_wghts * temp.'  ;
%             
%                 temp(1:nqpts) = ( localmati(1,1,:)+localmati(2,2,:) ) .* ( localmatj(1,1,:)+localmatj(2,2,:) );   
%  
%                 GradDivROM(i,j) = GradDivROM(i,j) + quad_wghts * temp.'  ;
%             
%         end
%     end
% end

tic
MassROM = zeros(N,N);
StiffROM = zeros(N,N);
GradDivROM = zeros(N,N);
% create mass, graddiv and stiffness matrices
for i=1:N
    for j=1:N
        MassROM(i,j) = PhiR(:,i)' * (MassMatrix * PhiR(:,j) );
        StiffROM(i,j) = PhiR(:,i)' * (StiffnessMatrix * PhiR(:,j) );
        GradDivROM(i,j) = PhiR(:,i)' * (GradDivMatrix * PhiR(:,j) );
    end
end
toc

% 
% create 3D matrix for the b(u,v,w) term
TriLinROM = zeros(N,N,N);
TriLinROM2 = zeros(N,N,N);


for i=1:N
    i
    for k=1:N
        k
        for j=1:N
          for itrg=1:size(elnode,1)
            % Description of triangle.
            triag_no = itrg ;
            cotri(1:3,1) = nodeco(elnode(triag_no, 1:3), 1) ;
            cotri(1:3,2) = nodeco(elnode(triag_no, 1:3), 2) ;
            Jmat = [(cotri(2,1) - cotri(1,1)), (cotri(3,1) - cotri(1,1)) ; ...
                (cotri(2,2) - cotri(1,2)) , (cotri(3,2) - cotri(1,2)) ] ;
            detJ = abs(Jmat(1,1)*Jmat(2,2) - Jmat(1,2)*Jmat(2,1));
            JInv = inv(Jmat) ;

            [quad_pts, quad_wghts] = feval('quad_75') ;
                        
            nqpts = size(quad_pts,1) ;    
            quad_wghts = detJ * quad_wghts ;
            
            basvals = basisValues(:,:,itrg);
            Gradtrue1 = gradBasisValues(:,:,:,itrg);
            
            Gradx(:,:) = Gradtrue1(:,1,:) ;
            Grady(:,:) = Gradtrue1(:,2,:) ;

            Vstart = [2*(elnode(triag_no,1:3) - 1) + 1 , 2*(elnode(triag_no,4:6) + nVert - 1) + 1 ] ;
            Vel1i = PhiR([Vstart], i) ;
            Vel1k = PhiR([Vstart], k) ;
            Vel1j = PhiR([Vstart], j) ;
                        
            Vstart = Vstart + 1 ;
            Vel2i = PhiR([Vstart], i) ;   
            Vel2k = PhiR([Vstart], k) ;   
            Vel2j = PhiR([Vstart], j) ;   

            vi1 = Vel1i.' * basvals ;
            vi2 = Vel2i.' * basvals ;
            vj1 = Vel1j.' * basvals ;
            vj2 = Vel2j.' * basvals ;
            vk1 = Vel1k.' * basvals ;
            vk2 = Vel2k.' * basvals ;
            
            vj1x = Vel1j.' * Gradx;
            vj1y = Vel1j.' * Grady;
            vj2x = Vel2j.' * Gradx;
            vj2y = Vel2j.' * Grady;

            vi1x = Vel1i.' * Gradx;
            vi1y = Vel1i.' * Grady;
            vi2x = Vel2i.' * Gradx;
            vi2y = Vel2i.' * Grady;

            temp = vi1 .* ( vk1 .* vj1x + vk2 .* vj1y) + vi2 .* (vk1 .* vj2x + vk2 .* vj2y);
            TriLinROM(i,j,k)=TriLinROM(i,j,k)+ quad_wghts * temp.';
            temp2 = vj1 .* ( vk1 .* vi1x + vk2 .* vi1y) + vj2 .* (vk1 .* vi2x + vk2 .* vi2y);
            TriLinROM2(i,j,k)=TriLinROM2(i,j,k)+ quad_wghts * (0.5 * temp - 0.5 * temp2 ).';
            

          end
        end

    end
end



% get lift/drag (phi_i,v_d) etc.
for i=1:N
    vdmass(i) = PhiR(:,i)' * (MassMatrix * CDvec);   %(phi_i,v_d)
    vlmass(i) = PhiR(:,i)' * (MassMatrix * CLvec);
    vdstiff(i) = PhiR(:,i)' * (StiffnessMatrix * CDvec);  %(grad phi_i,grad v_d)
    vlstiff(i) = PhiR(:,i)' * (StiffnessMatrix * CLvec);
end

% get b(phi_i,phi_j,v_d)

NLdrag = zeros(N,N);
NLlift = zeros(N,N);

for k=1:N
        k
        for j=1:N
          for itrg=1:size(elnode,1)
            % Description of triangle.
            triag_no = itrg ;
            cotri(1:3,1) = nodeco(elnode(triag_no, 1:3), 1) ;
            cotri(1:3,2) = nodeco(elnode(triag_no, 1:3), 2) ;
            Jmat = [(cotri(2,1) - cotri(1,1)), (cotri(3,1) - cotri(1,1)) ; ...
                (cotri(2,2) - cotri(1,2)) , (cotri(3,2) - cotri(1,2)) ] ;
            detJ = abs(Jmat(1,1)*Jmat(2,2) - Jmat(1,2)*Jmat(2,1));
            JInv = inv(Jmat) ;

            [quad_pts, quad_wghts] = feval('quad_75') ;
                        
            nqpts = size(quad_pts,1) ;    
            quad_wghts = detJ * quad_wghts ;
            
            basvals = basisValues(:,:,itrg);
            Gradtrue1 = gradBasisValues(:,:,:,itrg);
            
            Gradx(:,:) = Gradtrue1(:,1,:) ;
            Grady(:,:) = Gradtrue1(:,2,:) ;

            Vstart = [2*(elnode(triag_no,1:3) - 1) + 1 , 2*(elnode(triag_no,4:6) + nVert - 1) + 1 ] ;
            Vel1d = CDvec([Vstart]) ;
            Vel1l = CLvec([Vstart]) ;
            Vel1k = PhiR([Vstart], k) ;
            Vel1j = PhiR([Vstart], j) ;
                        
            Vstart = Vstart + 1 ;
            Vel2d = CDvec([Vstart]) ;   
            Vel2l = CLvec([Vstart]) ;   
            Vel2k = PhiR([Vstart], k) ;   
            Vel2j = PhiR([Vstart], j) ;   

            vl1 = Vel1l.' * basvals ;
            vl2 = Vel2l.' * basvals ;
            vd1 = Vel1d.' * basvals ;
            vd2 = Vel2d.' * basvals ;
            vk1 = Vel1k.' * basvals ;
            vk2 = Vel2k.' * basvals ;
            
            vj1x = Vel1j.' * Gradx;
            vj1y = Vel1j.' * Grady;
            vj2x = Vel2j.' * Gradx;
            vj2y = Vel2j.' * Grady;
            
            temp = vd1 .* ( vk1 .* vj1x + vk2 .* vj1y) + vd2 .* (vk1 .* vj2x + vk2 .* vj2y);
            NLdrag(j,k)=NLdrag(j,k)+ quad_wghts * temp.';

            temp = vl1 .* ( vk1 .* vj1x + vk2 .* vj1y) + vl2 .* (vk1 .* vj2x + vk2 .* vj2y);
            NLlift(j,k)=NLlift(j,k)+ quad_wghts * temp.';


          end
        end

end


save ROMtestSV35K_N16_166

function [quad_pts, quad_wghts] = quad_75
%
% This function contains weights and quadrature points
% which are exact for polynomials of degree five
%

quad_pts(1,:) = [(6 + sqrt(15))/21 , (9 - 2*sqrt(15))/21] ;
quad_pts(2,:) = [(6 + sqrt(15))/21 , (6 + sqrt(15))/21 ] ;
quad_pts(3,:) = [(9 - 2*sqrt(15))/21 , (6 + sqrt(15))/21 ] ;
quad_pts(4,:) = [(6 - sqrt(15))/21 , (6 - sqrt(15))/21 ] ;
quad_pts(5,:) = [(9 + 2*sqrt(15))/21 , (6 - sqrt(15))/21 ] ;
quad_pts(6,:) = [(6 - sqrt(15))/21 , (9 + 2*sqrt(15))/21 ] ;
quad_pts(7,:) = [1/3 , 1/3] ;


quad_wghts(1:3) = (155 + sqrt(15))/2400* ones(1,3) ;
quad_wghts(4:6) = (155 - sqrt(15))/2400* ones(1,3) ;
quad_wghts(7) = 270/2400 ;

end


