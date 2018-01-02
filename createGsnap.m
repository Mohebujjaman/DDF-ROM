clear all
load snapshotData35Kdt002SV_Re100
load ROMtestSV35K_N16_166  Snapshots MassROM StiffROM TriLinROM2 NLlift NLdrag vdmass vdstiff vlmass vlstiff GlobalV PhiR MassMatrix T dt nu BalanceTable nodeco GradDivMatrix elnode
% loop over one period timesteps (snapshots)
tic
n = 166
r = 8;
d = 16


wr = 0*Snapshots;
wd = 0*Snapshots;
Gsnap = zeros(r,n);

for i=1:n
   % display(['creating Gsnap, i=' num2str(i)])
    % Project (ur \cdot \nabla ur) into Xh
    
    % First get ur and ud
    vvv = Snapshots(:,1000+i);
    RHS = zeros(d,1);
    for j=1:d
        RHS(j) = vvv' * (MassMatrix * PhiR(:,j) );
    end
    A = MassROM(1:d,1:d);
    RHS = RHS - A(:,1)*1;
    A(1,:)=0;
    A(:,1)=0;
    A(1,1)=1;
    RHS(1)=1;
    ud = A  \ RHS;
    
    ur = ud(1:r);
    
    %% Now do the projection to get wr = P_Xh ( ur dot grad ur)
    
    % First put ur back in the FE basis
    Solnr = 0*GlobalV;
    for j=1:r
      Solnr = Solnr + ur(j)*PhiR(:,j);
    end
 %   GlobalV = Solnr;
    
    Solnd = 0*GlobalV;
    for j=1:d
      Solnd = Solnd + ud(j)*PhiR(:,j);
    end
 %   AllV(:,end)=Solnd;
    

    % Now do the projection
    RHSvec1 = zeros(NVU,1) ;
    RHSvec2 = zeros(NVU,1) ;
    Acoeff = MassMatrix;
    
    for itrg=1:size(elnode,1)

        % Set up unknown solution vector mapping to global velocity vector
          if strcmp(vel_bas_type, 'CtsQuad') == 1
                Vstart = [2*(elnode(itrg,1:3) - 1) + 1 , 2*(elnode(itrg,4:6) + nVert - 1) + 1];
                GlTrgVe = reshape([Vstart ; Vstart+1],[12,1]) ;            
          end
          localUnk=GlTrgVe;

          % velFun
          npts = size(quad_pts,1);
          ten1a=basisValues(:,:,itrg);         
          Gradtrue = gradBasisValues(:,:,:,itrg);
          Gradx(:,:) = Gradtrue(:,1,:) ;
          Grady(:,:) = Gradtrue(:,2,:) ;
          quad_wghts = quad_wghtsglobal(itrg,:);
          ten1=0*ten1a;
          for iq = 1:npts
            ten1(:,iq) = quad_wghts(iq) * ten1a(:,iq) ;  
          end
          nbas1 = size(ten1,1) ;
          
          Vstart = [2*(elnode(itrg,1:3) - 1) + 1 , 2*(elnode(itrg,4:6) + nVert - 1) + 1 ] ;
          Vel1 = Solnr(Vstart) ;
          Vstart = Vstart + 1 ;
          Vel2 = Solnr(Vstart) ;
   
          Velfun1(1,:) = Vel1.' * ten1a ;
          Velfun1(2,:) = Vel2.' * ten1a ;
          
          avec1x(1:npts) = Vel1.' * Gradx ;
          avec1y(1:npts) = Vel1.' * Grady ;
          avec2x(1:npts) = Vel2.' * Gradx ;
          avec2y(1:npts) = Vel2.' * Grady ;

          vv1(1,:) = Velfun1(1,:) .* avec1x + Velfun1(2,:) .* avec1y ;
          vv1(2,:) = Velfun1(1,:) .* avec2x + Velfun1(2,:) .* avec2y ;
          
          mat1 = ten1 * vv1(1,:).' ;
          mat2 = ten1 * vv1(2,:).' ;

          localmat = [ mat1 ; mat2 ] ; 

          [rdim cdim] = size(localmat) ;
          rowperm = [ ] ; colperm = [ ] ; 
          for ir = 1:nbas1
            rowperm = [rowperm ir:nbas1:rdim] ;
          end
          rhscvc1 = localmat(rowperm) ;
          
          % next one
          Vstart = [2*(elnode(itrg,1:3) - 1) + 1 , 2*(elnode(itrg,4:6) + nVert - 1) + 1 ] ;
          Vel1 = Solnd(Vstart) ;
          Vstart = Vstart + 1 ;
          Vel2 = Solnd(Vstart) ;
   
          Velfun1(1,:) = Vel1.' * ten1a ;
          Velfun1(2,:) = Vel2.' * ten1a ;
          
          avec1x(1:npts) = Vel1.' * Gradx ;
          avec1y(1:npts) = Vel1.' * Grady ;
          avec2x(1:npts) = Vel2.' * Gradx ;
          avec2y(1:npts) = Vel2.' * Grady ;

          vv1(1,:) = Velfun1(1,:) .* avec1x + Velfun1(2,:) .* avec1y ;
          vv1(2,:) = Velfun1(1,:) .* avec2x + Velfun1(2,:) .* avec2y ;
          
          mat1 = ten1 * vv1(1,:).' ;
          mat2 = ten1 * vv1(2,:).' ;

          localmat = [ mat1 ; mat2 ] ; 
          rhscvc2 = localmat(rowperm) ;
          
          
          
%           
%           
%           rhscvc1b = inner_prod_ten1_Vec(itrg,'quad_75','velGradvel', vel_bas_type) ;
%           rhscvc2b = inner_prod_ten1_Vec(itrg,'quad_75','velPGradvelP', vel_bas_type) ;
%    
%           norm(rhscvc1 - rhscvc1b)
%           norm(rhscvc2 - rhscvc2b)
%           
          RHSvec1(localUnk) = RHSvec1(localUnk) + rhscvc1;
          RHSvec2(localUnk) = RHSvec2(localUnk) + rhscvc2;
 
    end

    % These wr and wd are in Xh
    wr(:,i) = Acoeff \ RHSvec1;
    wd(:,i) = Acoeff \ RHSvec2;
    
    
    % We now need to filter wd with ROM filter (r), and project wr into ROM
    % space
    
    vvv = wd(:,i);
    RHS = zeros(r,1);
    for j=1:r
        RHS(j) = vvv' * (MassMatrix * PhiR(:,j) );
    end
    A = MassROM(1:r,1:r);
%     RHS = RHS - A(:,1)*1;
%     A(1,:)=0;
%     A(:,1)=0;
%     A(1,1)=1;
%     RHS(1)=1;
    wdr_rom(:,i) = A  \ RHS;
   
    vvv = wr(:,i);
    RHS = zeros(r,1);
    for j=1:r
        RHS(j) = vvv' * (MassMatrix * PhiR(:,j) );
    end
    A = MassROM(1:r,1:r);
%     RHS = RHS - A(:,1)*1;
%     A(1,:)=0;
%     A(:,1)=0;
%     A(1,1)=1;
%     RHS(1)=1;
    wr_rom(:,i) = A  \ RHS;
    
    % Now we can get Gsnap
    Gsnap(:,i)= ((wdr_rom(:,i) - wr_rom(:,i))' * MassROM(1:r,1:r))';
    
end
toc
    
save Gsnap_SV35K_r8_d16_N16_166
