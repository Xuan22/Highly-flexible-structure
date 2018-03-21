!---------------------------------------------------------------------
! Geometrically Exact FEM Shell Model
! 10/07/2017 modified for Banded Matrix
! 10/08/2017 modified for simplize input for functions
!---------------------------------------------------------------------
module Struct_Dyn
	use types
	use str_constants
	use global_var
	use ele_index
	use GaussQuad_mod
	use shapeFunction
	use shapeFunctionM
	use initCurvature
	use Strain
	use getTransformation
	use modfcn
	use ChangeInCurv
	use curvature_mod
	use Psi_mod
	use Phi_mod
	use massM_mod
	use Derivative
	use getZ_mod
	use getZm_mod
	use getR_mod
	use assembly_mod
	use applyboundary
	use MatInv
	use array_inv
	use addboundary
	implicit none
	real(dp), dimension(Nts+1) :: time
	real(dp), dimension(Mmax+1,Nmax+1) :: w_final
	real(dp), dimension(Mmax+1,Nmax+1) :: u_final
	real(dp), dimension(Mmax+1,Nmax+1) :: v_final

	real(dp), dimension(Nts) :: epsilon1
	real(dp), dimension(Nts) :: residual
	real(dp), dimension(Nts) :: wfc
	
	! getR
	real(dp), dimension(d*4,d*4,eleN) :: EMass
	real(dp), dimension(d*4,d*4,eleN) :: K_tilda
	real(dp), dimension(d*4,1,eleN) :: Kq1
	real(dp), dimension(d*4,1,eleN) :: Rr1
	
	real(dp), dimension(eleN) :: arrE
    real(dp), dimension(eleN) :: arrNu
    real(dp), dimension(eleN) :: arrRho
    real(dp), dimension(eleN) :: arrThck
    real(dp), dimension(eleN) :: arrDm
    real(dp), dimension(eleN,6) :: arrR1
	real(dp), dimension(eleN,6) :: r1
!!======================================================================================================================
contains
	subroutine StructDyn(cyc, omega, tStep, exterF2, exterF3, q2, qq, q_dot, q_ddot, deflection, q2_new, qq_new, q_dot_new, q_ddot_new)
		implicit none
		integer, intent(in) :: cyc
		real(dp), intent(in) :: omega
		real(dp), intent(in) :: tStep
		real(dp), dimension(Mmax,Nmax,Nts), intent(in) :: exterF2
		real(dp), dimension(Mmax,Nmax,Nts), intent(in) :: exterF3
		real(dp), dimension(eleN,d*4,Nts), intent(in) :: q2
		real(dp), dimension(nodN*d,Nts), intent(in) :: qq
		real(dp), dimension(nodN*d,Nts), intent(in) :: q_dot
		real(dp), dimension(nodN*d,Nts), intent(in) :: q_ddot

		real(dp), dimension(Mmax+1,Nmax+1,Nts,3), intent(out) :: deflection
		real(dp), dimension(eleN,d*4,Nts), intent(out) :: q2_new
		real(dp), dimension(nodN*d,Nts), intent(out) :: qq_new
		real(dp), dimension(nodN*d,Nts), intent(out) :: q_dot_new
		real(dp), dimension(nodN*d,Nts), intent(out) :: q_ddot_new
		
		integer, dimension(Mmax) :: ele1Col
		real(dp) :: aeroF
		integer  :: gdof
		integer :: indexDof
		integer :: indexU
		integer :: indexV
		integer :: indexW
		
		! Initialization
		deflection = 0._dp
		q2_new = 0._dp
		qq_new = 0._dp
		q_dot_new = 0._dp
		q_ddot_new = 0._dp
		
		u_final = 0._dp
		v_final = 0._dp
		w_final = 0._dp
!!----------------------------------------------------------------------------------------------------------------------
		! Element Conectivity (nel, index1)
		call sub_ele_index()
!		print *, 'eleN, nodN ', eleN, nodN
		
		!! Gauss Points and Weights
		call GaussQuad(gp, a0, a, x, Wx)
		call GaussQuad(gp, a0, b, y, Wy)
!	!	!===============================================================================================================
		! Start of Time Loop
		do tt = 1, Nts
!			print *, tt
			time(tt) = t_init +(tt-1)*tStep
		    EMass   = 0._dp
		    K_tilda = 0._dp
		    Kq1     = 0._dp                                               				! [K]{q}
		    Rr1     = 0._dp
		    ! pressure load
		    do j = 1, Mmax
			call element_index(j, coeff, Nx)
		    do i = 1, Nx
				ele1Col(1) = 1
				ele = 0
				do ni = 1, Mmax-1
					ele = ele + Nh(ni)*coeff(ni)
					ele1Col(ni+1) = ele+1;
				end do
				ele = 0
				do ni = 1, Mmax-1
					ele = ele +Nh(ni)*coeff(ni)
				end do
				ele = ele +i
!		        if (ele>=17 .AND. ele<=21) then      				! leading edge beam
!		            arrE(ele)    = Em2
!		           	arrNu(ele)   = nu2
!		           	arrRho(ele)  = rho2
!		           	arrThck(ele) = thick2
!		        else
		           	arrE(ele)    = Em1
		           	arrNu(ele)   = nu1
		           	arrRho(ele)  = rho1
		           	arrThck(ele) = thick1
!		        endif
		        if (ele==ele1Col(j)) then							! root beam
		           	arrE(ele)    = Em2
		           	arrNu(ele)   = nu2
		           	arrRho(ele)  = rho2
		           	arrThck(ele) = thick2
				endif
		        dm = arrRho(ele) *arrThck(ele) *a *b
				!-----------------------------------------------------------------------------------------
				r22 = 0._dp
		        cc1 = (pi/180.*55.) *dm /(a*b)
				if (cyc<=5) then
				    r33 =-cc1 *omega**2 *sin(omega*time(tt)-pi/2) *(i-0.5)*a *0.4
				elseif (cyc>5 .and.cyc<=10) then
			        r33 =-cc1 *omega**2 *sin(omega*time(tt)-pi/2) *(i-0.5)*a
		        elseif (cyc>10 .and.cyc<=15) then
				    r33 =-cc1 *omega**2 *sin(omega*time(tt)-pi/2) *(i-0.5)*a +exterF3(j,i,tt) *0.4
				else
				    r33 =-cc1 *omega**2 *sin(omega*time(tt)-pi/2) *(i-0.5)*a +exterF3(j,i,tt)
		        endif
!				r33 = 10._dp*sin(0.087*tt)
				r1(ele, :) = (/0._dp, r22, r33, 0._dp, 0._dp, 0._dp/)
				arrR1(ele,:) = r1(ele,:)
			end do
			end do
!	!	!	!===========================================================================================================
			do ele = 1, eleN
				E     = arrE(ele)
				nu    = arrNu(ele)
				rho   = arrRho(ele)
				thick = arrThck(ele)
				dm    = rho *thick *a *b
				r1(ele,:) = arrR1(ele,:)
		        !-------------------------------------------------------------------------------------------------------
				do jy = 1, size(y)
				do ix = 1, size(x)
		            call shapeFunc2(x(ix), y(jy), Ns, D_mat)                !! Ns is not used later
		            call shapeFuncM(x(ix), y(jy), Nm, Dm_mat)
		            do k = 1, d*4
		            	temp_q2(1, k) = q2(ele, k, tt)
		            end do
		            U_mat = matmul(D_mat, transpose(temp_q2))				    !! (24, 56) * (56, 1) ??
		            call initCurv(T_0, Kcap_10, Kcap_20)
		            call getStrain(U_mat, e1, e2)
		            call getTransform(U_mat, e1, e2, T1, gamma1)
		            call getTransformM(U_mat, T1m, gammam)
		            
		            nf    = 2
		            xf(1) = 0.0001_dp
					xf(2) = 0.0001_dp
					tol   = 1.0e-10
					lwa   = 30
		            call hybrd1(fcn, nf, xf, fvec, tol, info, wa, lwa)         ! g = fsolve(@(g) myfun(g,kg,gamma1), 0);
		            gam2d  = xf(2)
		            gam1d  = gamma1 -gam2d
		            Gam = reshape( (/cos(gam1d), sin(gam2d), 0._dp, &
		                        	 sin(gam1d), cos(gam2d), 0._dp, &
		                        	 0._dp,      0._dp,      1._dp/), (/3, 3/) )
		            call inv22(Gam,idenMat3,Gamma0)                		       ! Gamma = Gam \ eye(3);     !! inverse of matrix
		            Tr = matmul(Gamma0, T1)
					
		            call getChangeInCurv(x(ix), y(jy), T_0, q2(ele,:,tt), Tdx, Tdy)
		            call Curvature(Tr, Tdx, Tdy, Kcap_1, Kcap_2)
					
		            call getPsi(e1, e2, gam1d, gam2d, U_mat, Kcap_1, Kcap_2, psi)
		            call getPhi(E, nu, thick, Kcap_1, Kcap_2, phi)		                		                 		              
		            call getMassM(thick, rho, MassM)
					
		            call getDerivative(x(ix), y(jy), q2(ele,:,tt), T_0, C1sx, C1sy, C2sx, C2sy, C3sx, C3sy, C4sx, C4sy)
		            call getZ(e1, e2, gam1d, gam2d, gamma1, T1, Tdx, Tdy, Tr, T_0, Kcap_1, Kcap_2, &
							C1sx, C1sy, C2sx, C2sy, C3sx, C3sy, C4sx, C4sy, q2(ele,:,tt), x(ix), y(jy), Z_mat)
        			call getZm(T1m, U_mat, Zm)
		            call getR(r1(ele,:), gam1d, gam2d, gamma1, e1, e2, T1, Tr, R1_cap)
					
		            tmp1 = matmul(transpose(Dm_mat), transpose(Zm))
		            tmp2 = matmul(matmul(MassM, Zm), Dm_mat)
		            EMass(:,:,ele) = EMass(:,:,ele) + matmul(tmp1, tmp2) *Wx(ix) *Wy(jy)
					
		            tmp3 = matmul(transpose(D_mat), transpose(Z_mat))
		            tmp4 = matmul(matmul(Phi, Z_mat), D_mat)
		            K_tilda(:,:,ele) = K_tilda(:,:,ele) + matmul(tmp3, tmp4) *Wx(ix) *Wy(jy)
		            
		            tmp5 = matmul(transpose(D_mat), transpose(Z_mat))
		            tmp6 = matmul(Phi, Psi)
		            Kq1(:,:,ele) = Kq1(:,:,ele) + matmul(tmp5, tmp6) *Wx(ix) *Wy(jy)
		            
		            tmp7 = matmul(transpose(D_mat), R1_cap)
		            Rr1(:,:,ele) = Rr1(:,:,ele) + tmp7 *Wx(ix) *Wy(jy)
				end do
				end do
			end do
!	!   !   !===========================================================================================================	
			!! Assembly
			MT = 0._dp                                                                   
			KT = 0._dp
			KQ_cap = 0._dp
			RR_cap = 0._dp
			do ele = 1, eleN
			    call assembly(MT, EMass(:,:,ele),    index1(ele,:), MT) 
			    call assembly(KT, K_tilda(:,:,ele),  index1(ele,:), KT)
			    call assemblyr(KQ_cap, Kq1(:,1,ele), index1(ele,:), KQ_cap)
			    call assemblyr(RR_cap, Rr1(:,1,ele), index1(ele,:), RR_cap)
			end do
!-----------------------------------------------------------------------------------------------------------			
			A1_mat = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14/)
			A2_mat = (/1, 5, 9/)
			call applyboundfixpin(MT, A1_mat, A2_mat, newMT1)
			call applyboundfixpin(KT, A1_mat, A2_mat, newKT1)
			call applyboundfixpinr1(RR_cap, A1_mat, A2_mat, newRR1)
			call applyboundfixpinr1(KQ_cap, A1_mat, A2_mat, newKQ1)
			call applyboundfixpinr1(q_dot(:,tt),  A1_mat, A2_mat, q_dot1)
			call applyboundfixpinr1(q_ddot(:,tt), A1_mat, A2_mat, q_ddot1)
!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
			RR2 = newRR1 -newKQ1 -matmul(newMT1, q_ddot1) +matmul(newMT1, (1._dp/(2._dp*beta))*q_ddot1 +(1._dp/(beta*dt))*q_dot1)
			KT2 = newKT1 + (1._dp/(beta*dt**2))*newMT1
			
			! dQ_mat = inv(KT2) *RR2
			! Solver Inverse of Banded Stiffness Matrix
			do i = 1, NEQ
				KTT(i,NBW) = RR2(i,1)
				do j = 1, NEQ
					if(abs(i-j).LT.NHBW) then
						JJ = j +NHBW -i
						KTT(i,JJ) = KT2(i,j)
					end if
				end do
			end do
			call SLVUNSYM(KTT, NEQ, NEQ, NEQ, NHBW)
			do i = 1,NEQ										!JJ is matrix row
				dQ_mat(i,1) = KTT(i,NBW)
			end do
			
			dQ_dot  = alpha/(beta*dt) *dQ_mat
			dQ_ddot = (1._dp/(beta*dt**2)) *dQ_mat
			!-----------------------------------------------------------------------------------------------------------
			! Adding Deleted Boundary Terms to All Global Node DOF
			call addboundfixpin(dQ_mat, A2_mat, Qr_mat)
			call addboundfixpin(dQ_dot, A2_mat, Qr_dot)
			call addboundfixpin(dQ_ddot, A2_mat, Qr_ddot)
			qq_new(:,tt)     = qq(:,tt) + Qr_mat(:,1)
			q_dot_new(:,tt)  = q_dot(:,tt) + Qr_dot(:,1)
			q_ddot_new(:,tt) = q_ddot(:,tt) + Qr_ddot(:,1)
			!-----------------------------------------------------------------------------------------------------------
			! Find Corresponding Element Nodes DOF Value q2			
			do j = 1, Mmax
				call element_index(j, coeff, Nx)
				do i = 1, Nx
					ele = 0
					do ni = 1, Mmax-1
						ele = ele + Nh(ni) *coeff(ni)
					end do
					ele = ele +i
					do iter = 1, 2
					do k = 1, 14
						gdof = 0
						do ni = 1, Mmax-1
					    	gdof = gdof +(Nh(ni)+1) *d *coeff(ni)
					    end do
					    gdof = gdof +(i-1+iter-1)*d +k
						q2_new(ele, k+d*(iter-1), tt) = qq_new(gdof, tt);
				   	end do
					end do
					do iter = 3, 4
					do k = 1, 14
						gdof = 0
						do ni = 1, Mmax-1
					    	gdof = gdof + (Nh(ni)+1) *d *coeff(ni)
					    end do
					    gdof = gdof +(Nx+1)*d +(i-1+iter-3)*d +k
						q2_new(ele, k+d*(iter-1), tt) = qq_new(gdof, tt);
					end do
					end do
				end do !end i loop
			end do     !end j loop
			!-----------------------------------------------------------------------------------------------------------
			! Return X ,Y, Z Deflection
			do j = 1, Mmax+1
				call node_index(j, coeff, Nx)
				nNode = Nx +1
				do i = 1, nNode
					indexDof = 0
					do ni = 1, Mmax
						indexDof = indexDof + (Nh(ni)+1) *d *coeff(ni) 
					end do
					indexU = indexDof +(i-1)*d +1
					indexV = indexDof +(i-1)*d +5
					indexW = indexDof +(i-1)*d +9
					u_final(j,i) = qq_new(indexU,tt)
					v_final(j,i) = qq_new(indexV,tt)
					w_final(j,i) = qq_new(indexW,tt)
				end do !end i loop
			end do 	   !end j loop
	        
            deflection(:,:,tt,1) = u_final
            deflection(:,:,tt,2) = v_final
            deflection(:,:,tt,3) = w_final
		!---------------------------------------------------------------------------------------------------------------
		end do !end tt
!!======================================================================================================================	
	end subroutine StructDyn
end module Struct_Dyn



