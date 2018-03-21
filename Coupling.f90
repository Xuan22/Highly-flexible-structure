!! 03/28/2016 Xuan
module Coupling
	use types
	use str_constants
	use Struct_Dyn
	use ele_index
	use Aero_Dyn
	implicit none
contains
	subroutine CouplingAna(it, flag, omega, tStep, totalF)
		implicit none
		integer,intent(in) :: it
		integer,intent(in) :: flag
		real(dp),intent(in) :: omega
		real(dp),intent(in) :: tStep
		real(dp),dimension(1,Na),intent(out) :: totalF

		integer,parameter :: cycNum = 25
		integer :: cyc
	
		real(dp),dimension(eleN,d*4,Nts,cycNum+1) :: tq2
		real(dp),dimension(eleN,d*4,Nts) :: q2
		real(dp),dimension(eleN,d*4,Nts) :: q2_new
	
		real(dp),dimension(nodN*d,Nts, cycNum+1) :: tqq
		real(dp),dimension(nodN*d,Nts) :: qq
		real(dp),dimension(nodN*d,Nts) :: q_dot
		real(dp),dimension(nodN*d,Nts) :: q_ddot
	
		real(dp),dimension(nodN*d,Nts)::qq_new
		real(dp),dimension(nodN*d,Nts) :: q_dot_new
		real(dp),dimension(nodN*d,Nts) :: q_ddot_new
		!-------------------------------------------------------------------------------------------------------------------
		real(dp),dimension(Nts,cycNum) :: tip

		real(dp),dimension(cycNum,Na) :: time_totalF
		real(dp),dimension(cycNum,Na) :: time_LES
		real(dp),dimension(cycNum,Na) :: time_Wag

		real(dp),dimension(Mmax,Nmax,Nts) :: exterF2
		real(dp),dimension(Mmax,Nmax,Nts) :: exterF3
	
		real(dp),dimension(Mmax+1,Nmax+1,Nts,3) :: deflection
		real(dp),dimension(Mmax+1,Nmax+1,Nts) :: time_deform1
		real(dp),dimension(Mmax+1,Nmax+1,Nts) :: time_deform2
		real(dp),dimension(Mmax+1,Nmax+1,Nts) :: time_deform3
		real(dp),dimension(cycNum) :: tresd
		real(dp):: resd
	
		real(dp),dimension(Mmax+1,Nmax+1) :: xNode
		real(dp),dimension(Mmax+1,Nmax+1) :: yNode
	!!======================================================================================================================
		tip     = 0._dp
		q2      = 0._dp
		qq      = 0._dp
		q_dot   = 0._dp
		q_ddot  = 0._dp

		exterF2 = 0._dp
		exterF3 = 0._dp
		deflection = 0._dp

		tq2(:,:,:,1) = q2
		tqq(:,:,1)   = qq
		! cycles for only structure code
		do cyc = 1, 10
			print *, cyc, '--------'
			call StructDyn(cyc,omega,tStep,exterF2,exterF3,q2,qq,q_dot,q_ddot,deflection,q2_new,qq_new,q_dot_new,q_ddot_new)
			tip(:,cyc) = deflection(1,1,:,3)

			time_Deform1(:,:,:) = deflection(:,:,:,1)
			time_deform2(:,:,:) = deflection(:,:,:,2)
			time_deform3(:,:,:) = deflection(:,:,:,3)

			q2 = q2_new
			qq = qq_new
			q_dot  = q_dot_new
			q_ddot = q_ddot_new

			q_dot  = 0._dp
			q_ddot = 0._dp
			tq2(:,:,:,cyc+1) = q2
			tqq(:,:,cyc+1)   = qq
		
			tresd(cyc) = resd
		end do
		!-------------------------------------------------------------------------------------------------------------------
		!! Iterations
		if (cycNum>10) then
			do cyc = 11, cycNum
				print *, cyc, '--------'
				call StructDyn(cyc,omega,tStep,exterF2,exterF3,q2,qq,q_dot,q_ddot,deflection,q2_new,qq_new,q_dot_new,q_ddot_new)
				tip(:,cyc) = deflection(1,1,:,3)
!				print *,exterF3
				time_Deform1(:,:,:) = deflection(:,:,:,1)
				time_deform2(:,:,:) = deflection(:,:,:,2)
				time_deform3(:,:,:) = deflection(:,:,:,3)

				q2 = q2_new
				qq = qq_new
				q_dot  = q_dot_new
				q_ddot = q_ddot_new

				q_dot  = 0._dp
				q_ddot = 0._dp
				tq2(:,:,:,cyc+1) = q2
				tqq(:,:,cyc+1)   = qq
		
				tresd(cyc)=resd
		
				call AeroDyn(deflection, omega, tStep, totalF, exterF2, exterF3, L_LES, L_Wag)
				time_totalF(cyc,:) = totalF(1,:)
				time_LES(cyc,:) = L_LES(1,:)
				time_Wag(cyc,:) = L_Wag(1,:)
			end do
		end if
!!======================================================================================================================
		if(it==flag) then
		!---------------------------------------------------
			open(unit=100,file='time_deform1.dat')
			do j = 1, Mmax+1
				call node_index(j,coeff,Nx)
				do i = 1, Nx+1
					xNode(j,i) = (i-1)*a
					yNode(j,i) = (j-1)*b
					write(100,*) xNode(j,i), ' ', yNode(j,i), ' ', time_deform3(j,i,1)
				end do
			end do
			close(unit=100)	
		!---------------------------------------------------
			open(unit=200,file='time_deform5.dat')
			do j = 1, Mmax+1
				call node_index(j,coeff,Nx)
				do i = 1, Nx+1
					xNode(j,i) = (i-1)*a
					yNode(j,i) = (j-1)*b
					write(200,*) xNode(j,i), ' ', yNode(j,i), ' ', time_deform3(j,i,5)
				end do
			end do
			close(unit=200)	
		!-----------------------------------------------------------------------------------------------------------------------
			open(unit=500,file='tip1.dat')
			do cyc = 1, cycNum
				do tt = 1, Nts
					write(500,*) (cyc-1)*Nts +tt, '', tip(tt,cyc)
				end do
			end do
			close(unit=500)
		!---------------------------------------------------------------------------------------------------------------
			open(unit=600,file='totalLift.dat')
			do cyc = 1, cycNum
				do i = 1, Nts
					write(600,*) (cyc-1)*Nts+i, time_totalF(cyc,i)
				end do
			enddo
			close(unit=600)
		!---------------------------------------------------------------------------------------------------------------
			open(unit=610,file='finalLift.dat')
			do i = 1, Nts
				write(610,*) i, time_totalF(cycNum,i)
			end do
			close(unit=610)
		!---------------------------------------------------------------------------------------------------------------
			open(unit=620,file='LES.dat')
			do i = 1, Nts
				write(620,*) i, time_LES(cycNum,i)
			end do
			close(unit=620)
		!---------------------------------------------------------------------------------------------------------------
		end if
	end subroutine CouplingAna
end module Coupling



