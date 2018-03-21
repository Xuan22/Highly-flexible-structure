module ChangeInCurvk_mod
	use types
	use str_constants, only : a, b, k_50, k_10, k_610, k_40, k_620, k_20
	use global_var, only : ele, counter
	use InitCurvature
	use ChangeInT
	use curvature_mod
	implicit none
contains
subroutine ChangeInCurvK(T_0, q2ele, xi, yj, Tdx, Tdy, K1x_cap, K2x_cap, K1y_cap, K2y_cap)
		implicit none
		real(dp), dimension(3,3), intent(in) :: T_0
		real(dp), dimension(1,d*4), intent(in) :: q2ele
	
		real(dp), intent(in) :: xi, yj
		real(dp), dimension(3,3), intent(in) :: Tdx,Tdy
		
		real(dp), dimension(3,3) :: K1x_cap, K2x_cap, K1y_cap, K2y_cap 
		real(dp), dimension(3,3) :: temp1, temp2
		real(dp), dimension(3,3) :: Tpx, Tpy, Tmx, Tmy
		real(dp), dimension(3,3) :: K1px_cap, K2px_cap, K1py_cap, K2py_cap, K1mx_cap, K2mx_cap, K1my_cap, K2my_cap
		real(dp) :: temp3
		real(dp) :: dx, dy

		call initCurv(T_0, temp1, temp2)                   !! what is ~~ ?? //

		dx = a/100._dp
		dy = 0._dp
		call getChangeInT(xi,yj,q2ele,T_0,dx,dy,Tpx)
		call Curvature(Tpx,Tdx,Tdy,K1px_cap,K2px_cap)
		 
		dx = 0._dp
		dy = b /100._dp
		call getChangeInT(xi,yj,q2ele,T_0,dx,dy,Tpy)
		call Curvature(Tpy,Tdx,Tdy,K1py_cap,K2py_cap)
		 
		dx =-a /100._dp
		dy = 0._dp
		call getChangeInT(xi,yj,q2ele,T_0,dx,dy,Tmx)
		call Curvature(Tmx,Tdx,Tdy,K1mx_cap,K2mx_cap)
		 
		dx = 0._dp
		dy =-b /100._dp
		call getChangeInT(xi,yj,q2ele,T_0,dx,dy,Tmy)
		call Curvature(Tmy,Tdx,Tdy,K1my_cap,K2my_cap)

		K1x_cap  = (K1px_cap -K1mx_cap) /(2.*a/100.)
		K2x_cap  = (K2px_cap -K2mx_cap) /(2.*a/100.)

		K1y_cap = (K1py_cap -K1my_cap) /(2.*b/100.)
		K2y_cap = (K2py_cap -K2my_cap) /(2.*b/100.)

end subroutine ChangeInCurvK

end module ChangeIncurvK_mod


