!! Xuan 02/29/2016
module getZ_mod
	use types
	use str_constants, only : a, b, k_50, k_10, k_610, k_40, k_620, k_20
	use global_var, only : ele, counter
	use changeInCurvK_mod
	implicit none
contains
	subroutine getZ(e1d, e2d, gam1d, gam2d, gammad, T1d, Tddx, Tddy, Td, T0d, Kcap_1, Kcap_2, &
					C1sx, C1sy, C2sx, C2sy, C3sx, C3sy, C4sx, C4sy, q2ele, xi, yj, Z_mat)
		implicit none
		real(dp), intent(in) :: e1d, e2d, gam1d, gam2d, gammad
	
		real(dp), dimension(3,3), intent(in) :: T1d, Tddx, Tddy, Td, T0d
		real(dp), dimension(3,3), intent(in) :: Kcap_1, Kcap_2
		real(dp), intent(in) :: C1sx, C1sy, C2sx, C2sy, C3sx, C3sy, C4sx, C4sy
	
		real(dp), dimension(1, d*4), intent(in) :: q2ele
		real(dp), intent(in) :: xi, yj
	
		real(dp), dimension(12,24), intent(out) :: Z_mat                                      !! is dimension of Z right?? //
		
		real(dp) :: k_5, k_1, k_61, k_4, k_62, k_2
		real(dp), dimension(3,3) :: Kcap_10x, Kcap_20x, Kcap_10y, Kcap_20y
		real(dp) :: k_10dx, k_20dx, k_610dx, k_620dx, k_50dx, k_40dx
		real(dp) :: k_10dy, k_20dy, k_610dy, k_620dy, k_50dy, k_40dy
		real(dp) :: C0, C1, C2, C3, C4
		real(dp) :: C01, C02, C03, C04, C05, C06, C11, C12, C13, C14, C15, C16, C21, C22, C23, C24, C25, C26
		real(dp) :: C31, C32, C33, C34, C35, C36, C41, C42, C43, C44, C45, C46, C51, C52, C53, C54, C55, C56
		real(dp) :: C61, C62, C63, C64, C65, C66

		k_5  = Kcap_1(1,2)
		k_1  =-Kcap_1(1,3)
		k_61 =-Kcap_1(2,3)
		k_4  = Kcap_2(1,2)
		k_62 =-Kcap_2(1,3)
		k_2  =-Kcap_2(2,3)

		call changeInCurvK(T0d,q2ele,xi,yj,Tddx,Tddy,Kcap_10x, Kcap_20x, Kcap_10y, Kcap_20y)
		k_50dx  = Kcap_10x(1,2)
		k_10dx  =-Kcap_10x(1,3)
		k_610dx =-Kcap_10x(2,3)
		k_40dx  = Kcap_20x(1,2)
		k_620dx =-Kcap_20x(1,3)
		k_20dx  =-Kcap_20x(2,3)

		k_50dy  = Kcap_10y(1,2)
		k_10dy  =-Kcap_10y(1,3)
		k_610dy =-Kcap_10y(2,3)
		k_40dy  = Kcap_20y(1,2)
		k_620dy =-Kcap_20y(1,3)
		k_20dy  =-Kcap_20y(2,3)

		C0  = (1.+e1d) *cos(gam1d) + (1.+e2d) *cos(gam2d)
		C01 = (T1d(2,1) -sin(gammad)*T1d(1,1)) /cos(gammad)
		C02 = (T1d(2,2) -sin(gammad)*T1d(1,2)) /cos(gammad)
		C03 = (T1d(2,3) -sin(gammad)*T1d(1,3)) /cos(gammad)
		C04 = (T1d(1,1) -sin(gammad)*T1d(2,1)) /cos(gammad)
		C05 = (T1d(1,2) -sin(gammad)*T1d(2,2)) /cos(gammad)
		C06 = (T1d(1,3) -sin(gammad)*T1d(2,3)) /cos(gammad)

		C1 = cos(gam1d) /cos(gammad)*(1.+e2d)
		C2 = sin(gam2d) /cos(gammad)*(1.+e1d)
		C3 = sin(gam1d) /cos(gammad)*(1.+e2d)
		C4 = cos(gam2d) /cos(gammad)*(1.+e1d)

		C11 = (1.+e1d+(1.+e2d)*cos(gam2d)*cos(gam1d)) *T1d(1,1) - (1.+e2d)*sin(gam1d)*cos(gam2d)*C01
		C12 = (1.+e1d+(1.+e2d)*cos(gam2d)*cos(gam1d)) *T1d(1,2) - (1.+e2d)*sin(gam1d)*cos(gam2d)*C02
		C13 = (1.+e1d+(1.+e2d)*cos(gam2d)*cos(gam1d)) *T1d(1,3) - (1.+e2d)*sin(gam1d)*cos(gam2d)*C03
		C14 = (-(1.+e1d)*sin(gam2d)*sin(gam1d)) *T1d(2,1) - (1.+e1d)*sin(gam1d)*cos(gam2d)*C04
		C15 = (-(1.+e1d)*sin(gam2d)*sin(gam1d)) *T1d(2,2) - (1.+e1d)*sin(gam1d)*cos(gam2d)*C05
		C16 = (-(1.+e1d)*sin(gam2d)*sin(gam1d)) *T1d(2,3) - (1.+e1d)*sin(gam1d)*cos(gam2d)*C06

		C21 = (-(1.+e2d)*sin(gam2d)*sin(gam1d)) *T1d(1,1) - (1.+e2d)*sin(gam2d)*cos(gam1d)*C01
		C22 = (-(1.+e2d)*sin(gam2d)*sin(gam1d)) *T1d(1,2) - (1.+e2d)*sin(gam2d)*cos(gam1d)*C02
		C23 = (-(1.+e2d)*sin(gam2d)*sin(gam1d)) *T1d(1,3) - (1.+e2d)*sin(gam2d)*cos(gam1d)*C03
		C24 = (1.+e2d+(1.+e1d)*cos(gam2d)*cos(gam1d)) * T1d(2,1) - (1.+e1d)*sin(gam2d)*cos(gam1d)*C04      !! deleted unnecessary() of T1d()
		C25 = (1.+e2d+(1.+e1d)*cos(gam2d)*cos(gam1d)) * T1d(2,2) - (1.+e1d)*sin(gam2d)*cos(gam1d)*C05
		C26 = (1.+e2d+(1.+e1d)*cos(gam2d)*cos(gam1d)) * T1d(2,3) - (1.+e1d)*sin(gam2d)*cos(gam1d)*C06

		C31 = (1.+e2d)*cos(gam2d)*sin(gam1d)*(T1d(1,1)) + (1.+e2d)*cos(gam2d)*cos(gam1d)*C01
		C32 = (1.+e2d)*cos(gam2d)*sin(gam1d)*(T1d(1,2)) + (1.+e2d)*cos(gam2d)*cos(gam1d)*C02
		C33 = (1.+e2d)*cos(gam2d)*sin(gam1d)*(T1d(1,3)) + (1.+e2d)*cos(gam2d)*cos(gam1d)*C03
		C34 = (1.+e1d)*cos(gam1d)*sin(gam2d)*(T1d(2,1)) + (1.+e1d)*cos(gam2d)*cos(gam1d)*C04
		C35 = (1.+e1d)*cos(gam1d)*sin(gam2d)*(T1d(2,2)) + (1.+e1d)*cos(gam2d)*cos(gam1d)*C05
		C36 = (1.+e1d)*cos(gam1d)*sin(gam2d)*(T1d(2,3)) + (1.+e1d)*cos(gam2d)*cos(gam1d)*C06


		C41 = k_61/(2.*C0) *(2.*sin(gam1d)*T1d(1,1)+(((1.+e1d)*cos(gam1d)-&
							(1.+e2d)*cos(gam2d)+C0))*C01/(1+e1d)) -(C4*Tddx(3,1)+Td(3,1)*C4sx)-k_5*C2*Td(3,1)
		C42 = k_61/(2.*C0) *(2.*sin(gam1d)*T1d(1,2)+(((1.+e1d)*cos(gam1d)-&
							(1.+e2d)*cos(gam2d)+C0))*C02/(1+e1d)) -(C4*Tddx(3,2)+Td(3,2)*C4sx)-k_5*C2*Td(3,2)
		C43 = k_61/(2.*C0) *(2.*sin(gam1d)*T1d(1,3)+(((1.+e1d)*cos(gam1d)-&
							(1.+e2d)*cos(gam2d)+C0))*C03/(1.+e1d)) -(C4*Tddx(3,3)+Td(3,3)*C4sx)-k_5*C2*Td(3,3)
		C44 = k_61/(2.*C0) *(2.*sin(gam2d)*T1d(2,1)+(((1.+e1d)*cos(gam1d)-&
							(1.+e2d)*cos(gam2d)-C0))*C04/(1.+e2d)) -(C3*Tddx(3,1)+Td(3,1)*C3sx)-k_5*C1*Td(3,1)
		C45 = k_61/(2.*C0) *(2.*sin(gam2d)*T1d(2,2)+(((1.+e1d)*cos(gam1d)-&
							(1.+e2d)*cos(gam2d)-C0))*C05/(1+e2d)) -(C3*Tddx(3,2)+Td(3,2)*C3sx)-k_5*C1*Td(3,2)
		C46 = k_61/(2.*C0) *(2.*sin(gam2d)*T1d(2,3)+(((1.+e1d)*cos(gam1d)-&
							(1.+e2d)*cos(gam2d)-C0))*C06/(1+e2d)) -(C3*Tddx(3,3)+Td(3,3)*C3sx)-k_5*C1*Td(3,3)

		C51 = -k_62/(2.*C0) *(2.*sin(gam1d)*T1d(1,1)+(((1.+e1d)*cos(gam1d)-&
							 (1.+e2d)*cos(gam2d)+C0))*C01/(1+e1d)) -(C2*Tddy(3,1)+Td(3,1)*C2sy)-k_4*C4*Td(3,1)
		C52 = -k_62/(2.*C0) *(2.*sin(gam1d)*T1d(1,2)+(((1.+e1d)*cos(gam1d)-&
							(1.+e2d)*cos(gam2d)+C0))*C02/(1+e1d)) -(C2*Tddy(3,2)+Td(3,2)*C2sy)-k_4*C4*Td(3,2)
		C53 = -k_62/(2.*C0) *(2.*sin(gam1d)*T1d(1,3)+(((1.+e1d)*cos(gam1d)-&
							(1.+e2d)*cos(gam2d)+C0))*C03/(1+e1d)) -(C2*Tddy(3,3)+Td(3,3)*C2sy)-k_4*C4*Td(3,3)
		C54 = -k_62/(2.*C0) *(-2.*sin(gam2d)*T1d(2,1)+(((1.+e1d)*cos(gam1d)-&
							(1.+e2d)*cos(gam2d)-C0))*C04/(1+e2d)) -(C1*Tddy(3,1)+Td(3,1)*C1sy)-k_4*C3*Td(3,1)
		C55 = -k_62/(2.*C0) *(-2.*sin(gam2d)*T1d(2,2)+(((1.+e1d)*cos(gam1d)-&
							(1.+e2d)*cos(gam2d)-C0))*C05/(1+e2d)) -(C1*Tddy(3,2)+Td(3,2)*C1sy)-k_4*C3*Td(3,2)
		C56 = -k_62/(2.*C0) *(-2.*sin(gam2d)*T1d(2,3)+(((1.+e1d)*cos(gam1d)-&
							 (1.+e2d)*cos(gam2d)-C0))*C06/(1.+e2d)) -(C1*Tddy(3,3)+Td(3,3)*C1sy)-k_4*C3*Td(3,3)

		C61 = (k_2-k_1)/(2.*C0) *(2.*sin(gam1d)*T1d(1,1)+(((1.+e1d)*cos(gam1d) -(1.+e2d)*cos(gam2d)+C0))*C01/(1.+e1d))+&
								 (C2*Tddx(3,1)+Td(3,1)*C2sx) -(C4*Tddy(3,1)+Td(3,1)*C4sy) -k_4*C2*Td(3,1)-k_5*C4*Td(3,1)
		C62 = (k_2-k_1)/(2.*C0) *(2.*sin(gam1d)*T1d(1,2)+(((1.+e1d)*cos(gam1d) -(1.+e2d)*cos(gam2d)+C0))*C02/(1.+e1d)) +&
								 (C2*Tddx(3,2)+Td(3,2)*C2sx) -(C4*Tddy(3,2)+Td(3,2)*C4sy) -k_4*C2*Td(3,2)-k_5*C4*Td(3,2)
		C63 = (k_2-k_1)/(2.*C0) *(2.*sin(gam1d)*T1d(1,3)+(((1.+e1d)*cos(gam1d) -(1.+e2d)*cos(gam2d)+C0))*C03/(1.+e1d)) +&
								 (C2*Tddx(3,3)+Td(3,3)*C2sx) -(C4*Tddy(3,3)+Td(3,3)*C4sy) -k_4*C2*Td(3,3)-k_5*C4*Td(3,3)

		C64 = (k_2-k_1)/(2.*C0) *(-2.*sin(gam2d)*T1d(2,1)+(((1.+e1d)*cos(gam1d) -(1.+e2d)*cos(gam2d)-C0))*C04/(1.+e2d)) +&
								 (C3*Tddy(3,1)+Td(3,1)*C3sy) -(C1*Tddx(3,1)+Td(3,1)*C1sx) -k_4*C1*Td(3,1)-k_5*C3*Td(3,1)
		C65 = (k_2-k_1)/(2.*C0) *(-2.*sin(gam2d)*T1d(2,2)+(((1.+e1d)*cos(gam1d) -(1.+e2d)*cos(gam2d)-C0))*C05/(1.+e2d)) +&
								 (C3*Tddy(3,2)+Td(3,1)*C3sy) -(C1*Tddx(3,2)+Td(3,2)*C1sx) -k_4*C1*Td(3,2)-k_5*C3*Td(3,2)
		C66 = (k_2-k_1)/(2.*C0) *(-2.*sin(gam2d)*T1d(2,3)+(((1.+e1d)*cos(gam1d) -(1.+e2d)*cos(gam2d)-C0))*C06/(1.+e2d)) +&
								 (C3*Tddy(3,3)+Td(3,3)*C3sy) -(C1*Tddx(3,3)+Td(3,3)*C1sx) -k_4*C1*Td(3,3)-k_5*C3*Td(3,3)

		Z_mat(1,2)  = C11/C0
		Z_mat(1,8)  = C12/C0
		Z_mat(1,14) = C13/C0
		Z_mat(1,3)  = C14/C0
		Z_mat(1,9)  = C15/C0
		Z_mat(1,15) = C16/C0
		Z_mat(1,1)  = (C12*k_50-C13*k_10+C15*k_40-C16*k_620)/C0
		Z_mat(1,7)  =-(C11*k_50+C13*k_610+C14*k_40+C16*k_20)/C0
		Z_mat(1,13) = (C11*k_10+C12*k_610+C14*k_620+C15*k_20)/C0

		Z_mat(2,2)  = C21/C0
		Z_mat(2,8)  = C22/C0
		Z_mat(2,14) = C23/C0
		Z_mat(2,3)  = C24/C0
		Z_mat(2,9)  = C25/C0
		Z_mat(2,15) = C26/C0
		Z_mat(2,1)  = (C22*k_50-C23*k_10+C25*k_40-C26*k_620) / C0
		Z_mat(2,7)  =-(C21*k_50+C23*k_610+C24*k_40+C26*k_20) / C0
		Z_mat(2,13) = (C21*k_10+C22*k_610+C24*k_620+C25*k_20) / C0

		Z_mat(3,2)  = 2.*C31/C0
		Z_mat(3,8)  = 2.*C32/C0
		Z_mat(3,14) = 2.*C33/C0
		Z_mat(3,3)  = 2.*C34/C0
		Z_mat(3,9)  = 2.*C35/C0
		Z_mat(3,15) = 2.*C36/C0
		Z_mat(3,1)  = 2.*(C32*k_50-C33*k_10+C35*k_40-C36*k_620)/C0
		Z_mat(3,7)  =-2.*(C31*k_50+C33*k_610+C34*k_40+C36*k_20)/C0
		Z_mat(3,13) = 2.*(C31*k_10+C32*k_610+C34*k_620+C35*k_20)/C0

		Z_mat(4,2)  = C41+C3*Td(3,2)*k_40-C3*Td(3,3)*k_620-C4*Td(3,2)*k_50+C4*Td(3,3)*k_10
		Z_mat(4,8)  = C42-C3*Td(3,1)*k_40-C3*Td(3,3)*k_20+C4*Td(3,1)*k_50+C4*Td(3,3)*k_610
		Z_mat(4,14) = C43+C3*Td(3,1)*k_620+C3*Td(3,2)*k_20-C4*Td(3,1)*k_10-C4*Td(3,2)*k_610
		Z_mat(4,3)  = C44
		Z_mat(4,9)  = C45 
		Z_mat(4,15) = C46
		Z_mat(4,1)  = C42*k_50-C43*k_10+C45*k_40-C46*k_20+C3*Td(3,2)*k_40dx-C3*Td(3,3)*k_620dx-C4*Td(3,2)*k_50dx+C4*Td(3,3)*k_10dx
		Z_mat(4,7)  =-C41*k_50-C43*k_610-C44*k_40-C46*k_20-C3*Td(3,1)*k_40dx-C3*Td(3,3)*k_20dx+C4*Td(3,1)*k_50dx+C4*Td(3,3)*k_610dx
		Z_mat(4,13) = C41*k_10+C42*k_610+C44*k_620+C45*k_20-C3*Td(3,1)*k_620dx+C3*Td(3,2)*k_20dx-C4*Td(3,1)*k_10dx-C4*Td(3,2)*k_610dx

		Z_mat(4,5)  = C3*Td(3,1)
		Z_mat(4,11) = C3*Td(3,2)
		Z_mat(4,17) = C3*Td(3,3)
		Z_mat(4,4)  =-C4*Td(3,1)
		Z_mat(4,16) =-C4*Td(3,3)
		Z_mat(5,2)  = C51
		Z_mat(5,8)  = C52
		Z_mat(5,14) = C53
		Z_mat(5,3)  = C54 -C1*Td(3,2)*k_40 +C1*Td(3,3)*k_620 +C2*Td(3,2)*k_50 -C2*Td(3,3)*k_10
		Z_mat(5,9)  = C55 +C1*Td(3,1)*k_40 +C1*Td(3,3)*k_20 -C2*Td(3,1)*k_50 -C2*Td(3,3)*k_610
		Z_mat(5,15) = C56 -C1*Td(3,1)*k_620 -C1*Td(3,2)*k_20 +C2*Td(3,1)*k_10 +C2*Td(3,2)*k_610

		Z_mat(5,1) = C52*k_50 -C53*k_10 +C55*k_40 -C56*k_620 -C1*Td(3,2)*k_40dy +C1*Td(3,3)*k_620dy +&
				     C2*Td(3,2)*k_50dy -C2*Td(3,3)*k_10dx
		Z_mat(5,7) =-C51*k_50 -C53*k_610 -C54*k_40 -C56*k_20 +C1*Td(3,1)*k_40dy +C1*Td(3,3)*k_20dy -&
				     C2*Td(3,1)*k_50dy -C2*Td(3,3)*k_610dy
		Z_mat(5,13) = C51*k_10 +C52*k_610 +C54*k_620 +C55*k_20 -C1*Td(3,1)*k_620dy -C1*Td(3,2)*k_620dy +&
				      C2*Td(3,1)*k_10dy +C2*Td(3,2)*k_610dy

		Z_mat(5,5)  = C2*Td(3,1)
		Z_mat(5,11) = C2*Td(3,1)
		Z_mat(5,17) = C2*Td(3,3)
		Z_mat(5,6)  =-C1*Td(3,1)
		Z_mat(5,12) =-C1*Td(3,2)
		Z_mat(5,18) =-C1*Td(3,3)

		Z_mat(6,2)  = C61 -C1*Td(3,2)*k_40 +C1*Td(3,3)*k_620 +C2*Td(3,2)*k_50 -C2*Td(3,3)*k_10
		Z_mat(6,8)  = C62 +C1*Td(3,1)*k_40 +C1*Td(3,3)*k_20 -C2*Td(3,1)*k_50 -C2*Td(3,3)*k_610
		Z_mat(6,14) = C63 -C1*Td(3,1)*k_620 -C1*Td(3,2)*k_20 +C2*Td(3,1)*k_10 +C2*Td(3,2)*k_610
		Z_mat(6,3)  = C64 +C3*Td(3,2)*k_40 -C3*Td(3,3)*k_620 -C4*Td(3,2)*k_50 +C4*Td(3,3)*k_10
		Z_mat(6,9)  = C65 -C3*Td(3,1)*k_40 -C3*Td(3,3)*k_20 +C4*Td(3,1)*k_50 +C4*Td(3,3)*k_610
		Z_mat(6,15) = C66 +C3*Td(3,1)*k_620 +C3*Td(3,2)*k_20 -C4*Td(3,1)*k_10 -C4*Td(3,2)*k_610

		Z_mat(6,1)  = C62*k_50 -C63*k_10 +C65*k_40 -C66*k_620 -C1*Td(3,2)*k_40dx +C1*Td(3,3)*k_620dx +C2*Td(3,2)*k_50dx &
				 -C2*Td(3,3)*k_10dx+C3*Td(3,2)*k_40dy-C3*Td(3,3)*k_620dy-C4*Td(3,2)*k_50dy+C4*Td(3,3)*k_10dy
		Z_mat(6,7)  = -C61*k_50 -C63*k_610 -C64*k_40 -C66*k_20 +C1*Td(3,1)*k_40dx +C1*Td(3,3)*k_20dx -C2*Td(3,1)*k_50dx &
				  -C2*Td(3,3)*k_610dx-C3*Td(3,1)*k_40dy-C3*Td(3,3)*k_20dy+C4*Td(3,1)*k_50dy+C4*Td(3,3)*k_610dy
		Z_mat(6,13) = C61*k_10 +C62*k_610 +C64*k_620 +C65*k_20 -C1*Td(3,1)*k_620dx -C1*Td(3,2)*k_20dx +C2*Td(3,1)*k_10dx &
				  +C2*Td(3,2)*k_610dx+C3*Td(3,1)*k_620dy+C3*Td(3,2)*k_20dy-C4*Td(3,1)*k_10dy-C4*Td(3,2)*k_610dy

		Z_mat(6,5)   =-C1*Td(3,1) -C4*Td(3,1)
		Z_mat(6,11)  =-C1*Td(3,2) -C4*Td(3,2)
		Z_mat(6,17)  =-C1*Td(3,3) -C4*Td(3,3)
		Z_mat(6,6)   = C3*Td(3,1)
		Z_mat(6,12)  = C3*Td(3,2)
		Z_mat(6,18)  = C3*Td(3,3)
		Z_mat(6,4)   = C2*Td(3,1)
		Z_mat(6,10)  = C2*Td(3,2)
		Z_mat(6,16)  = C2*Td(3,3)
		Z_mat(7,20)  = 1._dp
		Z_mat(8,21)  = 1._dp
		Z_mat(9,23)  = 1._dp
		Z_mat(10,24) = 1._dp
		Z_mat(11,19) = 1._dp
		Z_mat(12,22) = 1._dp
	end subroutine getZ

end module getZ_mod


