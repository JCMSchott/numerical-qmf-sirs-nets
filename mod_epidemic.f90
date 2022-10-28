module mod_epidemic
	use geraRede
	use mod_tools
	implicit none
	
	contains
	
	
		!#################################################
		!	Modelo SIRS
		!#################################################		
					
    subroutine rk4_SIRS_grafo(this, t, dt, n1, alp, lam, mu, Si, Ii) !, f_Si, f_Ii
		    use types
		    !use geraRede
		
		    class(grafo) :: this
         	real(dp) :: t, dt
        	integer :: m1, n1
        	real(dp), intent(in) :: alp, lam, mu
        	real(dp) :: Si(n1), Ii(n1)
        	real(dp) :: k1_Si(n1), k2_Si(n1), k3_Si(n1), k4_Si(n1)
        	real(dp) :: k1_Ii(n1), k2_Ii(n1), k3_Ii(n1), k4_Ii(n1)
		 
  		!####################################################################################
			
		k1_Si = dt * f_Si(this, t, Si, Ii, n1, alp, lam)
		
		k1_Ii = dt * f_Ii(this, t, Si, Ii, n1, mu, lam)

		
		k2_Si = dt * f_Si(this, t + 0.5_dp * dt, Si + 0.5_dp * k1_Si, Ii + 0.5_dp * k1_Ii, n1, alp, lam)		
		
		k2_Ii = dt * f_Ii(this, t + 0.5_dp * dt, Si + 0.5_dp * k1_Si, Ii + 0.5_dp * k1_Ii, n1, mu, lam)



		k3_Si = dt * f_Si(this, t + 0.5_dp * dt, Si + 0.5_dp * k2_Si, Ii + 0.5_dp * k2_Ii, n1, alp, lam)

		k3_Ii = dt * f_Ii(this, t + 0.5_dp * dt, Si + 0.5_dp * k2_Si, Ii + 0.5_dp * k2_Ii, n1, mu, lam)



		k4_Si = dt * f_Si(this, t + dt, Si + k3_Si, Ii + k3_Ii, n1, alp, lam)
		
		k4_Ii = dt * f_Ii(this, t + dt, Si + k3_Si, Ii + k3_Ii, n1, mu, lam)

			
		Si = Si + (1_dp/6_dp) * (k1_Si + 2_dp * k2_Si + 2_dp * k3_Si + k4_Si)  		
  		
		
		Ii = Ii + (1_dp/6_dp) * (k1_Ii + 2_dp * k2_Ii + 2_dp * k3_Ii + k4_Ii)  		
  		!####################################################################################
  		! F_Ii
  	  	        		    
        end subroutine		


	function f_Si(this, t, Si, Ii, n_sitios, alp, lam)

		use types
		!use geraRede

		class(grafo) :: this
		real(dp), intent(in) :: t, alp, lam
		integer, intent(in) :: n_sitios
		real(dp) :: Si(n_sitios), Ii(n_sitios)
		real(dp) :: f_Si(n_sitios)
		real(dp) :: sumIi
		!###############################
		integer :: l1, l12, l2
		!###############################
		
		! Sitio i
		do l1 = 1, this%nodes
		    if(lista_de_clusters(l1) /= i_comp_gigante) cycle				
			f_Si(l1) = alp * (1.0_dp - Si(l1) - Ii(l1) )
			
			sumIi = 0.0_dp

			do l12 = this%aux(l1), this%aux(l1) + this%deg(l1) - 1
				l2 = this%listAdj(l12)
				sumIi = sumIi + Ii(l2)
			enddo
						
			f_Si(l1) = f_Si(l1) - lam * Si(l1) * sumIi
		enddo
	end function			

	!#################################################
	!	F_Ii
	!#################################################

	function f_Ii(this, t, Si, Ii, n_sitios, mu, lam)
		use types
		!use geraRede
		
		
		class(grafo) :: this
		real(dp), intent(in) :: t, mu, lam
		integer, intent(in) :: n_sitios
		real(dp) :: Si(n_sitios), Ii(n_sitios)
		real(dp) :: f_Ii(n_sitios)
		real(dp) :: sumIi
		!###############################
		integer :: l1, l12, l2, l3
		!###############################
		
		! Sitio i.
		do l1 = 1, this%nodes
			f_Ii(l1) = -mu * Ii(l1)
			
			if(lista_de_clusters(l1) /= i_comp_gigante) cycle
			
			sumIi = 0.0_dp
			do l12 = this%aux(l1), this%aux(l1) + this%deg(l1) - 1
				l2 = this%listAdj(l12)
				sumIi = sumIi +  Ii(l2)
			enddo
			
			f_Ii(l1) = f_Ii(l1) + lam * Si(l1) * sumIi

		enddo
	end function
     	
end module
