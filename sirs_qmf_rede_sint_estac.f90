module mod_epidemic
	use geraRede
	use mod_tools
	use types
	implicit none
	
	real(dp), allocatable :: Ii(:), Si(:)
	real(dp), allocatable :: k1_Ii(:), k1_Si(:)
	real(dp), allocatable :: k2_Ii(:), k2_Si(:)
	real(dp), allocatable :: k3_Ii(:), k3_Si(:)
	real(dp), allocatable :: k4_Ii(:), k4_Si(:)
	
	contains
	
	
		!#################################################
		!	Modelo SIRS
		!#################################################		

    subroutine aloca(this)
       class(grafo) :: this
       
       
       if(allocated(k1_Ii))then
          deallocate(Ii)
          deallocate(Si)
          
          deallocate(k1_Ii)
          deallocate(k2_Ii)
          deallocate(k3_Ii)
          deallocate(k4_Ii)
          deallocate(k1_Si)
          deallocate(k2_Si)
          deallocate(k3_Si)
          deallocate(k4_Si)
       endif          
          
       allocate(Ii(this%nodes))
       allocate(k1_Ii(this%nodes))
       allocate(k2_Ii(this%nodes))   
       allocate(k3_Ii(this%nodes))
       allocate(k4_Ii(this%nodes))       

       allocate(Si(this%nodes))   
       allocate(k1_Si(this%nodes))
       allocate(k2_Si(this%nodes))   
       allocate(k3_Si(this%nodes))
       allocate(k4_Si(this%nodes))
    
    end subroutine
					
    subroutine rk4_SIRS_grafo(this, t, dt, n1, alp, lam, mu) !, f_Si, f_Ii
		    use types
		    !use geraRede
		
		    class(grafo) :: this
         	real(dp) :: t, dt
        	integer :: n1
        	real(dp) :: alp, lam, mu		 
  		!####################################################################################
		
		
		!####################################################################################	
		!            f_Si(this, t, Si_1, Ii_1, n_sitios, alp, lam)
		k1_Si = dt * f_Si(this, t, Si, Ii, n1, alp, lam)
		
		!            f_Ii(this, t, Si, Ii, n_sitios, mu, lam)
		k1_Ii = dt * f_Ii(this, t, Si, Ii, n1, mu, lam)

		
		!####################################################################################
		
		!            f_Si(this, t              , Si_1               , Ii_1               , n_sitios, alp, lam)
		k2_Si = dt * f_Si(this, t + 0.5_dp * dt, Si + 0.5_dp * k1_Si, Ii + 0.5_dp * k1_Ii, n1, alp, lam)		
		
		!            f_Ii(this, t,               Si                 , Ii                 , n_sitios, mu, lam)
		k2_Ii = dt * f_Ii(this, t + 0.5_dp * dt, Si + 0.5_dp * k1_Si, Ii + 0.5_dp * k1_Ii, n1, mu, lam)


        !####################################################################################
		!            f_Si(this, t              , Si_1               , Ii_1               , n_sitios, alp, lam)        
		k3_Si = dt * f_Si(this, t + 0.5_dp * dt, Si + 0.5_dp * k2_Si, Ii + 0.5_dp * k2_Ii, n1, alp, lam)

        !            f_Ii(this, t,               Si                 , Ii                 , n_sitios, mu, lam)
		k3_Ii = dt * f_Ii(this, t + 0.5_dp * dt, Si + 0.5_dp * k2_Si, Ii + 0.5_dp * k2_Ii, n1, mu, lam)


       !####################################################################################
		!            f_Si(this, t     , Si_1      , Ii_1      , n_sitios, alp, lam)       
		k4_Si = dt * f_Si(this, t + dt, Si + k3_Si, Ii + k3_Ii, n1, alp, lam)
		
		!            f_Ii(this, t,      Si        , Ii        , n_sitios, mu, lam)
		k4_Ii = dt * f_Ii(this, t + dt, Si + k3_Si, Ii + k3_Ii, n1, mu, lam)



		!####################################################################################	
		Si = Si + (1.0_dp/6.0_dp) * (k1_Si + 2.0_dp * k2_Si + 2.0_dp * k3_Si + k4_Si)  		
  		
		
		Ii = Ii + (1.0_dp/6.0_dp) * (k1_Ii + 2.0_dp * k2_Ii + 2.0_dp * k3_Ii + k4_Ii)  		
  		!####################################################################################
  		! F_Ii
  	  	        		    
    end subroutine		


	function f_Si(this, t, Si_1, Ii_1, n_sitios, alp, lam)

		!use types
		!use geraRede

		class(grafo) :: this
		real(dp), intent(in) :: t, alp, lam
		integer, intent(in) :: n_sitios
		real(dp), intent(in) :: Si_1(n_sitios), Ii_1(n_sitios)
		real(dp) :: f_Si(n_sitios)
		real(dp) :: sumIi
		!###############################
		integer :: l1, l12, l2
		!###############################
		
		! Sitio i
		do l1 = 1, this%nodes
		    if(lista_de_clusters(l1) /= i_comp_gigante) cycle				
			f_Si(l1) = alp * (1.0_dp - Si_1(l1) - Ii_1(l1) )
			
			sumIi = 0.0_dp

			do l12 = this%aux(l1), this%aux(l1) + this%deg(l1) - 1
				l2 = this%listAdj(l12)
				sumIi = sumIi + Ii_1(l2)
			enddo
						
			f_Si(l1) = f_Si(l1) - lam * Si_1(l1) * sumIi
			
		enddo
	end function			

	!#################################################
	!	F_Ii
	!#################################################

	function f_Ii(this, t, Si_1, Ii_1, n_sitios, mu1, lam1)
		!use types
		!use geraRede
		
		
		class(grafo) :: this
		real(dp), intent(in) :: t, mu1, lam1
		integer, intent(in) :: n_sitios
		real(dp), intent(in) :: Si_1(n_sitios), Ii_1(n_sitios)
		real(dp) :: f_Ii(n_sitios)
		real(dp) :: sumIi
		!###############################
		integer :: l1, l12, l2, l3
		!###############################
		
		! Sitio i.
		do l1 = 1, this%nodes			
			if(lista_de_clusters(l1) /= i_comp_gigante) cycle
			f_Ii(l1) = -mu1 * Ii_1(l1)
			sumIi = 0.0_dp
			do l12 = this%aux(l1), this%aux(l1) + this%deg(l1) - 1
				l2 = this%listAdj(l12)
				sumIi = sumIi +  Ii_1(l2)
			enddo
			
			f_Ii(l1) = f_Ii(l1) + lam1 * Si_1(l1) * sumIi
			
			!write(*,*) "f_Ii = ", f_Ii(l1)

		enddo
	end function
     	
end module


program sirsPadrao

!####################################################
	use mod_epidemic
	use geraRede
	use mod_rndgen
	use types
	use mod_tools
!####################################################
! Modulos utilizados
	
	implicit none

!####################################################
	real(dp) :: t0, tf, t, dt
!####################################################
! Tempo discretizado e incremento

!########################################################	
	real(dp) :: Si_medio, Ii_medio, Ri_medio
	real(dp) :: Si_m, Ii_m, Ri_m
	real(dp) :: Si_m0, Ii_m0, Ri_m0
	integer ::  n_sitios(6), n_pts, semente
!########################################################
! Variaveis dinamicas

!########################################################
	real(dp) :: alp, lambda
	real(dp), parameter :: mu = 1.0_dp
!########################################################
! Taxas de transicao

!########################################################
	type(grafo_PL_UCM) :: rede
!########################################################
! Declara variaveis tipo grafo


!########################################################
	type(rndgen) :: gen
!########################################################
! Gerador de numeros pseudo-aleatorios

!########################################################
	real(dp) :: p0
	real(dp) :: p_sort
!########################################################
! Probabilidade inicial e numero do sorteio

!########################################################
	integer :: i1, i2, i3, i4, j1
!########################################################
! Variaveis auxiliares	

!########################################################
	character(len=30) :: rowfmt
!########################################################
! Tamanho do espaco requerido para printar cada entrada
! no arquivo

!########################################################
	integer, parameter :: sit_escol = 100
!########################################################
! Como sou estravagante, escolho 100 sitios para plotar
! seus estados

!########################################################
	character(len=500) :: caminho
!########################################################


	logical :: eh_absv1, eh_absv2
	
	real(dp) :: t_relax
	real(dp) :: err_test, Ii_mold, Ii_st
	!real(dp) :: exp_gama(3)
	real(dp) :: exp_gama
	real(dp) :: toll
	real(dp) :: lambda0, lambdaf, dlambda	
	integer :: nlambda
	integer :: j2, j3, j4
	character(len=8) :: N_char, gama_char, alp_char, tam_char
	real(dp) :: kmaxi
	integer :: repetiu
	!###################################################################
	character(len=1) :: flag_char
	integer :: flag1
    character(len=300) :: buffer
    integer :: n_args
    character(len=100) :: arquivo_char1, arquivo_rede
    !character(len=30) :: tam_char
    character(len=10) :: lambda_char
    integer :: label
    character(len=1000) :: local_arquivo, I_i_arquivo
    character(len=500) :: local
    integer :: sum_deg
    character(len=1000) :: lb_vs_I_i_arquivo, t_vs_I_im_arquivo
    real(dp) :: trec, trec1
    integer :: num_rec
    integer, allocatable :: tam_rede(:)
    integer :: m1
    
    integer(kind=8) :: t1, t2, taxa
    real(dp) :: te
    integer :: k_i
    real(dp) :: k_f    
    integer :: seed(10)
    type(rndgen) :: ger_inic
    character(len=5) :: indice
    real(dp), allocatable :: P_grau(:)
    character(len=500) ::arq_1
    real(dp), allocatable :: Ap_list(:)
    character(len=300) :: cwd, resultados
	!###################################################################
		
    call getcwd(cwd)
    
    cwd = trim(adjustl(cwd))
    
    resultados = 'Rst_QMF'
    
    resultados = trim(adjustl(resultados))
    	
	call system('mkdir -p '//trim(adjustl(resultados)))
	    
	local = trim(adjustl(trim(adjustl(cwd))//"/"//trim(adjustl(resultados))//"/"))
    !###################################################################
    
	t0 = 0.0_dp;	!tf = 1000.0_dp;
	
	dt = (10.0_dp)**(-2.0_dp)

	n_pts = 10**5

	tf = t0 + 1.0_dp * n_pts * dt
	
	trec = tf/1000.0_dp
	
	trec1 = trec    
    !###################################################################
!#######################   
   alp = 1.0_dp
   
   toll = (10.0_dp) ** (-5.0_dp)
!#######################################################################

   nlambda = 1000
   dlambda = 0.0125_dp/5.0_dp   
   lambda0 = 0.001_dp !5.0_dp * dlambda +0.5_dp
   
   lambda = lambda0
   lambdaf = 1.0_dp
!#######################   
   alp = 0.50_dp  
!#######################################################################
!  Indice da semente
   i4 = 1
!####################################################################### 
   
!################################################################################   
! DEFINO SEMENTE, GRAU MINIMO, TAMANHO DA REDE E EXP_GAMA
!#######################
   semente  = 947361823
!#######################
   k_i = 3   
!#######################   
   !tam_rede = 1000
   allocate(tam_rede(8))
   
   tam_rede(1) = 10**3; tam_rede(2) = 3*10**3
   tam_rede(3) = 10**4; tam_rede(4) = 3*10**4   
   tam_rede(5) = 10**5; tam_rede(6) = 3*10**5
   tam_rede(7) = 10**6; tam_rede(8) = 3*10**6   
!#######################   
   m1 = 3
   
   p0 = 10.0_dp/(1.0_dp * tam_rede(m1))
      
   !exp_gama = 2.3_dp
   !exp_gama = 2.7_dp
   exp_gama = 2.8_dp
   !exp_gama = 3.5_dp

!#############################################################################################
!#######################
   write(gama_char, '(f5.2)') exp_gama         
   gama_char = trim(adjustl(gama_char))
!#######################

   write(alp_char, '(f5.2)') alp
   alp_char = trim(adjustl(alp_char)) 
!##########################################################################################

if(.True.)then 
   if(tam_rede(m1) == 100)then
      tam_char = '0_1k'
   elseif(tam_rede(m1) == 1000)then
      tam_char = '1k'
   elseif(tam_rede(m1) == 3000)then
      tam_char = '3k'
   elseif(tam_rede(m1) == 5000)then
      tam_char = '5k'      
   elseif(tam_rede(m1) == 10000)then
      tam_char = '10k'
   elseif(tam_rede(m1) == 30000)then
      tam_char = '30k'
   elseif(tam_rede(m1) == 100000)then
      tam_char = '100k'
   elseif(tam_rede(m1) == 300000)then
      tam_char = '300k'
   elseif(tam_rede(m1) == 1000000)then
      tam_char = '1M'
   elseif(tam_rede(m1) == 3000000)then
      tam_char = '3M'
   elseif(tam_rede(m1) == 10000000)then
      tam_char = '10M'
   elseif(tam_rede(m1) == 30000000)then
      tam_char = '30M'
   else
      stop 'Escolha um tamanho de rede dentro do catalogo'
   endif
endif   


!##########################################################################################
!				Inicia grafo
!##########################################################################################
if(.True.)then

   if(exp_gama > 3.0_dp)then
      k_f = 2.0_dp * (1.0_dp * tam_rede(m1))**(0.5_dp)
      !k_f = (1.0_dp * tam_rede(m1))**(1_dp/(exp_gama - 1.0_dp))
      !call acha_cutoff_rigido(k_i, exp_gama, tam_rede(m1))
      !k_f = 1.0_dp * (1.0_dp * tam_rede(m1))**(1.0_dp/(exp_gama -1.0_dp))
   else
      k_f = 2.0_dp * (1.0_dp * tam_rede(m1))**(0.5_dp)
   endif 
 
   call ger_inic%init(semente)
   i2 = 1
   do i1 = 1, 1000
      if(mod(i1,100) > 0) cycle
      seed(i2)  = ger_inic%int(100000000,999999999)
      write(*,*) i1, seed(i2)
      i2 = i2+1      
   enddo 

!##########################################################################################
!				Inicia grafo
!##########################################################################################
 !  arquivo_rede_real='s00088.s838.net.edg'
 !  call rede%RedeReal(arquivo_rede_real, 111)
 !  close(111)
!#######################################################################   
  call rede%iniciaGrafo(tam_rede(m1))
!#######################################################################
  call rede%inicia(k_i, k_f, exp_gama, seed(i4))
!#######################################################################
  write(indice,'(I0)') i4
!#######################################################################
  indice = trim(adjustl(indice))
!####################################################################### 
  indice = trim(adjustl('_'//trim(adjustl(indice)))) 
!#######################################################################
  call rede%liga(seed(i4), .False.)
!#######################################################################  
  write(*,*) "Tamanho da rede eh ", rede%nodes  
!#######################################################################
   call sub_classifica_clusters(rede,.False., 000, 'sem_arquivo.dat')
   write(*,*) "Calculou a comp gigante, de indice: ", i_comp_gigante
   !stop "Concluiu a ligacao da rede"
!#######################################################################
   arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_clustering_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//trim(adjustl(indice))//'.dat'
   call clustering(rede,.True., 800, trim(adjustl(arq_1)))
   close(800)
!#######################################################################
   arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_knn_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//trim(adjustl(indice))//'.dat'
   call calcula_k_nn(rede, .True., 800, trim(adjustl(arq_1)))
   close(800)
!#######################################################################
   if(allocated(P_grau)) deallocate(P_grau)
!#######################################################################
   allocate(P_grau(rede%degMin:rede%degMax))
!#######################################################################
   P_grau = 0.0_dp
!#######################################################################   
   do i1 = 1, rede%nodes
      P_grau(rede%deg(i1)) = P_grau(rede%deg(i1)) + 1.0_dp
   enddo
!#######################################################################   
   P_grau = P_grau/(1.0_dp * rede%nodes)
!#######################################################################   
   arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_P_grau_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//trim(adjustl(indice))//'.dat'
   open(800, file=trim(adjustl(arq_1)), status='unknown')
!#######################################################################   
   do i1 = rede%degMin, rede%degMax
      write(800,*) i1, P_grau(i1)
   enddo
!#######################################################################   
   close(800)
   deallocate(P_grau)
!#######################################################################    
   write(*,*) "Gerou a rede"
   write(*,*) ""
!#######################################################################
   arquivo_char1 = trim(adjustl('lbd_vs_t_convQMF_tam_'//trim(adjustl(tam_char))//'_gm_'//trim(adjustl(gama_char))//'_alp_'//trim(adjustl(alp_char))))
   local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//trim(adjustl(indice))//'.dat'
   open(unit = 888, file = trim(adjustl(local_arquivo)), status = 'unknown')
!#######################################################################   		
   write(*,*) ""
   write(*,*)"O tamanho da rede eh: ", comp_gigante
   write(*,*) ""
   write(*,*)"O grau minimo da rede eh: ", rede%degMin
   write(*,*) ""
   write(*,*)"O grau maximo da rede eh: ", rede%degMax
   write(*,*) ""
   write(*,*)"O grau medio da rede eh: ", rede%degMean
   write(*,*) ""
            
   sum_deg = sum(rede%deg)
endif   
!##########################################################################################         
! DEFINE AS STRINGS
!#######################           
   tam_char = trim(adjustl(tam_char))

!#######################################################################
!   Aloca as listas de estados
				
!#######################################################################
! Setta o lambda inicial
			lambda = lambda0
!#######################################################################
            arquivo_char1 = trim(adjustl('lbd_vs_rhoQMF_tam_'//trim(adjustl(tam_char))//'_gm_'//trim(adjustl(gama_char))//'_alp_'//trim(adjustl(alp_char))))
            
            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//trim(adjustl(indice))//'.dat'
            t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
            open(unit = 13, file = trim(adjustl(t_vs_I_im_arquivo)), access='append', status = 'unknown')	




            call aloca(rede)
!#######################################################################
			do j4 = 1, nlambda
			   !write(lambda_char, '(f7.4)') lambda
			   if(lambda >= lambdaf) exit
!#######################################################################
               arquivo_char1 = trim(adjustl('t_vs_rhoQMF_tam_'//trim(adjustl(tam_char))//'_gm_'//trim(adjustl(gama_char))//'_alp_'//trim(adjustl(alp_char))))
            
               local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//trim(adjustl(indice))//'.dat'
               t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
               open(unit = 31, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')				
!#######################################################################
               arquivo_char1 = trim(adjustl('t_vs_drho_rho_qmf_tam_'//trim(adjustl(tam_char))//'_gm_'//trim(adjustl(gama_char))//'_alp_'//trim(adjustl(alp_char))))
            
               local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//trim(adjustl(indice))//'.dat'
               t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
               open(unit = 32, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')
!#######################################################################									
				!################################################
				! Porcentagem inicial de infectados
				!################################################

				do i1 = 1, rede%nodes
					Ii(i1) = p0
					Si(i1) = 1.0_dp - Ii(i1)				 				
				enddo

				!#######################################################
				! Setta as condicoes iniciais

				!#######################################################
!#######################################################################				
			
!#######################################################################
				Ii_m = p0
				Si_m = 1.0_dp - p0

				Ri_m = 1.0_dp - Ii_m - Si_m
                
                !write(*,*) "Ii_m = ", Ii_m, "Si_m = ", Si_m, "Ri_m = ", Ri_m
                
				repetiu = 0
				t = t0

				Ii_m0 = Ii_m
				Si_m0 = Si_m
					
				Ri_m0 = Ri_m
				!#######################################################
loop_t:			do i1 = 1, n_pts					
				    !###################################################

                    write(31,*) t, Ii_m
					!###################################################
					! Calcula as medias das variaveis dinamicas
					!##############################################################$$$$$$$$$$$$$$$$
					
					!    rk4_SIRS_grafo(this, t, dt, n1        , alp, lam   , mu)
					call rk4_SIRS_grafo(rede, t, dt, rede%nodes, alp, lambda, mu)	
					
    				Si_m = 0.0_dp
	    			Ii_m = 0.0_dp
	    			do j2 = 1, rede%nodes
		    		   if(lista_de_clusters(j2) /= i_comp_gigante) cycle
			    	   Ii_m = Ii_m + Ii(j2)
			    	   Si_m = Si_m + Si(j2)
			    	   !if(Si(j2) > 1.0_dp)then
			    	      !stop "Si > 1"
			    	   !elseif(Ii(j2) < 0.0_dp)then
			    	      !stop "Ii negativo"
			    	   !elseif( (1.0_dp - Ii(j2)-Si(j2)) < 0.0_dp)then
			    	      !write(*,*) "Tempo = ", t
			    	      !write(*,*) "i_lambda = ", j4
			    	      !write(*,*) "i_tempo = ", i1
			    	      !write(*,*) "sitio = ", j2			    	   
			    	      !write(*,*) "Ri = ", 1.0_dp - Ii(j2) - Si(j2)
			    	      !stop "Ri negativo"
			    	   !endif
				    enddo
				    Ii_m = 1.0_dp * Ii_m/comp_gigante
				    Si_m = 1.0_dp * Si_m/comp_gigante				
				    Ri_m = 1.0_dp - Si_m - Ii_m
				
                    write(32,*) t, abs(Ii_m-Ii_m0)/Ii_m0
											
               if(abs(Ii_m-Ii_m0)/Ii_m0**2.0_dp < toll)then
                  repetiu = repetiu + 1
                  if(repetiu > 1000)then
                     write(*,*) "Ii_m0 = ", Ii_m0, " Ii_m = ", Ii_m
                     write(*,*) "Conv"
                     exit loop_t
                  endif
               else                  
                  if(repetiu < 0)then
                     !write(*,*) "Ainda nao convergiu"                  
                     !write(*,*) "Repetiu = ", repetiu
                     repetiu = 0
                     !write(*,*) "Ii_m0 = ", Ii_m0, " Ii_m = ", Ii_m                    
                  endif
               endif
			   
			   Ii_m0 = Ii_m
			   
			   if(Ii_m < 1.0_dp/(1.0_dp *rede%nodes))then
                  write(*,*) "Abs sub"
                  exit loop_t
               endif		
			   
			   !##############################################################
               t = t + dt
               !##############################################################
			   ! Incrementa o tempo.			
               					
			   enddo loop_t
				
				close(31)
				close(32)
				close(33)
				!########################################################
	         !Rotinaprincipal
	            if (repetiu > 1000)then
	               write(13,*) lambda, Ii_m
	            endif
	                 					
				lambda = lambda + dlambda
				
				!########################################################
					
				!########################################################
				! Fecha os arquivos
		
			enddo	

			close(13)
			close(14)
			close(15)
!enddo tmn1
		!close(14)


	!######################################################################
	
	contains

   subroutine acha_cutoff_rigido(kming, gam_p, N_p)
      integer :: kming
      real(dp) :: gam_p
      integer :: N_p
      real(dp), parameter :: tol = 5d-5
      real(dp) :: gminus, gmais
      real(dp) :: kminus, kmais
      real(dp) :: k_p
      real(dp) :: gp
      real(dp) :: Ap1
      integer :: kl1, kl2
      integer, parameter :: N_it = 10**3
       
      !#################################################################
      !   Inicio
      !#################################################################           
      kminus = 1d0 * kming
      kmais = 1.5d0 * kming * (1d0 * N_p)**(1d0/gam_p)
      
      if(allocated(Ap_list)) deallocate(Ap_list)
      allocate(Ap_list(int(kming):(int(kmais))))
      
      gminus = g_func(kminus, gam_p, N_p)
      
      if(gminus >= 0d0) stop "Precisamos de um valor de kminus para que gminus < 0"
       
      Ap1 = 0d0
      do kl1 = kming, int(kmais)
         Ap1 = Ap1 + (1d0 * kl1)**(-gam_p)
         Ap_list(kl1) = Ap1
      enddo

      gmais = g_func(kmais, gam_p, N_p)
   
      if(gmais <= 0d0) stop "Precisamos de um valor para kmais, tal que gmais > 0"
      !#################################################################
      !   Execucao
      !#################################################################
      kl1 = 1
      do while(kl1 <= N_it)
         k_p = kminus + (kmais - kminus)/2d0
         gp = g_func(k_p, gam_p, N_p)
         if((gp == 0d0).or.((kmais-kminus)/2d0 <= tol))then
            k_f = k_p
            write(*,*) "Achou o grau max apos N = ", kl1, " iteracoes."
            exit
         endif
         kl1 = kl1 + 1
         if(gminus * gp > 0d0)then
            kminus = k_p
            gminus = gp
         else
            kmais = k_p
         endif
      enddo      
      !#################################################################
      !   Final
      !#################################################################
      deallocate(Ap_list)                   
   end subroutine
   
   function g_func(k_s, gam1, Nstr)
      real(dp) :: k_s
      real(dp) :: gam1
      integer :: Nstr
      real(dp) :: g_func
      real(dp) :: Ap
      !#################################################################
      !   Quando g_func for usada pela primeira vez,
      !   gbuffer = 0.0d0 e kmin = this%degMin,
      !   k_s recebe o valor que se quer aproximar de kc
      !   
      !#################################################################      
      !g_func = gbuffer
                  
      Ap = Ap_list(int(k_s))
      
      Ap = 1d0/Ap
      
      g_func = k_s - Ap**(1d0/gam1) * (1d0 * Nstr)**(1d0/gam1)
      
   end function 

end program
