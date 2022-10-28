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
	real(dp), allocatable :: Si(:), Ii(:), Ri(:)
	real(dp) :: Si_medio, Ii_medio, Ri_medio
	real(dp) :: Si_m, Ii_m, Ri_m
	real(dp) :: Si_m0, Ii_m0, Ri_m0
	integer ::  n_sitios(6), n_pts, semente
!########################################################
! Variaveis dinamicas

!########################################################
	real(dp) :: alp, lam, mu
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
	real(dp), parameter :: p0 = 0.01_dp
	real(dp) :: p_sort
!########################################################
! Probabilidade inicial e numero do sorteio

!########################################################
	integer :: i1, i2, i3, j1
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
	real(dp) :: exp_gama(3)
	real(dp), parameter :: toll = 5d-10
	real(dp) :: lam0, lamf, dlam	
	integer :: n_lam
	integer :: j2, j3, j4
	character(len=8) :: N_char, gama_char
	real(dp) :: kmaxi
	integer :: repetiu
	!###################################################################
	character(len=1) :: flag_char
	integer :: flag1
    character(len=300) :: buffer
    integer :: n_args
    character(len=100) :: arquivo_char1, arquivo_rede
    character(len=30) :: Ni_char
    character(len=10) :: lambda_char
    integer :: label
    character(len=1000) :: local_arquivo, I_i_arquivo
    character(len=500) :: local
    integer :: sum_deg
    character(len=1000) :: lb_vs_I_i_arquivo, t_vs_I_im_arquivo
    real(dp) :: trec, trec1
    integer :: num_rec
	!###################################################################
	    
	caminho = trim(adjustl("/home/jota/SIRS_QMF_Rede_Real/QMF_Rst_Rede_Real/"))
	local = trim(adjustl("/home/jota/SIRS_QMF_Rede_Real/QMF_Rst_Rede_Real/"))
!	caminho = trim(adjustl("/home/jota/SIRS_QMF/S_num/SIRS_padrao/Res_est/"))
    !###################################################################
    
	t0 = 0.0_dp;	!tf = 1000.0_dp;
	
	dt = 1d-3

	n_pts = 50000

	tf = t0 + 1.0_dp * n_pts * dt
	
	trec = 1.0_dp * tf/1000
	
	trec1 = trec
    
    !###################################################################
    	
	lam0 = 0.005_dp ! 0.4_dp !
	lamf = 1_dp
	n_lam = 1 !000
		
	mu = 1.0_dp; alp = 100.5_dp 
	dlam = 1.0_dp * (lamf - lam0)/n_lam	
    !###################################################################

!########################################################
! Define parametros


!########################################################
					
		
			!######################################
			semente = 995887671
			!######################################
			! Inicializa a semente

			!######################################
			eh_absv1 = .True. ; eh_absv2 = .True.
			!######################################
			! Chaves para garantir a unicidade
			! do ponto critico
			
!################################################################################   
   semente = 995887671   
!################################################################################   

!################################################################################         
! Coleta a chave 0 ou 1 e o arquivo da rede real
            n_args = iargc()

            if(n_args == 2)then
               call getarg(1, flag_char)
               read(flag_char,*) flag1
               if((flag1 /= 1).and.(flag1 /=0)) stop "Argumento invalido"
               call getarg(2,buffer)
	           arquivo_rede = trim(adjustl(buffer))
            else
               stop "Numero de argumentos invalido"
            endif
!################################################################################   
		    
            label = 100
            write(*,*) arquivo_rede
            open(unit=label, file=arquivo_rede,status='old')

!################################################################################
!   Gera a rede
            call rede%RedeReal(arquivo_rede, label)
            write(*,*) "Gerou a rede"
            write(*,*) ""
            close(label)
!################################################################################
!   Classifica os fragmentos da rede         
            call sub_classifica_clusters(rede, .False., 000, 'nenhum.dat')
            write(Ni_char,*) comp_gigante
            Ni_char = trim(adjustl(Ni_char))
!################################################################################
!   Soma os graus dos sitios da componente gigante
            sum_deg = 0
            do i1 = 1, rede%nodes
               if(lista_de_clusters(i1) /= i_comp_gigante)cycle
               sum_deg = sum_deg + rede%deg(i1)
            enddo
            sum_deg = sum(rede%deg)

            write(*,*) "Calculou a componente gigante"
            write(*,*) ""
		    
			rede%degMin = minval(rede%deg)		
			rede%degMax = maxval(rede%deg)
!#######################################################################
!   Aloca as listas de estados
			if(allocated(Si)) deallocate(Si)
				allocate(Si(rede%nodes))		

			if(allocated(Ii)) deallocate(Ii)
				allocate(Ii(rede%nodes))
				
!#######################################################################
! Setta o lambda inicial
			lam = lam0
!#######################################################################
            arquivo_char1 = trim(adjustl('lambda_vs_I_im_tam_'//trim(adjustl(Ni_char))//'_rede_'//trim(adjustl(arquivo_rede))))
            
            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
            t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
            open(unit = 13, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')	
!#######################################################################
            arquivo_char1 = trim(adjustl('t_vs_I_i_tam_'//trim(adjustl(Ni_char))//'_rede_'//trim(adjustl(arquivo_rede))))
            
            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
            t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
            open(unit = 14, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')	
!#######################################################################
            arquivo_char1 = trim(adjustl('t_vs_R_i_tam_'//trim(adjustl(Ni_char))//'_rede_'//trim(adjustl(arquivo_rede))))
            
            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
            t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
            open(unit = 15, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')
!#######################################################################
            arquivo_char1 = trim(adjustl('t_vs_S_i_tam_'//trim(adjustl(Ni_char))//'_rede_'//trim(adjustl(arquivo_rede))))
            
            local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
            t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
            open(unit = 16, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')
            
!#######################################################################
			do j4 = 1, n_lam
			   write(lambda_char, '(f7.4)') lam
!#######################################################################
               arquivo_char1 = trim(adjustl('t_vs_I_im_tam_'//trim(adjustl(Ni_char))//'_rede_'//trim(adjustl(arquivo_rede))//'_lbd_'//trim(adjustl(lambda_char))))
            
               local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
               t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
               open(unit = 31, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')				
!#######################################################################
               arquivo_char1 = trim(adjustl('t_vs_dIm_Im_tam_'//trim(adjustl(Ni_char))//'_rede_'//trim(adjustl(arquivo_rede))//'_lbd_'//trim(adjustl(lambda_char))))
            
               local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
               t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
               open(unit = 32, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')
!#######################################################################
               arquivo_char1 = trim(adjustl('t_vs_R_im_tam_'//trim(adjustl(Ni_char))//'_rede_'//trim(adjustl(arquivo_rede))//'_lbd_'//trim(adjustl(lambda_char))))
            
               local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
               t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
               open(unit = 33, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')				               
!#######################################################################
               arquivo_char1 = trim(adjustl('t_vs_S_im_tam_'//trim(adjustl(Ni_char))//'_rede_'//trim(adjustl(arquivo_rede))//'_lbd_'//trim(adjustl(lambda_char))))
            
               local_arquivo = trim(adjustl(local))//trim(adjustl(arquivo_char1))//'.dat'
               t_vs_I_im_arquivo = trim(adjustl(local_arquivo))
               open(unit = 34, file = trim(adjustl(t_vs_I_im_arquivo)), status = 'unknown')
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
				Ii_m = p0 !0.0_dp
				Si_m = 1.0_dp - p0

				Ri_m = 0.0_dp

				repetiu = 0
				t = t0
				num_rec = 1
				!#######################################################
loop_t:			do i1 = 1, n_pts
					
				    !###################################################
					
					Ii_m0 = Ii_m
					Si_m0 = Si_m
					
					Ri_m0 = Ri_m

                    write(31,*) t, Ii_m
                    write(33,*) t, Ri_m
                    write(34,*) t, Si_m
					!###################################################
					! Calcula as medias das variaveis dinamicas
										
					!###################################################		
					! Escreve nos arquivos.
					
					!##############################################################$$$$$$$$$$$$$$$$
					call rk4_SIRS_grafo(rede, t, dt, rede%nodes, alp, lam, mu, Si, Ii, f_SIRS_Si, f_SIRS_Ii)	
					!##############################################################$$$$$$$$$$$$$$$$
					! Chama o metodo de Runge-Kutta de ordem 4
					
    				Si_m = 0.0_dp
	    			Ii_m = 0.0_dp
	    			do j2 = 1, rede%nodes
		    		   if(lista_de_clusters(j2) /= i_comp_gigante) cycle
			    	   Ii_m = Ii_m + Ii(j2)
			    	   Si_m = Si_m + Si(j2)
			    	   if(Si(j2) > 1.0_dp)then
			    	      stop "Si > 1"
			    	   elseif(Ii(j2) < 0.0_dp)then
			    	      stop "Ii negativo"
			    	   elseif((1.0_dp - Ii(j2)-Si(j2)) < 0.0_dp)then
			    	      stop "Ri negativo"
			    	   endif
				    enddo
				    Ii_m = 1.0_dp * Ii_m/comp_gigante
				    Si_m = 1.0_dp * Si_m/comp_gigante				
				    Ri_m = 1.0_dp -Si_m -Ii_m
				
                    write(32,*) t, abs(Ii_m-Ii_m0)/Ii_m0
											
			   if(Ii_m < 0.0_dp)then
			      write(*,*) "Ii_m", Ii_m
			      write(*,*) "No tempo t = ", t
			      stop
			   elseif(Si_m > 1.0_dp)then
			      write(*,*) "Si_m", Si_m
			      write(*,*) "No tempo t = ", t
			      stop
			   !elseif(Ri_m < 0.0_dp)then
               !   write(*,*) "Ri_m", Ri_m
               !   write(*,*) "No tempo t = ", t
               !   stop
                  ! exit loop_t
               !elseif(abs(Ii_m-Ii_m0)/Ii_m0 < toll)then
               !   repetiu = repetiu + 1
               !   if(repetiu > 10)then
               !   exit loop_t
               !   endif
			   endif
			
			
!			   if(Ii_m < 1.0_dp/(5.0_dp *rede%nodes))then
!                  Ii_m = 0.0_dp
!                  Si_m = 1.0_dp
!                  Ri_m = 0.0_dp
!                  exit loop_t
!               endif		
					
			   !##############################################################
               t = t + dt
               !##############################################################
			   ! Incrementa o tempo.			
               if( (t < (trec + 1.5_dp * dt)) .and. (t > (trec - 1.5_dp * dt)))then
            
                   write(14,*) t
                   
                   write(14,*) lam
                          
                   do j1 = 1, size(Ii)
                      write(14,*) Ii(j1)
                   enddo
                   rewind(14)

                   do j1 = 1, size(Ii)
                      write(15,*) 1.0_dp -Ii(j1)-Si(j1)
                   enddo
                   rewind(15)

                   do j1 = 1, size(Ii)
                      write(16,*) Si(j1)
                   enddo
                   rewind(16)
                                             
                   trec = trec + trec1 
               endif
               					
			   enddo loop_t
				
				close(31)
				close(32)
				close(33)
				!########################################################
	         !Rotinaprincipal
	            write(13,*) lam, Ii_m
	                 					
				lam = lam + dlam
				
				!########################################################
					
				!########################################################
				! Fecha os arquivos
		
			enddo	
			close(13)
			close(14)
			close(15)
		!close(14)


	!######################################################################
	
	contains
	
	function f_SIRS_Si(this, t, Si, Ii, n_sitios, alp, lam)
		use types
		use geraRede
		!use mod_tools
		class(grafo) :: this
		real(dp), intent(in) :: t, alp, lam
		integer, intent(in) :: n_sitios
		real(dp), intent(in) :: Si(n_sitios), Ii(n_sitios)
		real(dp) :: f_SIRS_Si(n_sitios)
		real(dp) :: sumIi
		!###############################
		integer :: i1, i2, j1
		!###############################
		
		! Sitio i
		do i1 = 1, this%nodes
		    if(lista_de_clusters(i1) /= i_comp_gigante) cycle				
			sumIi = 0.0d0
			! Vizinhos do sitio i.
			do i2 = this%aux(i1), this%aux(i1) + this%deg(i1) - 1
				j1 = this%listAdj(i2)
				sumIi = sumIi + Ii(j1)
			enddo
			f_SIRS_Si(i1) = alp * (1.0d0 -Si(i1) -Ii(i1)) 				
			f_SIRS_Si(i1) = f_SIRS_Si(i1) -lam * Si(i1) * sumIi
		enddo
	end function			

	!#################################################
	!	F_Ii
	!#################################################

	function f_SIRS_Ii(this, t, Si, Ii, n_sitios, mu, lam)
		use types
		use geraRede
		!use mod_tools
		
		class(grafo) :: this
		real(dp), intent(in) :: t, mu, lam
		integer, intent(in) :: n_sitios
		real(dp), intent(in) :: Si(n_sitios), Ii(n_sitios)
		real(dp) :: f_SIRS_Ii(n)
		real(dp) :: sumIi
		!###############################
		integer :: i1, i2, j1
		!###############################
		
		! Sitio i.
		do i1 = 1, this%nodes
			if(lista_de_clusters(i1) /= i_comp_gigante) cycle
			! Vizinhos do sitio i.
			sumIi = 0.0d0
			do i2 = this%aux(i1), this%aux(i1) + this%deg(i1) - 1
				j1 = this%listAdj(i2)					
				sumIi = sumIi +  Ii(j1)
			enddo
			f_SIRS_Ii(i1) = -1.0d0 * mu * Ii(i1)
			
			f_SIRS_Ii(i1) = f_SIRS_Ii(i1) + lam * Si(i1) * sumIi
		enddo
	end function


end program
