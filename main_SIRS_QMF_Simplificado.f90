
module qmf
	use geraRede
	use mod_tools
	use types
	implicit none
	
	real(dp), allocatable :: Ii(:), Ri(:)
	real(dp), allocatable :: k1_Ii(:), k1_Ri(:)
	real(dp), allocatable :: k2_Ii(:), k2_Ri(:)
	real(dp), allocatable :: k3_Ii(:), k3_Ri(:)
	real(dp), allocatable :: k4_Ii(:), k4_Ri(:)
	real(dp) :: t, dt
	real(dp), allocatable :: v1(:)
	
	contains
	
	
		!#################################################
		!	Modelo SIRS
		!#################################################		
!=======================================================================
   subroutine condicao_inicial_QMF(this, p0, t0)
      class(grafo), intent(in) :: this      
      real(dp), intent(in) :: p0
      real(dp) :: t0
      integer :: l1
      !=================================================================
      t = t0
      do l1 = 1, this%nodes
         Ii(l1) = p0
         Ri(l1) = 0.01_dp
      enddo   
      !=================================================================
   end subroutine
!=======================================================================
    subroutine aloca_listas_QMF(this)
       class(grafo) :: this
       
       
       if(allocated(k1_Ii))then
          deallocate(Ii)
          deallocate(Ri)
          
          deallocate(k1_Ii)
          deallocate(k2_Ii)
          deallocate(k3_Ii)
          deallocate(k4_Ii)
          deallocate(k1_Ri)
          deallocate(k2_Ri)
          deallocate(k3_Ri)
          deallocate(k4_Ri)
       endif          
          
       allocate(Ii(this%nodes))
       allocate(k1_Ii(this%nodes))
       allocate(k2_Ii(this%nodes))   
       allocate(k3_Ii(this%nodes))
       allocate(k4_Ii(this%nodes))       

       allocate(Ri(this%nodes))   
       allocate(k1_Ri(this%nodes))
       allocate(k2_Ri(this%nodes))   
       allocate(k3_Ri(this%nodes))
       allocate(k4_Ri(this%nodes))
    
    end subroutine
					
    subroutine rk4_SIRS_QMF_grafo(this, t, dt, alp, lam, mu)
		    class(grafo) :: this
         	real(dp) :: t, dt
        	integer :: n1
        	real(dp) :: alp, lam, mu		 
  		!####################################################################################
		n1 = this%nodes
		
		!####################################################################################	
		!            f_Ri(this, t, Ri_1, Ii_1, n_sitios, alp, mu)
		k1_Ri = dt * f_Ri(this, t, Ri, Ii, n1, alp, mu)
		
		!            f_Ii(this, t, Ri, Ii, n_sitios, mu, lam)
		k1_Ii = dt * f_Ii(this, t, Ri, Ii, n1, lam, mu)

		
		!####################################################################################
		
		!            f_Ri(this, t              , Ri_1               , Ii_1               , n_sitios, alp, mu)
		k2_Ri = dt * f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, n1, alp, mu)		
		
		!            f_Ii(this, t,               Ri                 , Ii                 , n_sitios, mu, lam)
		k2_Ii = dt * f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k1_Ri, Ii + 0.5_dp * k1_Ii, n1, lam, mu)


        !####################################################################################
		!            f_Ri(this, t              , Ri_1               , Ii_1               , n_sitios, alp, mu)        
		k3_Ri = dt * f_Ri(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, n1, alp, mu)

        !            f_Ii(this, t,               Ri                 , Ii                 , n_sitios, mu, lam)
		k3_Ii = dt * f_Ii(this, t + 0.5_dp * dt, Ri + 0.5_dp * k2_Ri, Ii + 0.5_dp * k2_Ii, n1, lam, mu)


       !####################################################################################
		!            f_Ri(this, t     , Ri_1      , Ii_1      , n_sitios, alp, mu)       
		k4_Ri = dt * f_Ri(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, n1, alp, mu)
		
		!            f_Ii(this, t,      Ri        , Ii        , n_sitios, mu, lam)
		k4_Ii = dt * f_Ii(this, t + dt, Ri + k3_Ri, Ii + k3_Ii, n1, lam, mu)



		!####################################################################################	
		Ri = Ri + (1.0_dp/6.0_dp) * (k1_Ri + 2.0_dp * k2_Ri + 2.0_dp * k3_Ri + k4_Ri)  		
  		
		
		Ii = Ii + (1.0_dp/6.0_dp) * (k1_Ii + 2.0_dp * k2_Ii + 2.0_dp * k3_Ii + k4_Ii)  		
  		!####################################################################################
  		! F_Ii
  	  	        		    
    end subroutine		


	function f_Ri(this, t, Ri, Ii, n_sitios, alp, mu)

		!use types
		!use geraRede

		class(grafo) :: this
		real(dp), intent(in) :: t, alp, mu
		integer, intent(in) :: n_sitios
		real(dp), intent(in) :: Ri(n_sitios), Ii(n_sitios)
		real(dp) :: f_Ri(n_sitios)
		!###############################
		integer :: l1, l12, l2
		!###############################
		
		! Sitio i
		do l1 = 1, this%nodes
		    if(lista_de_clusters(l1) /= i_comp_gigante) cycle
			f_Ri(l1) = -alp * Ri(l1) + mu * Ii(l1)	
		enddo
		
	end function			

	!#################################################
	!	F_Ii
	!#################################################

	function f_Ii(this, t, Ri, Ii, n_sitios, lam, mu)
		!use types
		!use geraRede
		
		
		class(grafo) :: this
		real(dp), intent(in) :: lam,  mu, t 
		integer, intent(in) :: n_sitios
		real(dp), intent(in) :: Ii(n_sitios), Ri(n_sitios)
		real(dp) :: f_Ii(n_sitios)
		real(dp) :: sumIi
		!###############################
		integer :: l1, l12, l2, l3
		!###############################
		
		! Sitio i.
		do l1 = 1, this%nodes			
			if(lista_de_clusters(l1) /= i_comp_gigante) cycle
			f_Ii(l1) = -mu * Ii(l1)
			sumIi = 0.0_dp
			do l12 = this%aux(l1), this%aux(l1) + this%deg(l1) - 1
				l2 = this%listAdj(l12)
				sumIi = sumIi +  Ii(l2)
			enddo
			
			f_Ii(l1) = f_Ii(l1) + lam * (1.0_dp - Ri(l1) - Ii(l1) ) * sumIi
			
		enddo
	end function
     	
end module

module limiaresCampoMedio
   use geraRede
   use mod_rndgen
   implicit none
   
   integer, private, parameter :: dp = kind(0.0d0)
   
   real(dp), allocatable :: x(:), y(:)
   real(dp) :: Y4
   
   contains
   
   subroutine aloca_listas_e_matrizes(this)
      class(grafo) :: this
      if(allocated(x)) deallocate(x)
      allocate(x(this%nodes))

      if(allocated(y)) deallocate(y)
      allocate(y(this%nodes))
      
   end subroutine
   
   subroutine limiarQMF(this, lambda1, tole)
      class(grafo), intent(in) :: this
      real(dp), intent(inout) :: lambda1
      real(dp) :: tole
      !=================================================================
      real(dp) :: autovalor0, autovalor_old
      real(dp) :: Y4_0
      !=================================================================
      integer :: j1, j2, j3, j4
      integer(kind = 8) ::  j12, j23
      real(dp) :: xp, yp
      integer :: ipx, ipy
      type(rndgen) :: geni
      real(dp) :: soma
      integer :: vizim
      real(dp) :: a_l, a_2l
      integer :: semente2
      real(dp) :: erro               
      real(dp), allocatable :: v1(:)
      real(dp), allocatable :: x_prov(:)
      real(dp) :: x_norm, y_norm
      real(dp) :: theta
      integer :: num_cont
      integer(kind=8) :: cont
      
      call aloca_listas_e_matrizes(this)
      !=================================================================
      
      !=================================================================      
      call geni%init(semente2)     
      !=================================================================
	  ! Condicao inicial				
	  !x = 1.0d0
	  !xp = 1.0d0
      !ipx = 1
      
      !xp = 0.0d0
      !ipx = 1
      
      if(allocated(x_prov)) deallocate(x_prov)
      allocate(x_prov(this%nodes))
      
      do j1 = 1, this%nodes
         y(j1) = 1.0_dp * this%deg(j1) !geni%rnd() * 1.0d-3/(1.0d0 * this%nodes)
      enddo
      !=================================================================
      y_norm = sqrt(sum(y**2.0d0))
      cont = 0
      num_cont = 2
      !=================================================================  
!=======================================================================      
iter: do while(.True.)	
         !==============================================================
	     !Aqui comeca o algoritmo
         !==============================================================
         ! Y = A X
         !==============================================================
         y = y/y_norm         
         x_prov = y
         !==============================================================         
         do j3 = 1, num_cont
            x = y 
            !===========================================================         
            do j1 = 1, this%nodes
               !========================================================
               soma = 0.0_dp
               do j12 = this%aux(j1), this%aux(j1) + this%deg(j1) - 1
                  j2 = this%listAdj(j12)
                  soma = soma + x(j2)
               enddo
               y(j1) = soma
            enddo
         enddo

         y_norm = sqrt(sum(y**2.0d0)) 
         x_norm = sqrt(sum(x**2.0d0))
         
         Y4 = sum(((x/x_norm + y/y_norm)/2.0d0)**4.0d0)
                  
         autovalor0 = (dot_product(x_prov, y))**(1.0d0 /(1.0d0*num_cont))
                 
         erro = max( abs(Y4 - Y4_0), abs(autovalor0 - autovalor_old))
         
         Y4_0 = Y4
         autovalor_old = autovalor0
         
         if(erro < tole)then
            write(*,*) "LEV Jacobiana = ", autovalor0
            lambda1 = 1.0_dp/autovalor0
            write(*,*) "lambda_C = ", 1.0_dp/autovalor0
            !-----------------------------------------------------------
			!if(allocated(v1)) deallocate(v1)
			!allocate(v1(size(x)))
			
			!v1 = x
			
			!v1 = v1/(sum(v1**2.0_dp))**0.5_dp
			
			!Y4 = sum(v1**4.0_dp)

			!deallocate(v1)            
            !-----------------------------------------------------------                        
            exit iter
         endif
         
         !==============================================================
		 ! O ipy achado apos o produto matricial eh o indice ipx da
		 ! primeira maior componente do vetor x.
		 ! xp neste caso agora = 1.0d0, pois y ja foi normalizado
		 ! na linha anterior ao comando 'if(erro < tole)then'
		 !ipx = ipy
		 !xp = 1.0d0
		 !==============================================================
         !j1 = j1 + 1
         !==============================================================
      enddo iter   
      
      deallocate(x)
      deallocate(y)     
   end subroutine
   
end module

program main
   use geraRede
   use mod_rndgen
   use mod_tools
   use qmf      
   use limiaresCampoMedio
   implicit none
!#######################################################################   
   type(grafo_PL_UCM) :: rede
!#######################################################################   
   real(dp) :: dlamb
   real(dp) :: lamb0, lambdaf
   real(dp) :: lamb
   integer :: nlamb = 1000
   real(dp), parameter :: mu = 1.0_dp
   real(dp) :: alp
   character(len=10) :: alp_char2
   real(dp) :: rho, rho0
   real(dp) :: tol, tole
!#######################################################################
   !type(rndgen) :: gen1
   integer, allocatable :: seed(:)
   integer :: seed_i, seed_f
   integer :: resto
   integer :: seed1
   type(rndgen) :: ger_inic
!#######################################################################
   integer :: tam_rede
   real(dp) :: gama_exp
   integer :: grau_min
   real(dp) :: grau_max      
!#######################################################################   
   character(len=500) :: t_vs_Im
   character(len=500) :: lamb_vs_Im
   character(len=500) :: lamb_vs_Xi 
   character(len=1000) :: caminho
   character(len=500) :: arquivo_rede_real
   character(len=7) :: tam_char
   character(len=5) :: gama_char
   character(len=5) :: indice
!#######################################################################   
   integer :: i1, i12, i2, i3, i4, i5, i6, i7, i8, i9 !, j1, j2, j3
   logical :: T_vs
!#######################################################################
   integer :: sumDeg2
   real(dp) :: t0, tf
   real(dp) :: dt_m
   real(dp), allocatable :: Ap_list(:)
   real(dp), allocatable :: P_grau(:)
   character(len=500) ::arq_1
   integer :: ntempo
   
   integer :: n_it, n_sub, per_conv   
   real(dp) :: p0
   character(len=300) :: cwd, resultados, tipoCorte
   character(len=1000) :: local
   character(len=1000) :: nomeArquivo
   character(len=20) :: buffer
   !##############################################################
   
   real(dp) :: qm, q2m

   integer :: nargus, ind_lamb   
   character(len=3) :: char_ind_lamb
   integer :: niter_abs   

   real(dp) :: lbdC_PQMF
   integer :: ind_amostra
   integer :: divisor
   character(len=10) :: lamb_char
   real(dp) ::lamb_C
   real(dp) :: tempoEscrita, tempoEscrita0
   real(dp) :: lamb_Copia
   logical :: existe
   integer :: st
   logical :: usouCopia
   logical :: teveLeitura
   real(dp) :: lambda_Ultimo_Index
   character(len=10) :: teoriaCM
   logical :: aberta
   integer :: int_soCalculaIPR
   logical :: soCalculaIPR

   real(dp) :: Y4_old
   real(dp) :: sumIi2
   logical :: Ii_neg
   logical :: Ii_gt_1

   integer(kind=8) :: so_gasta   
!#######################################################################   
   seed1=947361823
!#######################################################################   
   if(allocated(seed)) deallocate(seed)
   !--------------------------------------------------------------------
   seed_i = 1
   !--------------------------------------------------------------------
   seed_f = seed_i + 9
   resto = seed_i - 1
   !--------------------------------------------------------------------
   allocate(seed(seed_i:seed_f))
!#######################################################################
   teoriaCM = trim(adjustl('QMF'))
   !------------------------------------------------------
   !tipoCorte ='_Rigido'
   tipoCorte ='_2sqrtN'
   !tipoCorte ='_sqrtN'
   !------------------------------------------------------
   resultados = 'Rst_QMF_Corte'//trim(adjustl(tipoCorte))
   !------------------------------------------------------
   resultados = trim(adjustl(resultados))
   !------------------------------------------------------
   call system('mkdir -p '//resultados)
   !------------------------------------------------------
   local = trim(adjustl(resultados))//"/"
   !------------------------------------------------------
   call entradaArgumentos()
   !------------------------------------------------------
!#######################################################################
   !====================================================================
   !  Se der ruim no dt, ele eh dividido por dois.
   !====================================================================   
   dt = 1.0d-1 
   !====================================================================
   ! A principio, t0 = 0, mas, se houver um estado salvo, muda
   ! para o t salvo no arquivo.
   !====================================================================
   t0 = 0.0_dp
   tf = 10000.0_dp
   
   tempoEscrita0 = 10.0_dp !tf/100.0d0
   
   ntempo = int( (tf- t0)/dt )
   tole = 1d-7
   per_conv = int( 10.0_dp/dt )
!#######################################################################

!#######################################################################
   call ger_inic%init(seed1)

   if( resto == 0)then
      i2 = 1
      do i1 = 1, 1000
         if(mod(i1,100) > 0) cycle
         seed(i2)  = ger_inic%int(100000000,999999999)
         write(*,*) i1, seed(i2)
         i2 = i2+1
      enddo
   elseif( resto > 0)then
      i2 = seed_i
      do i1 = 1, 1000
         if(mod(i1,100) /= resto )then !Antes, eu dava um cycle. Agora eh assim
            so_gasta = ger_inic%int(100000000,999999999)
         else 
            seed(i2)  = ger_inic%int(100000000,999999999)
            write(*,*) i1, seed(i2)
            i2 = i2+1
         endif      
      enddo
   endif
   
   local = trim(adjustl(trim(adjustl(local))//'tam_'//trim(adjustl(tam_char))//'/'))
   call system('mkdir -p '//trim(adjustl(local)) )
   
   local = trim(adjustl(trim(adjustl(local))//'gam_'//trim(adjustl(gama_char))//'/'))
   call system('mkdir -p '//trim(adjustl(local)) )   

   local = trim(adjustl(trim(adjustl(local))//'ams_'//trim(adjustl(indice))//'/'))
   
   call system('mkdir -p '//trim(adjustl(local)) )


!#######################################################################
!				Inicia grafo
!#######################################################################
 !  arquivo_rede_real='s00088.s838.net.edg'
 !  call rede%RedeReal(arquivo_rede_real, 111)
 !  close(111)
!#######################################################################
   call criaRedeEClassificaClusters(rede, tam_rede, grau_min, grau_max, gama_exp, seed(ind_amostra))

   !====================================================================
   write(*,*) "######################Dados da Rede######################"
   write(*,*) ""
   write(*,*) "Tamanho da rede ", rede%nodes, "."
   write(*,*) ""
   write(*,*) "Fracao correspondente aa componente gigante ", 100.0 * comp_gigante/rede%nodes,"%", "."
   write(*,*) ""
   write(*,*) "Expoente da distribuicao de graus da rede ", gama_exp, "."
   write(*,*) ""
   write(*,*) "Grau minimo ", rede%degMin, ".", " Grau maximo ", rede%degMax, "."
   write(*,*) ""
!#######################################################################
   write(alp_char2, '(f9.3)') alp

   local = trim(adjustl(trim(adjustl(local))//'alp_'//trim(adjustl(alp_char2))//trim(adjustl(tipoCorte))//'/'))

   call system('mkdir -p '//trim(adjustl(local)) )

   call kNN_e_clustering(rede)
!#######################################################################
   write(*,*) "######################Dados temporais######################"
   write(*,*) ""
   write(*,*) "Instante inicial ", t0, "."	
   write(*,*) ""
   write(*,*) "Resolucao ", dt, "."
   write(*,*) ""
   write(*,*) "Instante final ", tf, "."
   write(*,*) ""
   write(*,*) "Quantidade de pontos ", ntempo, "."
   write(*,*) ""
   write(*,*) "Tolerancia considerada no criterio de conv. ", tole, "."
!#######################################################################
   write(*,*) "###############Parametros Epidemicos#####################"
   write(*,*) "Probabilidade de sitio estar infectado ", p0, "."
   write(*,*) ""
   write(*,*) "Taxa de recuperacao ", mu, "."	
   write(*,*) ""
   write(*,*) "Taxa de enfraquecimento imunologico ", alp, "."
   write(*,*) ""
   write(*,*) "#########################################################"
!#######################################################################
   call calcula_P_grau(rede)
!#######################################################################   

!#######################################################################

!#######################################################################
   
   !====================================================================
   call aloca_listas_QMF(rede)
   !====================================================================
   teoriaCM = trim(adjustl(teoriaCM))
   usouCopia = .False.
   !====================================================================
   call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)
   !====================================================================
      
   !====================================================================
   if( nargus == 9)then
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//trim(adjustl(teoriaCM))         
      !=================================================================
      open(334, file=trim(adjustl(nomeArquivo))//'.dat', access='append', status='unknown')
      !=================================================================  
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_'))//trim(adjustl(teoriaCM))         
      !=================================================================
      open(335, file=trim(adjustl(nomeArquivo))//'.dat', access='append', status='unknown')

      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_'))//trim(adjustl(teoriaCM))         
      !=================================================================
      !====================================================================
      open(337, file=trim(adjustl(nomeArquivo))//'.dat', access='append', status='unknown')   
      !==================================================================== 
      
   elseif( nargus == 10)then
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(char_ind_lamb))       
      !====================================================================
      open(334, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
      !=================================================================
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(char_ind_lamb)) 
      !====================================================================
      open(335, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
      !=================================================================             
      !=================================================================
      nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(char_ind_lamb)) 
      !=================================================================             
      !====================================================================
      open(337, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')   
      !==================================================================== 
   endif
   stop "Calculou limiar. Nao vai rodar dinamica."
  
   
llbd: do while( lamb <= lambdaf )
      !============================================================
      ! tempo
      !============================================================
      write(lamb_char, '(f10.6)') lamb     
      n_it = 0
      n_sub = 0
      !============================================================
      if( nargus == 10)then
         nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_Y4_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
         !=================================================================
         inquire(unit=333, opened = aberta)
         if(aberta)close(333)
         open(333, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')   
         !=================================================================
      
         nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rho_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
         !=================================================================
         inquire(unit=336, opened = aberta)
         if(aberta)close(336)
         open(336, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')   
         !=================================================================
      elseif( nargus == 9)then
         nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_Y4_'))//trim(adjustl(teoriaCM)) 
         !=================================================================
         inquire(unit=333, opened = aberta)
         if(aberta)close(333)
         open(333, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')   
         !=================================================================
      
         nomeArquivo = trim(adjustl(local))//trim(adjustl('t_vs_rho_'))//trim(adjustl(teoriaCM)) 
         !=================================================================
         inquire(unit=336, opened = aberta)
         if(aberta)close(336)
         open(336, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')   
         !=================================================================
      endif           
           

      
      !=================================================================     
      if( .not. usouCopia)then 
         rho0 = p0
         
         rho = rho0
         
         Y4_old = 1.0_dp/(1.0_dp * comp_gigante)
         
         Y4 = Y4_old
         
         call condicao_inicial_QMF(rede, p0, 0.0_dp)              
         
         tempoEscrita = tempoEscrita0
      !=================================================================
      else
         usouCopia = .False.
      endif
      !=================================================================     
lt:   do while(t <= tf)
         
         call tempoDeEscrever()
         
         write(333, *) t, Y4
         
         write(336, *) t, rho
         
         !##############################################################
         call rk4_SIRS_QMF_grafo(rede, t, dt, alp, lamb, mu)
         !##############################################################
         
         rho = 0.0_dp
         
         Ii_neg = .False.
         Ii_gt_1 = .False.
         do i1 = 1, size(Ii)   ! Calcula rho
            if(lista_de_clusters(i1) /= i_comp_gigante) cycle
            if( Ii(i1) < 0.0_dp)then
               Ii_neg = .True.
               exit
            elseif( Ii(i1) > 1.0_dp )then
               Ii_gt_1 = .True.
               exit
            endif
            rho = rho + Ii(i1)
         enddo
         
         rho  = rho/( 1.0_dp * comp_gigante )
         
         sumIi2 = (sum(Ii**2.0_dp))**0.5_dp
         
		 Y4 = sum( (Ii/sumIi2)**4.0_dp )
         
         
         !##############################################################
         t = t + dt
         !##############################################################
         
         if( abs(rho - rho0)/rho0 < tole )then
            n_it = n_it + 1
            if( n_it >= int(2.0d0/dt) )then
               write(*,*) "lbd = ", lamb, " rho = ", rho
               exit lt
            endif
         else
            n_it = 0
         endif

         !if( ( abs(rho - rho0)/rho0 > 5.0d0 * tole ) .and. ( rho < 1.0_dp/(10.0_dp * rede%nodes) ) )then
         !   n_sub = n_sub + 1
         !   if ( n_sub >= int(5.0d0/dt) )then
         !      write(*,*) "Abs subcr com lambda = ", lamb, " e rho = ", rho, "."              
         !      exit lt
         !   endif
         !else
         !   n_sub = 0
         !endif
         
         !==============================================================
         !   Se deu ruim por causa de dt grande, a gente pega uma
         !   configuracao passada e abaixa o dt
         !   
         !==============================================================
                           
         if( (Ii_gt_1 == .True.) .or. (Ii_neg == .True.) )then                 
            write(*,*) "Dev. err. num. atliz. dt de ", dt, " para ", dt/2.0d0            
            usouCopia = .False.
            
            call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)
            
            dt = dt/2.0d0
            
            if( .not. usouCopia)then 
               rho0 = p0
         
               call condicao_inicial_QMF(rede, p0, 0.0_dp)                  
         
               tempoEscrita = tempoEscrita0
               !=========================================================
            else   
               usouCopia = .False.
            endif            
            
            cycle lt
         endif   
         !==============================================================
         !   Antes faziamos assim :/ ...
         !   Isso devia ir pra um museu!
         !==============================================================
         !stop "Ajuste dt para que rho deixe de ser negativo."
         !if( rho > 1.0_dp) stop "Ajuste dt para que rho deixe de ser > 1.0"
         !##############################################################
         rho0 = rho
         Y4_old = Y4
                  
      enddo lt
      !=================================================================
      close(333)
      close(336)
      !=================================================================
      if( n_it >= int(2.0d0/dt) )then
         write(334,*) lamb, rho   
         write(335,*) lamb, t
         write(337,*) lamb, Y4
      endif    
      !=================================================================
      lamb = lamb + dlamb
      !================================================================= 
      call system('rm -r '//trim(adjustl(local))//trim(adjustl('Copia*')) )
            
      !=================================================================       
        
      close(333)
   enddo llbd
   close(334)
   close(335)
   close(337)
!#######################################################################   

      if(allocated(k1_Ii))then
         deallocate(k1_Ii) 
         deallocate(k2_Ii)
         deallocate(k3_Ii)
         deallocate(k4_Ii)

         deallocate(k1_Ri) 
         deallocate(k2_Ri)
         deallocate(k3_Ri)
         deallocate(k4_Ri)
      endif            
	
	contains

   !====================================================================
      subroutine leAtehOfimUnidimensional(label, nomeArquivo, var1, teveLeitura)      
         integer, intent(in) :: label
         character(len=*) :: nomeArquivo
         real(dp), intent(out) :: var1
         logical, intent(inout) :: teveLeitura
         integer :: st, res
         
         inquire( file=trim(adjustl(nomeArquivo)), exist=res ) 
         if(res)then
            open(label, file=trim(adjustl(nomeArquivo)), status='old')
            write(*,*) "Arquivo jah existe"
            do
               read(label, *, iostat = st) var1
               write(*,*) st
               if( st /= 0) exit
               teveLeitura = .True.
            enddo
         else
            open(label, file=trim(adjustl(nomeArquivo)), status='new')
            write(*,*) "Arquivo teve que ser criado"
         endif
   !====================================================================         
      end subroutine


   !====================================================================
      subroutine leAtehOfimBidimensional(label, nomeArquivo, var1, var2, teveLeitura)      
         integer, intent(in) :: label
         character(len=*) :: nomeArquivo
         real(dp), intent(out) :: var1, var2
         logical, intent(inout) :: teveLeitura
         integer :: st, res
         
         inquire( file=trim(adjustl(nomeArquivo)), exist=res ) 
         if(res)then
            open(label, file=trim(adjustl(nomeArquivo)), status='old')
            write(*,*) "Arquivo jah existe"
            do
               read(label, *, iostat = st) var1, var2               
               if( st /= 0) exit
               teveLeitura = .True.
            enddo
         else
            open(label, file=trim(adjustl(nomeArquivo)), status='new')
            write(*,*) "Arquivo teve que ser criado"
         endif
   !====================================================================         
      end subroutine

 
      subroutine criaRedeEClassificaClusters(this, tam_rede1, grau_min1, grau_max1, gama_exp1, seme1)
        type(grafo_PL_UCM) :: this
        integer :: tam_rede1
        integer :: grau_min1
        real(dp) :: grau_max1
        real(dp) :: gama_exp1
        integer :: seme1
!#######################################################################   
        call this%iniciaGrafo(tam_rede1)
!#######################################################################
        call this%inicia(grau_min1, grau_max1, gama_exp1, seme1)
!#######################################################################
        call this%liga(seme1, .False.) 
!#######################################################################        
        call sub_classifica_clusters(this,.False., 000, 'sem_arquivo.dat') 
!#######################################################################
      end subroutine
!#######################################################################

      subroutine entradaArgumentos()
         nargus = iargc()

         if(nargus == 9)then
            !#############################
            ! Amostra
            !#############################
            write(*,*) "Recebeu 8 arqumentos"
            call getarg(1, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) ind_amostra
            write(*,*) "O indice da amostra eh: ", ind_amostra
            !#############################
            !	Tamanho da rede
            !#############################
            call getarg(2, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) tam_rede

            !#############################
            ! Tamanho	da rede
            !#############################
            call getarg(3, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grau_min

            !#############################
            ! Expoente Gama
            !#############################
            call getarg(4, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) gama_exp

            !#############################
            ! Lambda0
            !#############################
            call getarg(5, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lamb0

            !#############################
            ! Divisor que fornece dlambda
            !#############################
            call getarg(6, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) divisor
      
            write(*,*) "O valor do divisor de 0.0125 eh: ", divisor


            !#############################
            ! Lambdaf
            !#############################
            call getarg(7, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lambdaf
            write(*,*) "O valor de lambdaf eh: ", lambdaf
            !#############################
            ! Alfa
            !#############################
            call getarg(8, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) alp

            call getarg(9, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) int_soCalculaIPR
      
            if(int_soCalculaIPR == 0)then
               soCalculaIPR = .False.
            elseif(int_soCalculaIPR == 1)then
               soCalculaIPR = .True.
            else
               stop "Nao foi possivel identificar o valor de int_soCalculaIPR"
            endif                
            
         elseif(nargus == 10)then
            write(*,*) "Recebeu 9 arqumentos"
            !#############################
            ! Amostra
            !#############################
            call getarg(1, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) ind_amostra
            write(*,*) "O indice da amostra eh: ", ind_amostra
            !#############################
            !	Tamanho da rede
            !#############################
            call getarg(2, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) tam_rede

            !#############################
            ! Tamanho	da rede
            !#############################
            call getarg(3, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) grau_min

            !#############################
            ! Expoente Gama
            !#############################
            call getarg(4, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) gama_exp

            !#############################
            ! Lambda0
            !#############################
            call getarg(5, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lamb0

            !#############################
            ! Divisor que fornece dlambda
            !#############################
            call getarg(6, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) divisor
      
            write(*,*) "O valor do divisor de 0.0125 eh: ", divisor


            !#############################
            ! Lambdaf
            !#############################
            call getarg(7, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) lambdaf

            !#############################
            ! Alfa
            !#############################
            call getarg(8, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) alp
      
            !#############################
            ! Indice do ponto
            ! Quando fizermos
            ! paralelizacao burra.
            !#############################
            call getarg(9, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) ind_lamb   
            write(char_ind_lamb, '(I0)') ind_lamb
        
            char_ind_lamb = trim(adjustl(char_ind_lamb)) 

            call getarg(10, buffer)
            buffer = trim(adjustl(buffer))
            read(buffer,*) int_soCalculaIPR
      
            if(int_soCalculaIPR == 0)then
               soCalculaIPR = .False.
            elseif(int_soCalculaIPR == 1)then
               soCalculaIPR = .True.
            else
               stop "Nao foi possivel identificar o valor de int_soCalculaIPR"
            endif                    
         else
            stop "Forneca dados no arquivo 'sirs_estocastico_cluster.sh' "
         endif

!#######################################################################
         if( trim(adjustl(resultados)) == 'Rst_QMF_Corte_Rigido' )then
            call acha_cutoff_rigido(grau_min, gama_exp, tam_rede)
         elseif(trim(adjustl(resultados)) == 'Rst_QMF_Corte_sqrtN' )then
            grau_max = (1.0_dp * tam_rede)**(0.5_dp)
         elseif(trim(adjustl(resultados)) == 'Rst_QMF_Corte_2sqrtN' )then
            grau_max = 2.0_dp * (1.0_dp * tam_rede)**(0.5_dp)
         endif

         dlamb = 0.0125_dp/(1.0_dp * divisor)
!#######################################################################
         p0 = 10.0_dp/(1.0_dp * tam_rede) !1.0d-2
!#######################################################################
!  Indice da semente
!#######################################################################
         write(gama_char,'(f5.2)') gama_exp
!#######################################################################
         write(indice,'(I0)') ind_amostra 
!#######################################################################   
         if(tam_rede == 10)then
            tam_char = '10'
         elseif(tam_rede == 100)then
            tam_char = '100'
         elseif(tam_rede == 1000)then
            tam_char = '1k'
         elseif(tam_rede == 3000)then
            tam_char = '3k'
         elseif(tam_rede == 10000)then
            tam_char = '10k'
         elseif(tam_rede == 30000)then
            tam_char = '30k'
         elseif(tam_rede == 100000)then
            tam_char = '100k'
         elseif(tam_rede == 300000)then
            tam_char = '300k'
         elseif(tam_rede == 1000000)then
            tam_char = '1M'
         elseif(tam_rede == 3000000)then
            tam_char = '3M'
         elseif(tam_rede == 10000000)then
            tam_char = '10M'
         elseif(tam_rede == 30000000)then
            tam_char = '30M'
         elseif(tam_rede == 100000000)then
            tam_char = '100M'
         else
            stop 'Escolha um tamanho de rede dentro do catalogo'
         endif         
      end subroutine

!=======================================================================
   subroutine calcula_P_grau(this)      
      class(grafo) :: this
      integer :: i1
!#######################################################################
      if(allocated(P_grau)) deallocate(P_grau)
      allocate(P_grau(this%degMin:this%degMax))
      P_grau = 0_dp
      do i1 = 1, this%nodes
         P_grau(this%deg(i1)) = P_grau(this%deg(i1)) + 1.0_dp
      enddo
      P_grau = P_grau/(1.0_dp * this%nodes)
      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_P_grau_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'.dat'
      open(800, file=trim(adjustl(arq_1)), status='unknown')
      do i1 = this%degMin, this%degMax
         if(P_grau(i1) == 0.0_dp) cycle
         write(800,*) i1, P_grau(i1)
      enddo
      close(800)
      deallocate(P_grau)
!#######################################################################      
   end subroutine
!=======================================================================
   subroutine kNN_e_clustering(this)
      class(grafo) :: this

      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_knn_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'.dat'
      call calcula_k_nn(this, .True., 800, trim(adjustl(arq_1)))
      close(800)

      arq_1 = trim(adjustl(local))//trim(adjustl('grau_vs_clustering_'))//'tam_'//trim(adjustl(tam_char))//'_gam_'//trim(adjustl(gama_char))//'.dat'
      call clustering(this,.True., 800, trim(adjustl(arq_1)))
      close(800)
   end subroutine
!=======================================================================   

  recursive subroutine TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)
      real(dp) :: rho1
      logical ::  existe2
      integer ::  st1, st2
      integer :: ind_arquivo
      character(len=*) :: teoriaCM
      logical, intent(inout) :: usouCopia
      logical :: taAberta
      
      teoriaCM = trim(adjustl(teoriaCM))
      !=================================================================
      ! A principio, dizemos que nao usara a Copia.
      ! Apenas com os testes abaixo isso mudara. Ou nao.
      !=================================================================
      
      existe2 = .False.
      st1 = 0; st2 = 0
      
      existe = .False.
      
      if(nargus == 10)then

         nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_lambda'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))

         inquire( file=trim(adjustl(nomeArquivo))//'.dat', exist=existe )
                         
         if(existe)then
            !===========================================================
            usouCopia = .True.
            !===========================================================
            inquire(unit=776, opened=taAberta)
            if(taAberta) close(776)
            !===========================================================
            open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            read(776,*, iostat=st) lamb0
            close(776)
            !===========================================================
            if( st /= 0)then
               usouCopia = .False.
               call system('rm -r '//trim(adjustl(local))//'Copia*')               
               write(*,*) "Rotina deletou Copia_lambda.dat e se chamou sozinha"
               call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)
            endif
            !===========================================================
         else
            !===============================================================
            nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//trim(adjustl(teoriaCM))
            
            inquire(file=trim(adjustl(nomeArquivo))//'.dat', exist = existe2)
            teveLeitura = .False.
            !===========================================================                     
            if(existe2)then
               !========================================================
               inquire(unit=334, opened=taAberta)
               if(taAberta) close(334)
               
               open(334, file = trim(adjustl(nomeArquivo))//'.dat', status='old')
               read(334, *, iostat=st) lamb0, rho1
               close(334)
               
               if( st == 0 )then			                    
                  if( rho1 == 0.0_dp )then
                     nomeArquivo = trim(adjustl(local))//trim(adjustl('lamb_vs_Y4_'))//trim(adjustl(teoriaCM))                     
                     inquire(file=trim(adjustl(nomeArquivo))//'.dat', exist=existe2)
                     
                     st = -1
                     if(existe2)then
                        inquire(unit=444, opened=taAberta)
                        if(taAberta) close(444)                     
                        
                        open(444, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
                        read(444, *, iostat=st) lamb0, Y4
                        if(st == 0) teveLeitura = .True.
                        close(444)
                     endif
                     
                     if( st /= 0 )then                                                
			            if(soCalculaIPR)then
                           !============================================  
                           lamb0 = 0.0_dp                    ! Entrou de um jeito
                           call limiarQMF(rede, lamb0, tole)
                           !=============================================          
		   	               open(444, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
			               write(444,*) lamb0, Y4
			               close(444)
			               stop
			            endif
			         endif			      
			      endif
			   endif
			   close(334)               
               
               !========================================================
               !call leAtehOfimBidimensional(334, trim(adjustl(nomeArquivo))//'.dat', lamb0, rho1, teveLeitura)
               !========================================================
               !close(334)
               !========================================================
            endif

            if(.not.teveLeitura)then
               !============================================================
               nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//trim(adjustl(teoriaCM))//'_lamb_index_0'
               teveLeitura = .False.
               
               !========================================================
               inquire(unit=334, opened=taAberta)
               if(taAberta) close(334) 
               !========================================================
               open(334, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
               !========================================================
               !========================================================
               ! Damos voto de confianca
               !========================================================
               teveLeitura = .True.
               read(334, *, iostat=st2) lamb0, rho1
               close(334)
               !========================================================
               if(st2 /=0 ) teveLeitura = .False.    
               !========================================================               
               if( .not. teveLeitura)then        
                  !=====================================================  
                  lamb0 = 0.0_dp                    ! Entrou de um jeito
                  call limiarQMF(rede, lamb0, tole)
                  ! Saiu de outro.
                  !=====================================================
                  nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_'))//trim(adjustl(teoriaCM))//trim(adjustl('_lamb_index_0'))
                  
                  inquire(unit=444, opened=taAberta)
                  if(taAberta) close(444)    
                  
			      open(444, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
			      write(444,*) lamb0, Y4
   			      close(444)			      
                  !===================================================== 

                  nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//trim(adjustl(teoriaCM))//'_lamb_index_0'
                  open(unit=334, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
                  !=====================================================
                  write(334,*) lamb0, 0.0d0
                  close(334)
                  !=====================================================            
                  nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_'))//trim(adjustl(teoriaCM))//'_lamb_index_0.dat' 
                  !=====================================================
                  inquire(unit=335, opened=taAberta)
                  if(taAberta) close(335)
                  !=====================================================
                  open(335, file=trim(adjustl(nomeArquivo)), status = 'unknown')
                  write(335,*) lamb0, tf
                  !=====================================================
                  close(335)
               else
                  if(soCalculaIPR)then
                     !==================================================
                     lamb0 = 0.0_dp                    ! Entrou de um jeito
                     call limiarQMF(rede, lamb0, tole)
                     !==================================================
                     nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_Y4_'))//trim(adjustl(teoriaCM))//trim(adjustl('_lamb_index_0'))
                  
                     inquire(unit=444, opened=taAberta)
                     if(taAberta) close(444)    
                  
			         open(444, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
			         write(444,*) lamb0, Y4
   			         close(444)			      
                     !==================================================
			      endif                     
               endif                     
               endif             
            lamb0 = lamb0 + (1.0_dp * ind_lamb) * dlamb
         endif
         !==============================================================
         lambdaf = lamb0 
         lamb = lamb0       
      elseif(nargus == 9)then
         
         nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_lambda'))
 
         inquire( file=trim(adjustl(nomeArquivo))//'.dat', exist=existe ) 

         if(existe)then
            inquire(unit=776, opened=taAberta)
            if(taAberta) close(776)
            !===========================================================
            open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
            read(776,*, iostat=st1) lamb0
            close(776)
            !===========================================================
            if( st1 /= 0)then
               call system('rm -r '//trim(adjustl(local))//'Copia*')
               write(*,*) "Rotina deletou Copia_lambda.dat e se chamou sozinha"
               call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)
            endif  
            !===========================================================          
         else            
            nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//trim(adjustl(teoriaCM))//'_lamb_index_0'
            inquire( file = trim(adjustl(nomeArquivo))//'.dat', exist = existe2)
            
            !===========================================================
            ! Por padrao, essa variavel sera inicalizada assim
            ! A razao esta nas linhas onde lemos lbd_vs_rho_rDMP.dat.
            !===========================================================
            st2 = -1 
            
            if(existe2)then
               !========================================================
               inquire(unit=600, opened=taAberta)
               if(taAberta) close(600)
               !========================================================
               open(600, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
               read(600, *, iostat = st2) lambda_Ultimo_Index, rho1
               close(600)
               if( st2 == 0)then
lams:             do i1 = 1, 10
                     !==================================================
                     write(buffer, '(I0)') i1
                     buffer = trim(adjustl(buffer))
                     nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(buffer))
                     inquire(file=trim(adjustl(nomeArquivo))//'.dat', exist = existe2)
                     !==================================================
                     if(existe2)then
                       !========================================================
                       inquire(unit=600, opened=taAberta)
                       if(taAberta) close(600)
                       !========================================================                     
                        open(600, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
                        read(600, *, iostat = st2) lambda_Ultimo_Index, rho1
                        close(600)
                     else
                        exit lams   
                     endif
                     if(st2 /= 0)then
                        ind_arquivo = i1 - 1
                        write(buffer, '(I0)') ind_arquivo
                        buffer = trim(adjustl(buffer))
                        nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//trim(adjustl(teoriaCM))//'_lamb_index_'//trim(adjustl(buffer))
                        !========================================================
                        inquire(unit=600, opened=taAberta)
                        if(taAberta) close(600)
                        !========================================================                        
                        open(600, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
                        read(600, *, iostat= st2) lambda_Ultimo_Index, rho1
                        close(600)                        
                        exit lams
                     endif                   
                  enddo lams
               endif
            endif      
            !===========================================================                    
            teveLeitura = .False.   
            
            !===========================================================
            inquire(unit=334, opened=taAberta) 
            if( taAberta ) close(334)
            !===========================================================
            nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//trim(adjustl(teoriaCM))
            call leAtehOfimBidimensional(334, trim(adjustl(nomeArquivo))//'.dat', lamb0, rho1, teveLeitura)            
            
            close(334)
            
            if(.not. teveLeitura)then
               if(st2 /= 0)then              
                  !=====================================================
                  lamb0 = 0.0_dp
                  call limiarQMF(rede, lamb0, tole)
                  lamb_C = lamb0
                  !=====================================================
                  write(*,*) "lambda_C = ", lamb_C
                  
                  nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_rho_'))//trim(adjustl(teoriaCM))
                  
                  open(334, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
                  write(334,*) lamb_C, 0.0d0
                  close(334)
                  
                  nomeArquivo = trim(adjustl(local))//trim(adjustl('lbd_vs_tConv_'))//trim(adjustl(teoriaCM))
                  
                  !=====================================================
                  inquire(unit=335, opened=taAberta)            
                  if( taAberta ) close(335)
                  !=====================================================
                  
                  open(335, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
              
                  write(335,*) lamb_C, tf
                  close(335)
               else
                  lamb0 = lambda_Ultimo_Index  
               endif
            else
               if(st2 == 0)then
                  if(lambda_Ultimo_Index > lamb0) lamb0 = lambda_Ultimo_Index
               endif
            endif
            lamb0 = lamb0 + dlamb
            !===============================================================            
         endif       
      endif
!#######################################################################
   	  lamb = lamb0
!#######################################################################      
!=======================================================================
! Se nao existe Copia_lambda.dat, cai fora!
! Se existe, tenta a sorte vendo se o restante existe mesmo.
!=======================================================================
      if( .not. existe )then
         usouCopia = .False.
         rho0 = p0
         return
      endif
!=======================================================================
      
      if( nargus == 10)then
         !=========================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ii'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
                    
         inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )
                           
         if(existe)then
                    
            call system ('tar -xzf '//trim(adjustl(nomeArquivo))//'.tar.gz')

            !===========================================================            
            inquire(unit=777, opened=taAberta)
            if(taAberta)close(777)
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
            !===========================================================
            !===========================================================                        
            do i1 = 1, size(Ii)
               read(777, *, iostat=st) Ii(i1)               
               if(st /= 0)exit
            enddo
            !===========================================================
            close(777)
            !===========================================================
            if(st /= 0)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_*'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
               !==================================================
               call system('rm -r '//trim(adjustl(nomeArquivo))//'.*')
               call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)
               !==================================================            
            endif         
         else
            write(*,*) "Faltou Copia Ii" 
            nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_*'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            call system('rm -r '//trim(adjustl(nomeArquivo))//'*')                 
            call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)                       
         endif
         !===========================================================
         inquire(unit=777, opened=taAberta)
         if(taAberta) close(777)
         !===========================================================                             
         call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
         !===========================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !=========================================================
         nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ri'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
                    
         inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )
                           
         if(existe)then
                    
            call system ('tar -xzf '//trim(adjustl(nomeArquivo))//'.tar.gz')

            !===========================================================            
            inquire(unit=778, opened=taAberta)
            if(taAberta)close(778)
            !===========================================================
            open(778, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
            !===========================================================
            !===========================================================                        
            do i1 = 1, size(Ri)
               read(778, *, iostat=st) Ri(i1)               
               if(st /= 0)exit
            enddo
            !===========================================================
            close(778)
            !===========================================================
            if(st /= 0)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_*'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
               !==================================================
               call system('rm -r '//trim(adjustl(nomeArquivo))//'.*')
               call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)
               !==================================================            
            endif         
         else
            write(*,*) "Faltou Copia Ri" 
            nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_*'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            call system('rm -r '//trim(adjustl(nomeArquivo))//'*')                 
            call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)                       
         endif
         !===========================================================
         inquire(unit=778, opened=taAberta)
         if(taAberta) close(778)
         !===========================================================                             
         call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
         !===========================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@       
         !========================================================           
         nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_t'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
         inquire( file=trim(adjustl(nomeArquivo))//'.dat', exist=existe )

         !==============================================================
         if(existe)then
            !===========================================================
            inquire(unit=782, opened=taAberta)
            if(taAberta)close(782)
            !===========================================================
            open(782, file=trim(adjustl(nomeArquivo))//'.dat', status='old')                 
            read(782, *, iostat=st) t, dt
            close(782)
            if(st /= 0)then
               !========================================================               
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_*'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
               call system('rm -r '//trim(adjustl(nomeArquivo))//'*')
               call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)
               !==================================================            
            endif                                 
         else
            nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_*'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            call system('rm -r '//trim(adjustl(nomeArquivo))//'*')                
            call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)              
         endif                      
         !======================================================
            close(782)
         !======================================================  
      elseif( nargus == 9)then
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ii'))
                    
         inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )
                           
         if(existe)then
                    
            call system ('tar -xzf '//trim(adjustl(nomeArquivo))//'.tar.gz')

            !===========================================================            
            inquire(unit=777, opened=taAberta)
            if(taAberta)close(777)
            !===========================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
            !===========================================================
            !===========================================================                        
            do i1 = 1, size(Ii)
               read(777, *, iostat=st) Ii(i1)               
               if(st /= 0)exit
            enddo
            !===========================================================
            close(777)
            !===========================================================
            if(st /= 0)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_*'))
               !==================================================
               call system('rm -r '//trim(adjustl(nomeArquivo)))
               call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)
               !==================================================            
            endif         
         else
            write(*,*) "Faltou Copia Ii" 
            nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_*'))
            call system('rm -r '//trim(adjustl(nomeArquivo)))                 
            call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)                       
         endif
         !===========================================================
         inquire(unit=777, opened=taAberta)
         if(taAberta) close(777)
         !===========================================================                             
         call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
         !===========================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ri'))
                    
         inquire( file=trim(adjustl(nomeArquivo))//'.tar.gz', exist=existe )
                           
         if(existe)then
                    
            call system ('tar -xzf '//trim(adjustl(nomeArquivo))//'.tar.gz')

            !===========================================================            
            inquire(unit=778, opened=taAberta)
            if(taAberta)close(778)
            !===========================================================
            open(778, file=trim(adjustl(nomeArquivo))//'.dat', status='old')
            !===========================================================
            !===========================================================                        
            do i1 = 1, size(Ri)
               read(778, *, iostat=st) Ri(i1)               
               if(st /= 0)exit
            enddo
            !===========================================================
            close(778)
            !===========================================================
            if(st /= 0)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_*'))
               !==================================================
               call system('rm -r '//trim(adjustl(nomeArquivo)))
               call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)
               !==================================================            
            endif         
         else
            write(*,*) "Faltou Copia Ri" 
            nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_*'))
            call system('rm -r '//trim(adjustl(nomeArquivo)))                 
            call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)                       
         endif
         !===========================================================
         inquire(unit=778, opened=taAberta)
         if(taAberta) close(778)
         !===========================================================                             
         call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
         !===========================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
         nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_t.dat'))
         inquire( file=trim(adjustl(nomeArquivo)), exist=existe )
         !===========================================================
         if(existe)then
            !========================================================
            inquire(unit=782, opened=taAberta)
            if(taAberta)close(782)
            !========================================================
            open(782, file=trim(adjustl(nomeArquivo)), status='old')                 
            read(782, *, iostat = st) t, dt
            close(782)
            !========================================================
            if(st /= 0)then
               !=====================================================
               close(782)
               call system('rm -r '//trim(adjustl(local))//'Copia_*.dat')                  
               call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)
               !==================================================            
            endif                                             
         else             
            call system('rm -r '//trim(adjustl(local))//'Copia*')
            write(*,*) "Faltou Copia tempo"               
            call TestaSeTemCopiaEstadoSistema(teoriaCM, usouCopia)
         endif
         !======================================================
         close(782)
         !======================================================     
      endif   
      !=================================================================
      tempoEscrita = t + tempoEscrita0
      !=================================================================
      rho0 = 0.0d0
      do i1 = 1, size(Ii)   ! Calcula rho
         if(lista_de_clusters(i1) /= i_comp_gigante) cycle
         rho0 = rho0 + Ii(i1)
      enddo
      rho0 = rho0/(1.0d0 * comp_gigante) 
      usouCopia = .True.        
      !=================================================================
   end subroutine       

   subroutine TempoDeEscrever()
         logical :: aberta                
         if( t >= tempoEscrita)then
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_lambda'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_lambda'))
            endif  
            !===========================================================         
            inquire(unit=776, opened=aberta)
            if(aberta)close(776)
            
            open(776, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            write(776,*) lamb
            close(776)
            !===========================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ii'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ii'))
            endif           
                 
            !======================================================
            open(777, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            !======================================================
            do i1 = 1, size(Ii)
               write(777,*) Ii(i1)
            enddo
            !======================================================
            close(777)                 
            !======================================================
            call system('tar -czf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !======================================================
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ri'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_Ri'))
            endif           

            open(778, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
                 
            do i1 = 1, size(Ri)
               write(778,*) Ri(i1)
            enddo
            close(778)
            !===========================================================
            call system('tar -czf '//trim(adjustl(nomeArquivo))//'.tar.gz '//trim(adjustl(nomeArquivo))//'.dat')
            call system('rm '//trim(adjustl(nomeArquivo))//'.dat')
            !===========================================================

            !===========================================================                 
            if(nargus == 10)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_t'))//'_lamb_index_'//trim(adjustl(char_ind_lamb))
            elseif(nargus == 9)then
               nomeArquivo = trim(adjustl(local))//trim(adjustl('Copia_t'))
            endif           
            !======================================================
            open(782, file=trim(adjustl(nomeArquivo))//'.dat', status='unknown')
            write(782,*) t, dt
            close(782)                 
            !======================================================
            tempoEscrita = tempoEscrita + tempoEscrita0
         !=========================================================
         endif
              
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
      
      Ap = 1.0_dp/Ap
      
      g_func = k_s - Ap**(1.0_dp/gam1) * (1.0_dp * Nstr)**(1.0_dp/gam1)
      
   end function                          

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
      kminus = 1.0_dp * kming
      kmais = 1.5_dp * kming * (1.0_dp * N_p)**(1.0_dp/gam_p)
      
      if(allocated(Ap_list)) deallocate(Ap_list)
      allocate(Ap_list(int(kming):(int(kmais))))
      
      gminus = g_func(kminus, gam_p, N_p)
      
      if(gminus >= 0.0_dp) stop "Precisamos de um valor de kminus para que gminus < 0"
       
      Ap1 = 0.0_dp
      do kl1 = kming, int(kmais)
         Ap1 = Ap1 + (1.0_dp * kl1)**(-gam_p)
         Ap_list(kl1) = Ap1
      enddo

      gmais = g_func(kmais, gam_p, N_p)
   
      if(gmais <= 0.0_dp) stop "Precisamos de um valor para kmais, tal que gmais > 0"
      !#################################################################
      !   Execucao
      !#################################################################
      kl1 = 1
      do while(kl1 <= N_it)
         k_p = kminus + (kmais - kminus)/2.0_dp
         gp = g_func(k_p, gam_p, N_p)
         if((gp == 0.0_dp).or.((kmais-kminus)/2.0_dp <= tol))then
            grau_max = k_p
            write(*,*) "Achou o grau max apos N = ", kl1, " iteracoes."
            exit
         endif
         kl1 = kl1 + 1
         if(gminus * gp > 0.0_dp)then
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
   
  
end program
