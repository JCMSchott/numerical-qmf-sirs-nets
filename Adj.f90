program hashimoto
!  use mod_tubos
   use geraRede
   use mod_rndgen
   use types
   
   implicit none

! Calcula o autovalor principal da matriz B_hsh
! 
 
   integer, allocatable :: Adj(:,:), Adj2(:,:)
   integer :: col_hashi
   integer :: i1, i2, i3, j1, j2, j3, j4
   type(grafo_PL_UCM) :: rede
   real(dp) :: gama
   integer :: grau_min, grau_max
   
   type(rndgen) :: geni
   
   integer ::  n_arestas, n_sitios(10), semente, semente2

   real(dp) :: lam

   real(dp), parameter :: p0 = 0.3d0
   real(dp) :: p_sort

   integer :: i4

   real(dp) :: exp_gama(3)
   character(len=30) :: rowfmt
   character(len=500) :: caminho
   character(len=8) :: N_char, gama_char
   integer :: viz_1, viz_2, viz_3, viz_4
   integer :: n_tubos, n_nucleo
   integer :: ipx, ipy
   real(dp) :: xp, yp
   real(dp) :: erro
   real(dp), parameter :: tol = 5.0d-8
   real(dp) :: mu
   integer, parameter :: n_it = 10000
   real(dp), allocatable :: x(:), y(:)
   integer :: n_Adj
   character(len=5) :: n_charac(10)	
   real(dp) :: kmaxi

   caminho = trim(adjustl("/home/jota/SIRS_MP/PLA_Tubos/Rst/"))


	!#########################################################
	exp_gama(1) = 2.30d0; exp_gama(2) = 2.70d0; exp_gama(3) = 3.50d0
	!#########################################################
	! Valores de gamma P \propto k^(-exp_gama)	
	


	n_sitios(1) = 1000; n_sitios(2) = 3000;
	n_sitios(3) = 10000; n_sitios(4) = 30000;
	n_sitios(5) = 100000; n_sitios(6) = 300000;
	n_sitios(7) = 1000000; n_sitios(8) = 3000000;
	n_sitios(9) = 10000000; n_sitios(10) = 30000000;	

	n_charac(1) = trim(adjustl('1k')) ; n_charac(2) = trim(adjustl('3k'))
   n_charac(3) = trim(adjustl('10k')); n_charac(4) = trim(adjustl('30k'))
   n_charac(5) = trim(adjustl('100k')); n_charac(6) = trim(adjustl('300k'))
   n_charac(7) = trim(adjustl('1M')); n_charac(8) = trim(adjustl('3M'))
   n_charac(9) = trim(adjustl('10M')); n_charac(10) = trim(adjustl('30M'))   
!########################################################
! Define parametros

!########################################################

lg:	do j2 = 3, 3 !size(exp_gama)
		
		write(gama_char, '(F5.2)') exp_gama(j2)
		
		open(unit=13, file=trim(adjustl(trim(adjustl(caminho))//'N_vs_lambdaQMF_'//'gama_'//trim(adjustl(gama_char))//'.dat')), status='unknown')

		open(unit=14, file=trim(adjustl(trim(adjustl(caminho))//'N_vs_n_Adj_'//'gama_'//trim(adjustl(gama_char))//'.dat')), status='unknown')

	
ls:		do j3 = 1, size(n_sitios)
			write(N_char, '(I0)') n_sitios(j3)
		   						
			!######################################
			semente = 995887671
			semente2 = semente
			
			call geni%init(semente2)
			
		       n_nucleo = n_sitios(j3)

!                      write(*,*) " "
!                      write(*,*) "Iniciando rotina |PLA_Tubos"
!                      write(*,*) " "
 
             call rede%iniciaGrafo(n_sitios(j3))

   
             if(exp_gama(j2) < 3.0_dp)then
                kmaxi = (1.0_dp * n_sitios(j3))**(0.5_dp)
             else
                kmaxi = (1.0_dp * n_sitios(j3))**(1.0_dp/(exp_gama(j2)-1.0_dp))
             endif

             call rede%inicia(3, kmaxi, exp_gama(j2), semente)

		       call rede%liga(semente, .False.)
				
		       rede%degMin = minval(rede%deg)		
		       rede%degMax = maxval(rede%deg)                      

!		       call PL_Tubos(rede, n_sitios(j3), 3, sqrt(1.0d0 * n_sitios(j3)), exp_gama(j2), semente, 0.25d0, 0.0d0, .True., 2.0d0)

!                       write(*,*) " "
!                       write(*,*) "Conluiu rotina |PLA_Tubos"
!                       write(*,*) " "



!			n_tubos = 0
!			do j1 = 1, rede%nodes
!				if(rede%deg(j1) == 2) n_tubos = n_tubos + 1
!			enddo
				
			rede%degMin = minval(rede%deg)		
			rede%degMax = maxval(rede%deg)
			
			n_arestas = sum(rede%deg)
			
			!########################################################
			n_sitios(j3) = rede%nodes

!			write(*,*) 'O tamanho real da rede e: ', n_sitios(j3)


			!Cria a matriz de Hashimoto
			
			!########################################################
			! Conta o tamanho necessario para a locar B_hsh

			n_Adj = 0
			
!siti:	do i1 = 1, n_sitios(j3)
!vizi  	   do i2 = rede%aux(i1), rede%aux(i1) + rede%deg(i1) -1
!               n_Adj = n_Adj + 1
!			   enddo vizi
!			enddo siti

         n_Adj = sum(rede%deg)
         
			write(14,*) rede%nodes, n_Adj
         !#######################################################

			if(allocated(Adj)) deallocate(Adj)
			
			allocate(Adj(n_Adj,2))
 
			j4 = 0

sit1:			do i1 = 1, n_sitios(j3)
viz1:			   do i2 = rede%aux(i1), rede%aux(i1) + rede%deg(i1) -1
                  viz_1 = rede%listAdj(i2)
	   		      j4 = j4 + 1
	   		      Adj(j4,1) = i1
	   		      Adj(j4,2) = viz_1
	   		            
	   		      if(Adj(j4,1) == 0)then
	   		         write(*,*) 'Elemento nulo em Adj'
					      stop
	   		      endif
	   		            
	   		      if(Adj(j4,2) == 0)then
	   		         write(*,*) 'Elemento nulo em Adj'
	   		         stop
	   		      endif
			   enddo viz1
			enddo sit1
                       
                       !#######################################################
                                                                     
                       !########################################################
			
			if(allocated(x)) deallocate(x)
			allocate(x(n_sitios(j3)))
			
			xp = 0.0d0
			do i2 = 1, size(x)   
			   x(i2) = 2.0d0 * geni%rnd()
			   if(abs(x(i2)) > abs(xp))then
			      xp = x(i2)
			      ipx = i2
			   endif
			enddo
			x = x/xp
			
			if(allocated(y)) deallocate(y)
			allocate(y(n_sitios(j3)))
		        

		   i1 = 0
		        
iter:		do while(i1 < n_it)	
			   !Aqui comeca o algoritmo
			   ipy = 0
			   yp = 0.0d0
			   y = 0.0d0
			   do j4 = 1, n_Adj
			         
			   !Adj(i2,i3)
			   !Adj(j4,1) = i2
	   		!Adj(j4,2) = i3
			         
			   y(Adj(j4,1)) = y(Adj(j4,1)) + x(Adj(j4,2))
			      
			   i2 = Adj(j4, 1)
			   if(abs(y(i2)) > abs(yp))then
			      yp = y(i2)
			      ipy = i2
	         endif
		   enddo
			mu = y(ipy)
			   !#######################
			   ! Se yp = 0, para tudo!
			   if(yp == 0.0d0)then
               do i2 = 1, size(x)   
			         x(i2) = 2.0d0 * geni%rnd()
			         if(abs(x(i2)) > xp)then
			            xp = x(i2)
			            ipx = i2
			         endif
			      enddo
			      x = x/xp
			      i1 = 0
			      write(*,*) 'Autovalor nulo encontrado'
			      cycle iter			   
			   endif
			   
			   
			   erro = abs(maxval(x - y/yp))
			   
			   if(erro < tol)then
			   	lam = 1.0d0/mu
			   	write(13,*) n_sitios(j3), lam
			   	write(*,*) 'Autovalor principal ', mu, ' encontrado com erro relativo ', erro, '. '
			   	write(*,*) " "
			   	write(*,*) 'Foram necessarios ', i1, ' passos.'
			   	write(*,*) " "
			   	cycle ls
			   endif
			   x = y/yp
			   !######################
			   i1 = i1 + 1
			enddo iter
			
			
			
			
			
			                       
			!########################################################
			! Aqui corrigimos, caso a rede tenha numero inferior
			! a n_sitios(j3)
			
			
		
			!write(*,*) "O tamanho do grafo indireto eh: ", n_arestas
			!########################################################
			! Gera o grafo que vamos usar	

			!######################################
			
			!######################################
			! Transforma em character
			
!			open(unit=14, file=trim(adjustl(trim(adjustl(caminho))//'N_vs_Ip_'//trim(adjustl(N_char))//'_gama_'//trim(adjustl(gama_char))//'.dat')), status='unknown')

!			open(unit=15, file=trim(adjustl(trim(adjustl(caminho))//'I_Medio_vs_lambda_N_'//trim(adjustl(N_char))//'_gama_'//trim(adjustl(gama_char))//'.dat')), status='old')

!			open(unit=16, file=trim(adjustl(trim(adjustl(caminho))//'I_MedioNucleo_vs_lambda_N_'//trim(adjustl(N_char))//'_gama_'//trim(adjustl(gama_char))//'.dat')), status='old')
			
		!	write(13,*) n_sitios(j3), lam
			
			deallocate(rede%deg)
			deallocate(rede%aux)
			deallocate(rede%listAdj)
			if(allocated(rede%matriz)) deallocate(rede%matriz)
			
		enddo ls
		close(13)
		close(14)
	enddo lg

!########################################################
	stop
end program
