program qswalk
  use accuracy
  use inout
  use calculations
  use debug
  use f95_lapack
  use omp_lib
  implicit none
  
  !Declarations
  integer :: ii, jj, kk, dim, steps
  complex(dp), allocatable :: tmp(:,:), tmp2(:,:), tmp3(:,:), tmp4(:,:)
  complex(dp), allocatable :: density(:,:), comm(:,:), tmp5(:,:), sigma(:,:)
  complex(dp), allocatable :: lindL(:,:), lindR(:,:), lindD(:,:)
  real(dp), allocatable :: hamilton(:,:), ecurrIN(:), ecurrOUT(:)
  real(dp), allocatable :: pp(:,:), ecurrDEPH(:), entropy(:)
  real(dp) :: dt, tt, dephstr, extin, extout, aa
  !Declarations needed for Expokit-Routine ZGPADM
  integer :: ideg, m, ldh, lwsp, iexph, ns, iflag
  complex(dp), allocatable :: wsp(:), ipiv(:), UU(:,:), UUa(:,:)
  
  !Read Hamiltonian + other Input
  call read_hamilton(steps,dim,hamilton,density,dt,&
       &dephstr,extin,extout)
  !Expand density matrix by 2 dimensions which act as baths
  call calc_expand(dim,density,hamilton)

  !Settings for Expokit-Routine ZGPADM
  allocate(UU(dim,dim))
  allocate(UUa(dim,dim))
  allocate(tmp(dim,dim))
  ideg = 8
  m = dim
  ldh = dim
  lwsp = 4 * m**2 + ideg + 2
  allocate(wsp(lwsp))
  allocate(ipiv(m))

  !Initial probabilities
  allocate(pp(dim,steps))
  !$OMP PARALLEL DO
  do ii = 1, dim
     pp(ii,1) = real(density(ii,ii))
  end do
  !$OMP END PARALLEL DO
  
  !Calculate Unitary matrix U = exp(- i*H*dt)
  tmp = cmplx(0.0_dp, - 1.0_dp, dp) * hamilton
  call ZGPADM(ideg,m,dt,tmp,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)
  do ii = 1, dim
     !$OMP PARALLEL DO
     do jj = 1, dim
        UU(ii,jj) = wsp(iexph + (ii-1) * m + (jj-1))
     end do
     !$OMP END PARALLEL DO
  end do
  !Adjoint of U
  UUa = conjg(transpose(UU))
  deallocate(tmp)
  

  !Open Files
  open(11,file="realstats.dat",status="replace",&
       &action="write",form="formatted")
  open(12,file="imagstats.dat",status="replace",&
       &action="write",form="formatted")
  !File for positivity check
  open(20,file="positive.dat",status="replace",&
       &action="write",form="formatted")
  !File for the relative entropy
  open(21,file="entropy.dat",status="replace",&
       &action="write",form="formatted")
  !Calculation
  !=============================================!
  !Additional allocations
  allocate(comm(dim,dim))
  allocate(tmp(dim,dim))
  allocate(tmp2(dim,dim))
  allocate(ecurrIN(steps))
  allocate(ecurrOUT(steps))
  allocate(ecurrDEPH(steps))
  allocate(tmp3(dim,dim))
  allocate(tmp4(dim,dim))
  allocate(tmp5(dim,dim))
  allocate(lindD(dim,dim))
  allocate(lindL(dim,dim))
  allocate(lindR(dim,dim))
  allocate(entropy(steps))
  allocate(sigma(dim,dim))

  !$OMP PARALLEL NUM_THREADS(3)
  do kk = 1, steps
     !$OMP SINGLE
     tt = kk * dt
     !Refresh baths
     density(1,1) = cmplx( 0.5_dp, 0.0_dp, dp )
     density(dim,dim) = 0.0_dp
     !$OMP END SINGLE
     
 
     !$OMP SECTIONS
     !Calculate Lindblad-Operators for In-/Output
     !Input
     !$OMP SECTION
     call calc_leftlindblad(density,tmp,tmp2,lindL)
     call calc_energycurr(hamilton,tmp,lindL,ecurrIN,kk,dim,extin)
     !Output
     !$OMP SECTION
     call calc_rightlindblad(density,tmp3,tmp4,lindR,dim)
     call calc_energycurr(hamilton,tmp3,lindR,ecurrOUT,kk,dim,extout)
     !$OMP SECTION
     call calc_dephlindblad(density,tmp5,dim,lindD)
     call calc_energycurr(hamilton,tmp5,lindD,ecurrDEPH,kk,dim,dephstr)
     !$OMP END SECTIONS

     !Evolution     
     !$OMP SINGLE
     density = density + dt * ( extin * lindL + extout * lindR&
          & + dephstr * lindD)
     density = matmul(UU,matmul(density,UUa))
     !$OMP END SINGLE
     
     !===========================================!
     !Write new density matrix into file
     !call write_density(density,t)
     
     !$OMP SECTIONS
     !$OMP SECTION
     !write(11,*) tt, real(density(2:dim-1,2:dim-1))
     !$OMP SECTION
     !write(12,*) tt, aimag(density(2:dim-1,2:dim-1))
     !$OMP END SECTIONS NOWAIT
     
     !New probabilities
     !$OMP DO
     do ii = 1, dim
        pp(ii,kk) = real(density(ii,ii))
     end do
     !$OMP END DO

     !Calculate Entropy
     !$OMP SINGLE
     call calc_rel_entropy(entropy,kk,dim,density,sigma,tmp)
     !$OMP END SINGLE NOWAIT
     !Check positivity of the density Matrix for debugging purposes(optional)
     !call check_positivity(density,dim,t)
     !Stop condition(optional)
     !if(pp(dim,kk) > 1.0_dp/(2.0_dp * dim)) then
     !   write(*,"(A,F8.3)") "Probability reached at t = ", t
     !   stop
     !end if
  end do
  !$OMP END PARALLEL
  !================================================!
  !Close files
  close(11,status="keep")
  close(12,status="keep")
  close(20,status="keep")
  close(21,status="keep")
  !Write probabilities into files
  call write_prob(pp,steps,dim,dt)
  !Write final density matrix into file
  call write_lastdens(density,dim)
  !Write relative entropy into file
  call write_entropy(dt,steps,entropy)
  !Write Energy currents into file
  call write_energycurr(dt,ecurrIN,ecurrOUT,ecurrDEPH,steps)
  !(Optional) Calculate currents
  call calc_curr(dim,steps,hamilton,pp,dt)
  !(Optional) Write Voltage on screen
  aa = real(density(2,2)) - real(density(dim-1,dim-1))
  write(*,"(A,F30.15)") "Voltage: ", aa

end program qswalk
