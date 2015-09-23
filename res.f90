!Calculates effective resistance depending on a Parameter PP
program resistance
  use accuracy
  use inout
  use calculations
  use f95_lapack
  use omp_lib
  implicit none

  !Declarations
  complex(dp), allocatable :: inidens(:,:), density(:,:), tmp(:,:)
  real(dp) :: dt, dephstr, extin, extout, voltage, resist, bb, sample
  real(dp) :: eIN, eOUT, eDEPH
  real(dp), allocatable :: hamilton(:,:), array(:,:)
  integer :: ii, jj, num, steps, dim
  character(len=1) :: AA
  !Declarations for Expokit-Routine ZGPADM
  integer :: ideg, m, ldh, lwsp, iexph, ns, iflag
  complex(dp), allocatable :: wsp(:), ipiv(:), UU(:,:), UUa(:,:)

  !Read Hamiltonian + other Input
  call read_hamilton(steps,dim,hamilton,inidens,dt,&
       &dephstr,extin,extout)
  !Expand density matrix by 2 dimensions which act as baths
  call calc_expand(dim,inidens,hamilton)
  
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
  !Calculate Unitary matrix U = exp(- i*H*dt)
  tmp = cmplx(0.0_dp, - 1.0_dp, dp) * hamilton
  call ZGPADM(ideg,m,dt,tmp,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)
  do ii = 1, dim
     do jj = 1, dim
        UU(ii,jj) = wsp(iexph + (ii-1) * m + (jj-1))
     end do
  end do
  !Adjoint of U
  UUa = conjg(transpose(UU))
  deallocate(tmp)

  allocate(density(dim,dim))
  !Read number of iterations and desired sampling
  write(*,"(A)",advance="no") "Number of Iterations:"
  read(*,*) num
  allocate(array(num,5))
  write(*,"(A)",advance="no") "Simple sampling or logarithmic?(S/L)"
  read(*,*) AA
  write(*,"(A)",advance="no") "Sample from 0 to:"
  read(*,*) sample
  
  !Quantum Walks with changing parameter PP
  if(AA == 'S') then
     !Simple sampling
     !$OMP PARALLEL DO PRIVATE(density,dephstr,voltage,resist,eIN,eOUT,eDEPH)
     do ii = 1, num
        write(*,*) ii
        dephstr = (ii-1) * sample / num
        density = inidens
        call calc_qswalk(steps,dim,UU,UUa,density,dt,&
             &dephstr,extin,extout,voltage,eIN,eOUT,eDEPH,hamilton)
        !Write effective resistance
        resist = voltage
        array(ii,1) = dephstr
        array(ii,2) = resist
        array(ii,3) = eIN
        array(ii,4) = eOUT
        array(ii,5) = eDEPH
     end do
     !$OMP END PARALLEL DO
  else if(AA == 'L') then
     !Logarithmic sampling
     do ii = num, 1, -1
        bb = - ii * 8.0_dp / num
        dephstr = exp(bb)
        density = inidens
        call calc_qswalk(steps,dim,UU,UUa,density,dt,&
             &dephstr,extin,extout,voltage,eIN,eOUT,eDEPH,hamilton)
        !Effective resistance
        resist = voltage
        jj = num + 1 - ii
        array(jj,1) = dephstr
        array(jj,2) = resist
        array(jj,3) = eIN
        array(jj,4) = eOUT
        array(jj,5) = eDEPH
     end do
  else
     write(*,*) "Wrong specifier!"
     stop
  end if
  !Write to file
  open(11,file="resi.dat",form="formatted",&
       &action="write",status="replace")
  do ii = 1, num
     write(11,*) array(ii,:)
  end do
  close(11,status="keep")

!End
contains


  !Contains a simplified version of the normal program with the only
  !purpose of calculating the steady-state voltage
  subroutine calc_qswalk(steps,dim,UU,UUa,density,dt,&
       &dephstr,extin,extout,voltage,eIN,eOUT,eDEPH,hamilton)
    !Declarations
    complex(dp), intent(in) :: UU(:,:), UUa(:,:)
    complex(dp), intent(inout) :: density(:,:)
    real(dp), intent(in) :: dephstr, extin, extout, dt, hamilton(:,:)
    real(dp), intent(out) :: voltage, eIN, eOUT, eDEPH
    integer, intent(in) :: steps, dim
    integer :: kk
    complex(dp), allocatable :: tmp(:,:), tmp2(:,:), tmp3(:,:)
    complex(dp), allocatable :: tmp4(:,:), tmp5(:,:)
    complex(dp), allocatable :: lindL(:,:), lindR(:,:), lindD(:,:)
    real(dp) :: tt
    !Additional allocations
    allocate(tmp(dim,dim))
    allocate(tmp2(dim,dim))
    allocate(tmp3(dim,dim))
    allocate(tmp4(dim,dim))
    allocate(tmp5(dim,dim))
    allocate(lindD(dim,dim))
    allocate(lindL(dim,dim))
    allocate(lindR(dim,dim))
    !Calculation
    do kk = 1, steps
       !==========================================!
       tt = kk * dt
       !Refresh baths
       density(1,1) = cmplx( 1.0_dp, 0.0_dp, dp )
       density(dim,dim) = 0.0_dp
       !Calculate Lindblad-Operators for In-/Output and Dephasing
       call calc_leftlindblad(density,tmp,tmp2,lindL)
       call calc_rightlindblad(density,tmp3,tmp4,lindR,dim)
       call calc_dephlindblad(density,tmp5,dim,lindD)
       !Evolution     
       density = density + dt * ( extin * lindL + extout * lindR&
            & + dephstr * lindD)
       density = matmul(UU,matmul(density,UUa))   
    end do
    !================================================!    
    !Output
    voltage = real(density(2,2)) - real(density(dim-1,dim-1))
    !Steady-State Energy Currents
    call calc_ecurr(hamilton,lindL,eIN,dim,extin)
    call calc_ecurr(hamilton,lindR,eOUT,dim,extout)
    call calc_ecurr(hamilton,lindD,eDEPH,dim,dephstr)
  end subroutine calc_qswalk



  !Calculates energy currents caused by an operator L(rho)
  subroutine calc_ecurr(hamilton,lind,ecurr,dim,str)
    complex(dp), intent(in) :: lind(:,:)
    real(dp), intent(in) :: hamilton(:,:), str
    real(dp), intent(out) :: ecurr
    integer, intent(in) :: dim
    complex(dp), allocatable :: tmp(:,:)
    integer :: ii
    
    allocate(tmp(dim,dim))
    tmp = matmul(hamilton, lind)
    ecurr = 0.0_dp
    !Calculate the Trace
    do ii = 1, dim
       ecurr = ecurr + real(tmp(ii,ii))
    end do
    ecurr = str * ecurr
    deallocate(tmp)
  end subroutine calc_ecurr
end program resistance
