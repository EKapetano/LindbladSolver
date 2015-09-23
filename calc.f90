!> Module containing most of the calculations
!!
!! This module contains the routines for most of the calculations, f.e.
!! the calculation of the Lindblad-Operators and the Relative Entropy.
module calculations
  use accuracy
  use f95_lapack
  use omp_lib
  implicit none
  !Definition for Pi
  real(dp), parameter :: PI = 4.0_dp * atan(1.0_dp)
  !Contains some of the Calculations
contains


  !> Calculates the commutator - i * [H,rho]
  !!
  !! \details This routine calculates the commutator -i * [H,rho] from
  !! a given Hamiltonian H and density matrix Rho.
  !!
  !!
  !! \param hamilton  Contains the system's Hamiltonian matrix
  !! \param density  Contains the density matrix
  !! \param comm  Contains the commutator -i * [H,rho] on exit
  subroutine calc_commutator(comm, hamilton, density)
    complex(dp), intent(out) :: comm(:,:)
    complex(dp), intent(in) :: density(:,:)
    real(dp), intent(in) :: hamilton(:,:)

    comm = 0.0_dp
    comm = cmplx(0.0_dp, - 1.0_dp, dp) * &
         &( matmul(hamilton,density)&
         & - matmul(density,hamilton))
  end subroutine calc_commutator

  !> Expands the Hamiltonian and the density matrix by 2 levels
  !!
  !! \details This routine expands the density matrix and the system's
  !! Hamiltonian by 2 levels which will later act as the particle baths.
  !!
  !!
  !! \param dim  Initial dimension of the system. Is raised by 2 on exit.
  !! \param density  The system's initial density matrix. Contains 2 empty
  !! additional rows and columns on exit.
  !! \param hamilton  The system's Hamilton matrix. Contains 2 empty
  !! additional rows and columns on exit.
  subroutine calc_expand(dim,density,hamilton)
    complex(dp), allocatable, intent(inout) :: density(:,:)
    real(dp), allocatable, intent(inout) :: hamilton(:,:)
    complex(dp), allocatable :: tmp(:,:)
    real(dp), allocatable :: tmp2(:,:)
    integer, intent(inout) :: dim

    !$OMP PARALLEL
    !$OMP SECTIONS
    !$OMP SECTION
    allocate(tmp(dim,dim))
    tmp = density
    deallocate(density)
    !$OMP SECTION
    allocate(tmp2(dim,dim))
    tmp2 = hamilton
    deallocate(hamilton)
    !$OMP END SECTIONS
    !$OMP SINGLE
    dim = dim + 2
    !$OMP END SINGLE
    !$OMP SECTIONS
    !$OMP SECTION
    allocate(density(dim,dim))
    density = 0.0_dp
    density(2:dim-1,2:dim-1) = tmp
    deallocate(tmp)
    !$OMP SECTION
    allocate(hamilton(dim,dim))
    hamilton = 0.0_dp
    hamilton(1,1) = 1.0_dp
    hamilton(dim,dim) = 1.0_dp
    hamilton(2:dim-1,2:dim-1) = tmp2
    deallocate(tmp2)
    !$OMP END SECTIONS
    !$OMP END PARALLEL
  end subroutine calc_expand

  !> Calculates the Lindblad-Operator for Dephasing
  !!
  !!
  !! \details This routine calculates the Lindblad-Operator which is used
  !! to model dephasing.
  !!
  !! \param density  The system's current density matrix.
  !! \param tmp  Work-Array.
  !! \param dim  Dimension of the system.
  !! \param lindD  Contains the Dephasing-Operator on exit.
  !!
  !! \note The Work-Array should have the same dimension as the
  !! density matrix.
  subroutine calc_dephlindblad(density,tmp,dim,lindD)
    complex(dp), intent(in) :: density(:,:)
    complex(dp), intent(out) :: tmp(:,:), lindD(:,:)
    integer, intent(in) :: dim
    integer :: ii
    
    lindD = 0.0_dp
    do ii = 1, dim
       tmp = 0.0_dp
       tmp(ii,ii) = cmplx(1.0_dp, 0.0_dp, dp)
       lindD = lindD + 2.0_dp * density(ii,ii) * tmp - &
            &matmul(tmp,density) - matmul(density,tmp)
    end do    
  end subroutine calc_dephlindblad


  !> Calculates the Lindblad-Operator for adding particles
  !!
  !!
  !! \details This routine calculates the Lindblad-Operator which swaps
  !! particles from the particle source to the system.
  !!
  !! \param density  The system's current density matrix.
  !! \param tmp  Work-Array.
  !! \param tmp2  Work-Array.
  !! \param lindL  Contains the desired Operator on exit.
  !!
  !! \note The Work-Arrays should have the same dimension as the
  !! density matrix.
  subroutine calc_leftlindblad(density,tmp,tmp2,lindL)
    complex(dp), intent(inout) :: density(:,:)
    complex(dp), intent(out) :: tmp(:,:), tmp2(:,:), lindL(:,:)

    tmp = 0.0_dp
    tmp(1,1) = cmplx(1.0_dp, 0.0_dp, dp)
    tmp2 = 0.0_dp
    tmp2(2,2) = cmplx(1.0_dp, 0.0_dp, dp)
    lindL = ( 2.0_dp * density(1,1) * tmp2 - matmul(tmp,density)&
         & - matmul(density,tmp) )
  end subroutine calc_leftlindblad

  !> Calculates the Lindblad-Operator for deleting particles
  !!
  !!
  !! \details This routine calculates the Lindblad-Operator which swaps
  !! the particles out of the system once they reach the end of the circuit.
  !!
  !! \param density  The system's current density matrix.
  !! \param tmp  Work-Array.
  !! \param tmp2  Work-Array.
  !! \param dim  Dimension of the system.
  !! \param lindR  Contains the desired Operator on exit.
  !!
  !! \note The Work-Arrays should have the same dimension as the
  !! density matrix.
  subroutine calc_rightlindblad(density,tmp,tmp2,lindR,dim)
    complex(dp), intent(inout) :: density(:,:)
    complex(dp), intent(out) :: tmp(:,:), tmp2(:,:), lindR(:,:)
    integer, intent(in) :: dim

    tmp = 0.0_dp
    tmp(dim-1,dim-1) = cmplx(1.0_dp, 0.0_dp, dp)
    tmp2 = 0.0_dp
    tmp2(dim,dim) = cmplx(1.0_dp, 0.0_dp, dp)
    lindR =( 2.0_dp * density(dim-1,dim-1) * tmp2 - matmul(tmp,density)&
         & - matmul(density,tmp) )
  end subroutine calc_rightlindblad
  

  !> Calculates Energy-Currents caused by a Lindblad-Operator
  !!
  !!
  !! \details This routine calculates the Energy-Current in the system
  !! caused by a given Lindblad-Operator.
  !!
  !! \param hamilton  Contains the system's Hamilton matrix.
  !! \param tmp  Work-Array.
  !! \param lind  Contains the Lindblad-Operator which causes the
  !! Energy-Current.
  !! \param kk  Contains the current timestep.
  !! \param dim  Dimension of the system.
  !! \param cc  Contains the strength of the Lindblad-Operator.
  !! \param ecurr Contains the Energy-Current caused by the Lindblad-Operator
  !! at the timestep kk
  !!
  !! \note The Work-Array should have the same dimension as the
  !! density matrix.
  subroutine calc_energycurr(hamilton,tmp,lind,ecurr,kk,dim,cc)
    complex(dp), intent(in) :: lind(:,:)
    complex(dp), intent(out) :: tmp(:,:)
    real(dp), intent(in) :: hamilton(:,:), cc
    real(dp), intent(out) :: ecurr(:)
    integer, intent(in) :: kk, dim
    integer :: ii
    
    tmp = matmul(hamilton, lind)
    ecurr(kk) = 0.0_dp
    !Calculate the Trace
    do ii = 1, dim
       ecurr(kk) = ecurr(kk) + real(tmp(ii,ii))
    end do
    ecurr(kk) = cc * ecurr(kk)
  end subroutine calc_energycurr


  !> Calculates and writes Particle-Currents
  !!
  !!
  !! \details This routine calculates the particle currents from the
  !! systems Hamiltonian and the previously calculated populations.
  !! After the calculation, the currents are written into a file.
  !!
  !! \param dim  Dimension of the system.
  !! \param steps  Total number of timesteps of the simulation.
  !! \param hamilton  The system's Hamilton matrix.
  !! \param dt  Length of a timestep.
  !! \param pp  Contains the previously calculated populations.
  subroutine calc_curr(dim,steps,hamilton,pp,dt)
    integer, intent(in) :: dim, steps
    real(dp), intent(in) :: hamilton(:,:), pp(:,:), dt
    real(dp), allocatable :: current(:,:), tmp(:,:), pp2(:,:)
    integer :: ii, kk, aa
    !Calculate reduced hamiltonian/probabilities
    allocate(tmp(dim-2,dim-2))
    allocate(pp2(dim-2,steps))
    tmp(1:dim-2,1:dim-2) = hamilton(2:dim-1,2:dim-1)
    pp2(1:dim-2,1:steps) = pp(2:dim-1,1:steps)
    aa = dim - 2
    !Calculate Currents
    allocate(current(aa-1, steps))
    do kk = 1, steps
       !$OMP PARALLEL
       !$OMP DO
       do ii = 1, aa-1
          current(ii,kk) = (- tmp(ii+1,ii) * pp2(ii,kk)&
               & + tmp(ii,ii+1) * pp2(ii+1,kk))
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    end do
    current = current / (pp2(1,steps) - pp2(aa,steps))
    !Write to File
    open(11,file="current.dat",form="formatted",&
         &action="write",status="replace")
    do kk = 1, steps
       write(11,*) kk * dt, current(:,kk)
    end do
    close(11,status="keep")
  end subroutine calc_curr
  
  
  !> Calculates the relative entropy between a density matrix and it's
  !> decohered counterpart
  !!
  !!
  !! \details This subroutine calculates the decohered density matrix
  !! in position basis. Afterwards, the relative entropy between the
  !! density matrix and it's decohered counterpart is being determined.
  !! In order to obtain the eigenvalues of the density matrix, a
  !! LAPACK95-Routine is used.
  !!
  !! \param kk  Current timestep of the simulation.
  !! \param dim  Dimension of the system.
  !! \param density  Current density matrix of the system.
  !! \param sigma  Work-Array. Contains the decohered density matrix on exit.
  !! \param tmp  Work-Array.
  !!
  !! \note The Work-Arrays should have the same dimension as the
  !! density matrix.
  subroutine calc_rel_entropy(entropy,kk,dim,density,sigma,tmp)
    complex(dp), intent(in) :: density(:,:)
    complex(dp), intent(out) :: sigma(:,:), tmp(:,:)
    real(dp), intent(out) :: entropy(:)
    complex(dp), allocatable :: tmp2(:,:)
    integer, intent(in) :: kk, dim
    real(dp) :: aa, bb
    integer :: ii, jj
    !Declarations for LAPACK95-Routine LA_HEEV
    real(dp), allocatable :: eigval(:)
    allocate(eigval(dim))
    
    !First Term: Tr(rho * ln(rho))
    tmp = density
    !Eigenvalues
    call LA_HEEV(tmp,eigval)
    aa = 0.0_dp
    do ii = 2, dim-1
       aa = aa + eigval(ii) * log(eigval(ii))
    end do
    !Second Term: Tr(rho * ln(sigma))
    sigma = density
    !Calculate decohered matrix
    do ii = 1, dim
       do jj = 1, dim
          if(ii == jj) then
             !Don't change diagonals
          else
             !Set offdiagonals to 0
             sigma(ii,jj) = 0.0_dp
          end if
       end do
    end do
    do ii = 1, dim
       sigma(ii,ii) = cmplx(real(sigma(ii,ii)) + 0.001_dp,0.0_dp,dp)
       sigma(ii,ii) = log(sigma(ii,ii))
    end do
    allocate(tmp2(dim-2,dim-2))
    tmp2 = matmul(density(2:dim-1,2:dim-1),sigma(2:dim-1,2:dim-1))

    !Trace
    bb = 0.0_dp
    do ii = 1, dim-2
       bb = bb + abs(tmp2(ii,ii))
    end do

    !Relative Entropy
    entropy(kk) = abs(aa - bb)
  end subroutine calc_rel_entropy


end module calculations
