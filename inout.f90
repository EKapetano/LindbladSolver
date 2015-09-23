!> Module containing input/output-related routines
module inout
  use accuracy
  use omp_lib
  implicit none
  
contains

  !> Reads Hamiltonian, initial state and external input
  !! 
  !!
  !! \details This routine reads the system's Hamiltonian, dimension,
  !! initial state, in-/output rates, dephasing strength, time discretization
  !! and the desired number of timesteps from the hamilton.inp file.
  !!
  !! \param steps  Number of timesteps.
  !! \param dim  Dimension of the system.
  !! \param hamilton  The system's Hamiltonian matrix.
  !! \param density  Initial density matrix.
  !! \param dt  Desired time-discretization (length of a timestep).
  !! \param dephstr  Desired Dephasing strength.
  !! \param extin  Input rate for the external particle current.
  !! \param extout  Output rate for the external particle current.
  !!
  !! \note The hamilton.inp.example - file contains an example showing how
  !! an input-file looks.
  subroutine read_hamilton(steps,dim,hamilton,density,dt,&
       &dephstr,extin,extout)
    integer, intent(out) :: steps, dim
    real(dp), allocatable, intent(out) :: hamilton(:,:)
    complex(dp), allocatable, intent(out) :: density(:,:)
    real(dp), intent(out) :: dt, dephstr, extin, extout
    real(dp), allocatable :: tmp(:,:), tmp2(:,:)
    character(len=1) :: old
    integer :: ii, jj

    open(10, file="hamilton.inp", status="old",&
         &form="formatted", action="read")
    read(10,*) steps
    read(10,*) dim
    !Allocate arrays
    allocate(hamilton(dim,dim))
    allocate(density(dim,dim))
    !Hamiltonian
    do ii = 1, dim
       read(10,*) hamilton(ii,:)
    end do
    read(10,*) dt
    read(10,*) dephstr
    read(10,*) extin, extout
    read(10,*) old
    !Read initial density matrix
    if(old == 'Y') then
       allocate(tmp(dim,dim))
       allocate(tmp2(dim,dim))
       !Read final matrix from last run
       open(11,file="lastdens.dat",status="old",&
            &form="formatted",action="read")
       do ii = 1, dim
          read(11,*) tmp(ii,:)
       end do
       do ii = 1, dim
          read(11,*) tmp2(ii,:)
       end do
       close(11,status="keep")
       !Density matrix
       do ii = 1, dim
          do jj = 1, dim
             density(ii,jj) = cmplx(tmp(ii,jj), tmp2(ii,jj),dp)
          end do
       end do
    else
       !Read density matrix directly from input-file
       do ii = 1, dim
          read(10,*) density(ii,:)
       end do
    end if
    close(10, status="keep")
  end subroutine read_hamilton

  !> Writes Energy-Currents into file
  !!
  !!
  !! \details This routine writes the Energy-Currents, which are caused
  !! by external particle currents and dephasing, into a file.
  !!
  !! \param dt  Time-discretization which has been used for the run.
  !! \param ecurrIN  Contains the Energy-Current caused by adding particles.
  !! \param ecurrOUT  Contains the Energy-Current caused by erasing particles.
  !! \param ecurrDEPH  Contains the Energy-Current caused by dephasing.
  !! \param steps  Number of timesteps which has been used for the run.
  subroutine write_energycurr(dt,ecurrIN,ecurrOUT,ecurrDEPH,steps)
    real(dp), intent(in) :: dt, ecurrIN(:), ecurrOUT(:), ecurrDEPH(:)
    integer, intent(in) :: steps
    integer :: ii
    

    open(11,file="ecurr.dat",action="write",&
         &form="formatted",status="replace")
    do ii = 1, steps
       write(11,*) ii * dt, ecurrIN(ii), ecurrOUT(ii), ecurrDEPH(ii)
    end do
    close(11,status="keep")
  end subroutine write_energycurr

  !> Writes the relative entropy into a file
  !! 
  !! 
  !! \details This routine writes the relative entropy between the 
  !! density matrix and it's coherence-free counterpart into a file.
  !!
  !! \param dt Time-discretization which has been used for the run.
  !! \param steps  Number of timesteps which has been used for the run.
  !! \param entropy  Contains the previously calculated relative entropies.
  subroutine write_entropy(dt,steps,entropy)
    real(dp), intent(in) :: dt, entropy(:)
    integer, intent(in) :: steps
    integer :: ii

    open(11,file="entropy.dat",action="write",&
         &form="formatted",status="replace")
    do ii = 1, steps
       write(11,*) ii * dt, entropy(ii)
    end do
    close(11,status="keep")
  end subroutine write_entropy


  !> Writes last density matrix into a file for later use
  !!
  !!
  !! \details This routine writes the density matrix at the end of the run
  !! into a file. By specifying a variable in the hamilton.inp - file, the
  !! next run of the code will start from that matrix.
  !!
  !! \param density  Density-matrix at the end of the run.
  !! \param dim  Dimension of the system.
  subroutine write_lastdens(density,dim)
    complex(dp), intent(in) :: density(:,:)
    integer, intent(in) :: dim
    integer :: ii
    open(14,file="lastdens.dat",status="replace",&
         &form="formatted",action="write")
    !Real part
    do ii = 2, dim-1
       write(14,*) real(density(ii,2:dim-1))
    end do
    !Imaginary part
    do ii = 2, dim-1
       write(14,*) aimag(density(ii,2:dim-1))
    end do
    close(14,status="keep")

  end subroutine write_lastdens
 

  !> Write probabilities/populations into files
  !!
  !!
  !! \details This routine writes the particle populations in the system
  !! (the diagonal elements of the density matrix) into a file as columns.
  !!
  !! \param pp  Contains the previously calculated populations as rows.
  !! \param steps  Number of timesteps used in the run.
  !! \param dim  The system's dimension.
  !! \param dt  Time-discretization which has been used in the run.
  subroutine write_prob(pp,steps,dim,dt)
    real(dp), intent(in) :: pp(:,:), dt
    real(dp), allocatable :: pp2(:,:)
    integer, intent(in) :: steps, dim
    integer :: ii

    !Write as columns
    allocate(pp2(steps+1,dim))
    pp2 = transpose(pp)
    open(12, file="histo.dat", status="replace",&
         &form="formatted", action="write")
    do ii = 1, steps
       write(12,*) (ii-1) * dt, pp2(ii,2:dim-1)
    end do
    close(12, status="keep")
  end subroutine write_prob



end module inout
