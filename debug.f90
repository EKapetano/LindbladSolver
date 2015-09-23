!> Contains some routines which might be useful for debugging
module debug
  use accuracy
  use calculations
  implicit none

  contains
    !> Checks the positivity of the density matrix
    !!
    !! \details  This routine checks if a given density matrix fulfills
    !! the positivity condition and writes the result in a file.
    !!
    !! \param density  The system's current density matrix.
    !! \param dim  The system's dimension.
    !! \param tt  Current time.
    !!
    !! \note  The routine writes a 1 into the file when the condition
    !! is violated, otherwise it returns Zero.
    subroutine check_positivity(density,dim,tt)
      complex(dp), intent(in) :: density(:,:)
      real(dp), intent(in) :: tt
      integer, intent(in) :: dim
      !Declarations needed for LA_GEEVX
      complex(dp), allocatable :: tmp3(:,:), eigval(:), VL(:,:), VR(:,:)
      real(dp), allocatable :: SCAL(:), RCONDE(:), RCONDV(:)
      integer :: INFO, ll, ii, ILO, IHI
      real(dp) :: ABNRM
      character(len=1) :: BALANC
      !Settings and allocations
      allocate(tmp3(dim,dim))
      allocate(eigval(dim))
      allocate(VR(dim,dim))
      allocate(VL(dim,dim))
      allocate(SCAL(dim))
      allocate(RCONDE(dim))
      allocate(RCONDV(dim))
      BALANC = 'N'
      tmp3 = density
      !Calculate Eigenvalues with LAPACK95-Routine
      call LA_GEEVX( tmp3, eigval, VL, VR, BALANC, ILO, &
           IHI, SCAL, ABNRM, RCONDE, &
           RCONDV, INFO )
      !Check positivity in the computed Eigenvalues
      !ll = 0 for if positive, else 1
      ll = 0
      do ii = 1, dim
         if(real(eigval(ii)) < 0.0_dp) then
            ll = 1
         end if
      end do
      !Write
      write(20,*) tt, ll
    end subroutine check_positivity

    !> Shows Time-Evolution-Operators and checks if they're unitary
    !!
    !!
    !! \param UU  Time-evolution-operator as a matrix.
    !! \param UUa  Adjoint of UU.
    !! \param dim  The system's dimension.
    !! \param tt  Current time.
    !!
    !! \note  This routine was only meant to be used for dim = 3.
    subroutine write_unitary(UU,UUa,dim,tt)
      complex(dp), intent(in) :: UU(:,:), UUa(:,:)
      real(dp), intent(in) :: tt
      integer, intent(in) :: dim
      complex(dp) :: tmp(3,3)
      integer  :: ii

      tmp = matmul(UU,UUa)
      !Write Operators
      write(*,"(A,F8.2)") "Unitary Operator U at t = ", tt
      do ii = 1, dim
         write(*,"(3(F8.4,F8.4))") UU(ii,:)
      end do
      write(*,*) "Completeness check:"
      do ii = 1, dim
         write(*,"(3(F8.2,F8.2))") tmp(ii,:)
      end do
    end subroutine write_unitary
  end module debug
  
