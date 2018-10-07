module diag_mod

  use ieee_arithmetic
  use parallel_mod
  use mesh_mod
  use data_mod
  use types_mod
  use log_mod
  use params_mod

  implicit none

  private

  public diag_init
  public diag_run
  public diag_final
  public diag_total_energy
  public diag_type
  public diag

  type diag_type
    real total_mass
    real total_energy
    real, allocatable :: vor(:,:)
    real, allocatable :: div(:,:)
  end type diag_type

  type(diag_type) diag

contains

  subroutine diag_init()

    if (.not. allocated(diag%vor)) call parallel_allocate(diag%vor)
    if (.not. allocated(diag%div)) call parallel_allocate(diag%div)

    call log_notice('Diag module is initialized.')

  end subroutine diag_init

  subroutine diag_run(state)

    type(state_type), intent(in) :: state

    real vm1, vp1, um1, up1
    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vm1 = state%v(i-1,j)
        vp1 = state%v(i+1,j)
        um1 = state%u(i,j-1) * mesh%full_cos_lat(j-1)
        up1 = state%u(i,j+1) * mesh%full_cos_lat(j+1)
        diag%vor(i,j) = 0.5 * ((vp1 - vm1) / coef%full_dlon(j) - (up1 - um1) / coef%full_dlat(j))
        um1 = state%u(i-1,j)
        up1 = state%u(i+1,j)
        vm1 = state%v(i,j-1) * mesh%full_cos_lat(j-1)
        vp1 = state%v(i,j+1) * mesh%full_cos_lat(j+1)
        diag%div(i,j) = 0.5 * ((up1 - um1) / coef%full_dlon(j) - (vp1 - vm1) / coef%full_dlat(j))
      end do
    end do

    diag%total_mass = 0.0
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        diag%total_mass = diag%total_mass + mesh%full_cos_lat(j) * mesh%dlon * mesh%dlat * state%gd(i,j)
      end do
    end do
    diag%total_mass = diag%total_mass * radius**2

    if (ieee_is_nan(diag%total_mass)) then
      call log_error('Total mass is NaN!')
    end if

    diag%total_energy = diag_total_energy(state)

    if (ieee_is_nan(diag%total_energy)) then
      call log_error('Total energy is NaN!')
    end if

  end subroutine diag_run

  subroutine diag_final()

    if (allocated(diag%vor)) deallocate(diag%vor)
    if (allocated(diag%div)) deallocate(diag%div)

  end subroutine diag_final

  real function diag_total_energy(state) result(res)

    type(state_type), intent(in) :: state

    integer i, j

    res = 0.0
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + (state%iap%u(i,j)**2 + state%iap%v(i,j)**2 + (state%gd(i,j) + static%ghs(i,j))**2) * mesh%full_cos_lat(j)
      end do
    end do

  end function diag_total_energy

end module diag_mod
