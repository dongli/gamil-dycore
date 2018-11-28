module weno_mod

  use params_mod
  use log_mod
  use parallel_mod
  use types_mod
  use data_mod

  private

  public weno_init
  public weno_zonal_momentum_advection_operator
  public weno_meridional_momentum_advection_operator
  public weno_final

  real, parameter :: eps = 1.0d-6
  real, parameter :: u_max = 20.0d0
  real, parameter :: v_max = 20.0d0

  ! Interpolation coefficients
  real, allocatable :: c(:,:)
  ! Smooth indicator coefficients
  real, allocatable :: cb(:)
  ! Optimal stencil weights
  real, allocatable :: wo(:)
  
  ! Positive fluxes at cell interfaces for each stencil
  real, allocatable :: fp_u(:,:), fp_v(:,:)
  ! Negative fluxes at cell interfaces for each stencil
  real, allocatable :: fn_u(:,:), fn_v(:,:)
  ! Total fluxes at cell interfaces
  real, allocatable :: f_u(:,:), f_v(:,:)

contains

  subroutine weno_init()

    if (.not. allocated(c)) then
      select case (weno_order)
      case (2)
        allocate(c(2,2))
        allocate(cb(1))
        allocate(wo(2))
        c(:,1) = [- 1.0 / 2.0,   3.0 / 2.0]
        c(:,2) = [  1.0 / 2.0,   1.0 / 2.0]
        wo(:)  = [  1.0 / 3.0,   2.0 / 3.0]
        cb(1)  = 1.0
      case (3)
        allocate(c(3,3))
        allocate(cb(2))
        allocate(wo(3))
        c(:,1) = [  1.0 /  3.0, - 7.0 /  6.0,  11.0 /  6.0]
        c(:,2) = [- 1.0 /  6.0,   5.0 /  6.0,   1.0 /  3.0]
        c(:,3) = [  1.0 /  3.0,   5.0 /  6.0, - 1.0 /  6.0]
        wo(:)  = [  1.0 / 10.0,   6.0 / 10.0,   3.0 / 10.0]
        cb(:)  = [ 13.0 / 12.0,   1.0 /  4.0]
      end select

      call parallel_allocate(fp_u, half_lon=.true.)
      call parallel_allocate(fn_u, half_lon=.true.)
      call parallel_allocate(f_u,  half_lon=.true.)
      call parallel_allocate(fp_v, half_lat=.true.)
      call parallel_allocate(fn_v, half_lat=.true.)
      call parallel_allocate(f_v,  half_lat=.true.)
    end if

  end subroutine weno_init

  subroutine weno_zonal_momentum_advection_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real u
    integer i, j

    ! U
    ! Calculate splitted flux at cell centers by using Lax-Friedrichs splitting.
    ! The splitting is for numerical stability by respecting upwind.
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        fp_u(i,j) = 0.5 * (state%u(i,j) + u_max) * state%iap%u(i,j)
        fn_u(i,j) = 0.5 * (state%u(i,j) - u_max) * state%iap%u(i,j)
      end do
    end do
    call parallel_fill_halo(fp_u, left_halo=.true., right_halo=.true.)
    call parallel_fill_halo(fn_u, left_halo=.true., right_halo=.true.)

    ! V
    ! Calculate splitted flux at cell centers by using Lax-Friedrichs splitting.
    ! The splitting is for numerical stability by respecting upwind.
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        u = 0.25 * (state%u(i-1,j) + state%u(i-1,j+1) + state%u(i,j) + state%u(i,j+1))
        fp_v(i,j) = 0.5 * (u + u_max) * state%iap%v(i,j)
        fn_v(i,j) = 0.5 * (u - u_max) * state%iap%v(i,j)
      end do
    end do
    call parallel_fill_halo(fp_v, left_halo=.true., right_halo=.true.)
    call parallel_fill_halo(fn_v, left_halo=.true., right_halo=.true.)

    select case (weno_order)
    case (2)
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          call weno_2nd_order_pass(fp_u(i-1,j), fp_u(i,j  ), fp_u(i+1,j), &
                                   fn_u(i,j  ), fn_u(i+1,j), fn_u(i+2,j), &
                                   f_u(i,j))
        end do
      end do

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          call weno_2nd_order_pass(fp_v(i-1,j), fp_v(i,j  ), fp_v(i+1,j), &
                                   fn_v(i,j  ), fn_v(i+1,j), fn_v(i+2,j), &
                                   f_v(i,j))
        end do
      end do
    case (3)
      ! Stencils:
      !
      ! 3rd order:
      !
      !   For positive flux:
      ! 
      !                               |___________S3__________|
      !                               |                       |
      !                       |___________S2__________|       |
      !                       |                       |       |
      !               |___________S1__________|       |       |
      !            ...|---o---|---o---|---o---|---o---|---o---|...
      !               |  i-2  |  i-1  |   i   |  i+1  |  i+2  |
      !                                      -|
      !                                     i+1/2
      !
      !   For negative flux:
      !
      !                       |___________S1__________|
      !                       |                       |
      !                       |       |___________S2__________|
      !                       |       |                       |
      !                       |       |       |___________S3__________|
      !                    ...|---o---|---o---|---o---|---o---|---o---|...
      !                       |  i-1  |   i   |  i+1  |  i+2  |  i+3  |
      !                                       |+
      !                                     i+1/2
    end select
    call parallel_fill_halo(f_u, left_halo=.true.)
    call parallel_fill_halo(f_v, left_halo=.true.)

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_adv_lon(i,j) = (f_u(i,j) - f_u(i-1,j) - &
          (state%u(i+1,j) - state%u(i-1,j)) * state%iap%u(i,j) * 0.25) / coef%half_dlon(j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%v_adv_lon(i,j) = (f_v(i,j) - f_v(i-1,j) - &
          (state%u(i,j) + state%u(i,j+1) - state%u(i-1,j) - state%u(i-1,j+1)) * state%iap%v(i,j) * 0.25) / coef%half_dlon(j)
      end do
    end do

  end subroutine weno_zonal_momentum_advection_operator

  subroutine weno_meridional_momentum_advection_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real v
    integer i, j

    ! U
    ! Calculate splitted flux at cell centers by using Lax-Friedrichs splitting.
    ! The splitting is for numerical stability by respecting upwind.
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        v = 0.25 * (state%v(i,j-1) + state%v(i,j) + state%v(i+1,j-1) + state%v(i+1,j))
        fp_u(i,j) = 0.5 * (v + v_max) * state%iap%u(i,j)
        fn_u(i,j) = 0.5 * (v - v_max) * state%iap%u(i,j)
      end do
    end do
    call parallel_fill_halo(fp_u, top_halo=.true., bottom_halo=.true.)
    call parallel_fill_halo(fn_u, top_halo=.true., bottom_halo=.true.)

    ! V
    ! Calculate splitted flux at cell centers by using Lax-Friedrichs splitting.
    ! The splitting is for numerical stability by respecting upwind.
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        fp_v(i,j) = 0.5 * (state%v(i,j) + v_max) * state%iap%v(i,j)
        fn_v(i,j) = 0.5 * (state%v(i,j) - v_max) * state%iap%v(i,j)
      end do
    end do
    call parallel_fill_halo(fp_v, top_halo=.true., bottom_halo=.true.)
    call parallel_fill_halo(fn_v, top_halo=.true., bottom_halo=.true.)

    select case (weno_order)
    case (2)
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          call weno_2nd_order_pass(fp_u(i,j-1), fp_u(i,j  ), fp_u(i,j+1), &
                                   fn_u(i,j  ), fn_u(i,j+1), fn_u(i,j+2), &
                                   f_u(i,j))
        end do
      end do

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          call weno_2nd_order_pass(fp_v(i,j-1), fp_v(i,j  ), fp_v(i,j+1), &
                                   fn_v(i,j  ), fn_v(i,j+1), fn_v(i,j+2), &
                                   f_v(i,j))
        end do
      end do
    end select
    call parallel_fill_halo(f_u, bottom_halo=.true.)
    call parallel_fill_halo(f_v, bottom_halo=.true.)

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_adv_lat(i,j) = (f_u(i,j) - f_u(i,j-1) - &
          (state%v(i-1,j) + state%v(i,j) - state%v(i-1,j-1) - state%v(i,j-1)) * state%iap%u(i,j) * 0.25) / coef%full_dlat(j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_adv_lat(i,j) = (f_v(i,j) - f_v(i,j-1) - &
          (state%v(i,j+1) - state%v(i,j-1)) * state%iap%v(i,j) * 0.25) / coef%half_dlat(j)
      end do
    end do

  end subroutine weno_meridional_momentum_advection_operator

  subroutine weno_final()

    if (allocated(c))  deallocate(c)
    if (allocated(wo)) deallocate(wo)
    if (allocated(fp_u)) deallocate(fp_u)
    if (allocated(fn_u)) deallocate(fn_u)
    if (allocated(fp_u)) deallocate(fp_u)
    if (allocated(fn_v)) deallocate(fn_v)
    if (allocated(f_u))  deallocate(f_u)
    if (allocated(f_v))  deallocate(f_v)

  end subroutine weno_final

  subroutine weno_2nd_order_pass(fp1, fp2, fp3, fn2, fn3, fn4, f)

    real, intent(in) :: fp1
    real, intent(in) :: fp2
    real, intent(in) :: fp3
    real, intent(in) :: fn2
    real, intent(in) :: fn3
    real, intent(in) :: fn4
    real, intent(out) :: f

    real fs(2)  ! Total fluxes at cell interfaces for each stencil
    real b(2)   ! Smooth indicators
    real w(2)   ! Stencil weights

    ! Positive flux at cell interfaces
    !
    !            |_______S2______|
    !    |_______S1______|
    ! ...|---o---|---o---|---o---|...
    !    |   1   |   2   |   3   |
    !                   -|
    !                    *
    ! - Calculate flux at interfaces for each stencil.
    fs(1) = c(1,1) * fp1 + c(2,1) * fp2
    fs(2) = c(1,2) * fp2 + c(2,2) * fp3
    ! - Calculate smooth indicators for each stencil regarding cell centers.
    b(1) = (fp2 - fp1)**2
    b(2) = (fp3 - fp2)**2
    ! - Calculate stencil linear combination weights considering smooth indicators.
    w(:) = wo(:) / (eps + b(:))**2
    w(:) = w(:) / sum(w(:))
    f = sum(w(:) * fs(:))
    ! Negative flux at cell interfaces
    !
    !              |_______S2______|
    !                      |_______S1______|
    !           ...|---o---|---o---|---o---|...
    !              |   2   |   3   |   4   |
    !                      |+
    !                      *
    !
    ! - Calculate flux at interfaces for each stencil.
    fs(1) = c(1,1) * fn4 + c(2,1) * fn3
    fs(2) = c(1,2) * fn3 + c(2,2) * fn2
    ! - Calculate smooth indicators for each stencil regarding cell centers.
    b(1) = (fn3 - fn4)**2
    b(2) = (fn2 - fn3)**2
    ! - Calculate stencil linear combination weights considering smooth indicators.
    w(:) = wo(:) / (eps + b(:))**2
    w(:) = w(:) / sum(w(:))
    f = f + sum(w(:) * fs(:))

  end subroutine weno_2nd_order_pass

end module weno_mod