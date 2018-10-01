module diffusion_mod

  use parallel_mod
  use params_mod
  use mesh_mod
  use types_mod
  use data_mod
  use diag_mod
  use filter_mod

  implicit none

  private

  public diffusion_init
  public divergence_diffusion
  public ordinary_diffusion
  public diffusion_final

  real, allocatable :: ud(:,:)
  real, allocatable :: vd(:,:)
  real, allocatable :: gdd(:,:)

  real, allocatable :: u(:,:)
  real, allocatable :: v(:,:)
  real, allocatable :: gd(:,:)

contains

  subroutine diffusion_init()

    if (.not. allocated(ud))  call parallel_allocate(ud, half_lon=.true.)
    if (.not. allocated(vd))  call parallel_allocate(vd, half_lat=.true.)
    if (.not. allocated(gdd)) call parallel_allocate(gdd)
    if (.not. allocated(u))  call parallel_allocate(u, half_lon=.true.)
    if (.not. allocated(v))  call parallel_allocate(v, half_lat=.true.)
    if (.not. allocated(gd)) call parallel_allocate(gd)

  end subroutine diffusion_init

  subroutine divergence_diffusion(dt, diag, state)

    real, intent(in) :: dt
    type(diag_type), intent(in) :: diag
    type(state_type), intent(inout) :: state

    real, parameter :: vd = 1.0e5
    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%u(i,j) = state%u(i,j) + dt * vd * (diag%div(i+1,j) - diag%div(i,j)) / coef%full_dlon(j)
      end do
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%iap%u(i,j) = state%u(i,j) * 0.5 * (state%iap%gd(i,j) + state%iap%gd(i+1,j))
      end do
    end do
    do j = parallel%half_lon_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%v(i,j) = state%v(i,j) + dt * vd * (diag%div(i,j+1) - diag%div(i,j)) / radius / mesh%dlat
      end do
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%v(i,j) = state%v(i,j) * 0.5 * (state%iap%gd(i,j) + state%iap%gd(i,j+1))
      end do
    end do

    call parallel_fill_halo(state%iap%u(:,:), all_halo=.true.)
    call parallel_fill_halo(state%iap%v(:,:), all_halo=.true.)
    call parallel_fill_halo(state%u(:,:), all_halo=.true.)
    call parallel_fill_halo(state%v(:,:), all_halo=.true.)

  end subroutine divergence_diffusion

  subroutine ordinary_diffusion(dt, state)

    real, intent(in) :: dt
    type(state_type), intent(inout) :: state

    real sp, np
    integer i, j, order, sign

    u(:,:) = state%u(:,:)
    v(:,:) = state%v(:,:)
    gd(:,:) = state%gd(:,:)

    ! Scalar diffusion:
    !
    ! 2nd order:
    !
    !   ∂ F         1      ∂² F        1      ∂          ∂ F
    !   --- = 𝞶 ---------- ---- + 𝞶 --------- --- cos(φ) ---
    !   ∂ t     a² cos²(φ) ∂ λ²     a² cos(φ) ∂ φ        ∂ φ
    !
    !
    ! Vector diffusion:
    !
    ! 2nd order:
    !
    !            u        2 sin𝞿   ∂ v
    ! ∇² u - --------- + --------- ---
    !        a² cos²𝞿    a² cos²𝞿  ∂ 𝞴
    !
    !            v        2 sin𝞿   ∂ u
    ! ∇² v - --------- - --------- ---
    !        a² cos²𝞿    a² cos²𝞿  ∂ 𝞴

    do order = 1, diffusion_order / 2
      ! H
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          gdd(i,j) = sign * (gd(i+1,j) - 2 * gd(i,j) + gd(i-1,j)) / coef%full_dlon(j)**2 + &
                     ((gd(i,j+1) - gd(i,j  )) * mesh%half_cos_lat(j  ) - &
                      (gd(i,j  ) - gd(i,j-1)) * mesh%half_cos_lat(j-1)) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j)
        end do
      end do
      if (parallel%has_south_pole) then
        j = parallel%full_lat_south_pole_idx
        sp = 0.0
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          sp = sp + gd(i,j+1) - gd(i,j)
        end do
        call parallel_zonal_sum(sp)
        gdd(:,j) = sp * mesh%half_cos_lat(j) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j) / mesh%num_full_lon
      end if
      if (parallel%has_north_pole) then
        j = parallel%full_lat_north_pole_idx
        np = 0.0
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          np = np - (gd(i,j) - gd(i,j-1))
        end do
        call parallel_zonal_sum(np)
        gdd(:,j) = np * mesh%half_cos_lat(j-1) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j) / mesh%num_full_lon
      end if
      ! U
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          ud(i,j) = (u(i+1,j) - 2 * u(i,j) + u(i-1,j)) / coef%full_dlon(j)**2 + &
                    ((u(i,j+1) - u(i,j  )) * mesh%half_cos_lat(j  ) - &
                     (u(i,j  ) - u(i,j-1)) * mesh%half_cos_lat(j-1)) / coef%full_dlat(j)**2 * mesh%full_cos_lat(j) !- &
                    ! (u(i,j) - mesh%full_sin_lat(j) * (v(i+1,j-1) + v(i+1,j) - v(i,j-1) - v(i,j)) / mesh%dlon) / &
                    !  radius**2 / mesh%full_cos_lat(j)**2
        end do
      end do
      ! V
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          vd(i,j) = (v(i+1,j) - 2 * v(i,j) + v(i-1,j)) / coef%half_dlon(j)**2
        end do
      end do
      do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          vd(i,j) = vd(i,j) + ((v(i,j+1) - v(i,j  )) * mesh%full_cos_lat(j+1) - &
                               (v(i,j  ) - v(i,j-1)) * mesh%full_cos_lat(j  )) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j) !- &
                              ! (v(i,j) + mesh%half_sin_lat(j) * (u(i,j) + u(i,j+1) - u(i-1,j) - u(i-1,j+1)) / mesh%dlon) / &
                              !  radius**2 / mesh%half_cos_lat(j)**2
        end do
      end do
      if (parallel%has_south_pole) then
        j = parallel%half_lat_south_pole_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          vd(i,j) = vd(i,j) + (v(i,j+1) - v(i,j)) * mesh%full_cos_lat(j+1) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j) !- &
                              ! (v(i,j) + mesh%half_sin_lat(j) * (u(i,j+1) - u(i-1,j+1)) / mesh%dlon) / radius**2 / mesh%half_cos_lat(j)**2
        end do
      end if
      if (parallel%has_north_pole) then
        j = parallel%half_lat_north_pole_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          vd(i,j) = vd(i,j) - (v(i,j) - v(i,j-1)) * mesh%full_cos_lat(j) / coef%half_dlat(j)**2 * mesh%half_cos_lat(j) !- &
                              ! (v(i,j) + mesh%half_sin_lat(j) * (u(i,j) - u(i-1,j)) / mesh%dlon) / radius**2 / mesh%half_cos_lat(j)**2
        end do
      end if
      if (order /= diffusion_order / 2) then
        gd(:,:) = gdd(:,:)
        u(:,:) = ud(:,:)
        v(:,:) = vd(:,:)
      end if
    end do

    ! Do FFT filter.
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (filter_full_zonal_tend(j)) then
        call filter_array_at_full_lat(j, gdd(:,j))
        call filter_array_at_full_lat(j, ud(:,j))
      end if
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      if (filter_half_zonal_tend(j)) then
        call filter_array_at_half_lat(j, vd(:,j))
      end if
    end do

    sign = (-1)**(diffusion_order / 2 + 1)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%gd(i,j) = state%gd(i,j) + sign * dt * diffusion_coef * gdd(i,j)
      end do
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%u(i,j) = state%u(i,j) + sign * dt * diffusion_coef * ud(i,j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%v(i,j) = state%v(i,j) + sign * dt * diffusion_coef * vd(i,j)
      end do
    end do

    call parallel_fill_halo(state%gd, all_halo=.true.)
    call parallel_fill_halo(state%u, all_halo=.true.)
    call parallel_fill_halo(state%v, all_halo=.true.)

    call iap_transform(state)

  end subroutine ordinary_diffusion

  subroutine diffusion_final()

    if (allocated(ud))  deallocate(ud)
    if (allocated(vd))  deallocate(vd)
    if (allocated(gdd)) deallocate(gdd)
    if (allocated(u))  deallocate(u)
    if (allocated(v))  deallocate(v)
    if (allocated(gd)) deallocate(gd)

  end subroutine diffusion_final

end module diffusion_mod
