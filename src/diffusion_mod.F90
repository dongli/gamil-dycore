module diffusion_mod

  use parallel_mod
  use params_mod
  use mesh_mod
  use types_mod
  use data_mod
  use reduce_mod
  use filter_mod

  implicit none

  private

  public diffusion_init
  public diffusion_run
  public diffusion_final

  real, allocatable :: ud_lon(:,:)
  real, allocatable :: ud_lat(:,:)
  real, allocatable :: vd_lon(:,:)
  real, allocatable :: vd_lat(:,:)
  real, allocatable :: gdd_lon(:,:)
  real, allocatable :: gdd_lat(:,:)

contains

  subroutine diffusion_init()

    if (.not. allocated(ud_lon)) call parallel_allocate(ud_lon, half_lon=.true., extended_halo=.true.)
    if (.not. allocated(ud_lat)) call parallel_allocate(ud_lat, half_lon=.true., extended_halo=.true.)
    if (.not. allocated(vd_lon)) call parallel_allocate(vd_lon, half_lat=.true., extended_halo=.true.)
    if (.not. allocated(vd_lat)) call parallel_allocate(vd_lat, half_lat=.true., extended_halo=.true.)
    if (.not. allocated(gdd_lon)) call parallel_allocate(gdd_lon, extended_halo=.true.)
    if (.not. allocated(gdd_lat)) call parallel_allocate(gdd_lat, extended_halo=.true.)

  end subroutine diffusion_init

  subroutine diffusion_run(dt, state)

    real, intent(in) :: dt
    type(state_type), intent(inout) :: state

    real reduced_tend(parallel%full_lon_start_idx:parallel%full_lon_end_idx)
    real sp, np
    integer i, j, k

    !
    ! ∂ F         1      ∂² F        1      ∂         ∂ F
    ! --- = 𝞶 ---------- ---- + 𝞶 --------- --- cos(φ) ---
    ! ∂ t     a² cos²(φ) ∂ λ²     a² cos(φ) ∂ φ        ∂ φ
    !

    ! U
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (full_reduce_factor(j) == 1) then
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          ud_lon(i,j) = (state%u(i+1,j) - 2 * state%u(i,j) + state%u(i-1,j)) / (0.5 * coef%full_dlon(j))**2
        end do
      else
        ud_lon(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = (full_reduced_state(j)%u(i+1,k,2) - 2 * full_reduced_state(j)%u(i,k,2) + full_reduced_state(j)%u(i-1,k,2)) / &
              (0.5 * coef%full_dlon(j) * full_reduce_factor(j))**2
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, ud_lon(:,j))
        end do
        call parallel_overlay_inner_halo(ud_lon(:,j), left_halo=.true.)
      end if
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        ud_lat(i,j) = ((state%u(i,j+1) - state%u(i,j)) * mesh%half_cos_lat(j) - &
                       (state%u(i,j) - state%u(i,j-1)) * mesh%half_cos_lat(j-1)) / &
                      (0.5 * coef%full_dlat(j))**2 * mesh%full_cos_lat(j)
      end do
    end do

    ! Do FFT filter.
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (filter_full_zonal_tend(j)) then
        call filter_run(ud_lon(:,j))
        call filter_run(ud_lat(:,j))
      end if
    end do

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%u(i,j) = state%u(i,j) + dt * diffusion_coef * (ud_lon(i,j) + ud_lat(i,j))
      end do
    end do

    call parallel_fill_halo(state%u, all_halo=.true.)

    ! V
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      if (half_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          vd_lon(i,j) = (state%v(i+1,j) - 2 * state%v(i,j) + state%v(i-1,j)) / (0.5 * coef%half_dlon(j))**2
        end do
      else
        vd_lon(:,j) = 0.0
        do k = 1, half_reduce_factor(j)
          do i = reduced_start_idx_at_half_lat(j), reduced_end_idx_at_half_lat(j)
            reduced_tend(i) = (half_reduced_state(j)%v(i+1,k,2) - 2 * half_reduced_state(j)%v(i,k,2) + half_reduced_state(j)%v(i-1,k,2)) / &
              (0.5 * coef%half_dlon(j) * half_reduce_factor(j))**2
          end do
          call append_reduced_tend_to_raw_tend_at_half_lat(j, k, reduced_tend, vd_lon(:,j))
        end do
        call parallel_overlay_inner_halo(vd_lon(:,j), left_halo=.true.)
      end if
    end do
    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vd_lat(i,j) = ((state%v(i,j+1) - state%v(i,j)) * mesh%full_cos_lat(j+1) - &
                       (state%v(i,j) - state%v(i,j-1)) * mesh%full_cos_lat(j)) / &
                      (0.5 * coef%half_dlat(j))**2 * mesh%half_cos_lat(j)
      end do
    end do
    if (parallel%has_south_pole) then
      j = parallel%half_lat_start_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vd_lat(i,j) = (state%v(i,j+1) - state%v(i,j)) * mesh%full_cos_lat(j+1) / &
          (0.5 * coef%half_dlat(j))**2 * mesh%half_cos_lat(j)
      end do
    end if
    if (parallel%has_north_pole) then
      j = parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        vd_lat(i,j) = (state%v(i,j) - state%v(i,j-1)) * mesh%full_cos_lat(j) / &
          (0.5 * coef%half_dlat(j))**2 * mesh%half_cos_lat(j)
      end do
    end if

    ! Do FFT filter.
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      if (filter_half_zonal_tend(j)) then
        call filter_run(vd_lon(:,j))
        call filter_run(vd_lat(:,j))
      end if
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%v(i,j) = state%v(i,j) + dt * diffusion_coef * (vd_lon(i,j) + vd_lat(i,j))
      end do
    end do

    call parallel_fill_halo(state%v, all_halo=.true.)

    ! H
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (full_reduce_factor(j) == 1) then
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          gdd_lon(i,j) = (state%gd(i+1,j) - 2 * state%gd(i,j) + state%gd(i-1,j)) / (0.5 * coef%full_dlon(j))**2
        end do
      else
        gdd_lon(:,j) = 0.0
        do k = 1, full_reduce_factor(j)
          do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
            reduced_tend(i) = (full_reduced_state(j)%gd(i+1,k,2) - 2 * full_reduced_state(j)%gd(i,k,2) + full_reduced_state(j)%gd(i-1,k,2)) / &
              (0.5 * coef%full_dlon(j) * full_reduce_factor(j))**2
          end do
          call append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, gdd_lon(:,j))
        end do
        call parallel_overlay_inner_halo(gdd_lon(:,j), left_halo=.true.)
      end if
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        gdd_lat(i,j) = ((state%gd(i,j+1) - state%gd(i,j)) * mesh%half_cos_lat(j) - &
                        (state%gd(i,j) - state%gd(i,j-1)) * mesh%half_cos_lat(j-1)) / &
                       (0.5 * coef%full_dlat(j))**2 * mesh%full_cos_lat(j)
      end do
    end do
    if (parallel%has_south_pole) then
      j = parallel%full_lat_start_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + state%gd(i,j+1) - state%gd(i,j)
      end do
      call parallel_zonal_sum(sp)
      sp = sp * mesh%half_cos_lat(j) / (0.5 * coef%full_dlat(j))**2 * mesh%full_cos_lat(j) / mesh%num_full_lon
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        gdd_lat(i,j) = sp
      end do
    end if
    if (parallel%has_north_pole) then
      j = parallel%full_lat_end_idx
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np + state%gd(i,j) - state%gd(i,j-1)
      end do
      call parallel_zonal_sum(np)
      np = np * mesh%half_cos_lat(j-1) / (0.5 * coef%full_dlat(j))**2 * mesh%full_cos_lat(j) / mesh%num_full_lon
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        gdd_lat(i,j) = np
      end do
    end if

    ! Do FFT filter.
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (filter_full_zonal_tend(j)) then
        call filter_run(gdd_lon(:,j))
        call filter_run(gdd_lat(:,j))
      end if
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%gd(i,j) = state%gd(i,j) + dt * diffusion_coef * (gdd_lon(i,j) + gdd_lat(i,j))
      end do
    end do

    call parallel_fill_halo(state%gd, all_halo=.true.)

    call iap_transform(state)

  end subroutine diffusion_run

  subroutine diffusion_final()

    if (allocated(ud_lon)) deallocate(ud_lon)
    if (allocated(ud_lat)) deallocate(ud_lat)
    if (allocated(vd_lon)) deallocate(vd_lon)
    if (allocated(vd_lat)) deallocate(vd_lat)
    if (allocated(gdd_lon)) deallocate(gdd_lon)
    if (allocated(gdd_lat)) deallocate(gdd_lat)

  end subroutine diffusion_final

end module diffusion_mod
