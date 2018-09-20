module types_mod

  use mesh_mod
  use parallel_mod

  implicit none

  private

  public allocate_data
  public deallocate_data
  public iap_transform
  public copy_state
  public average_state
  public coef_type
  public state_type
  public static_type
  public iap_type
  public tend_type

  type coef_type
    ! Coriolis coefficient at full meridional grids
    real, allocatable :: cori(:)
    ! Curvature coefficient at full meridional grids 
    real, allocatable :: curv(:)
    ! Zonal difference coefficient at full/half meridional grids
    real, allocatable :: full_dlon(:)
    real, allocatable :: half_dlon(:)
    ! Meridional difference coefficient at full/half meridional grids
    real, allocatable :: full_dlat(:)
    real, allocatable :: half_dlat(:)
  end type coef_type

  ! IAP transformed variables
  type iap_type
    real, allocatable :: u_a(:,:)
    real, allocatable :: v_a(:,:)
    real, allocatable :: u_c(:,:)
    real, allocatable :: v_c(:,:)
    real, allocatable :: gd(:,:)
  end type iap_type

  type state_type
    real, allocatable :: u_a(:,:)
    real, allocatable :: v_a(:,:)
    real, allocatable :: u_c(:,:)
    real, allocatable :: v_c(:,:)
    real, allocatable :: gd(:,:) ! Geopotential depth
    type(iap_type) iap
  end type state_type

  type static_type
    real, allocatable :: ghs(:,:) ! Surface geopotential
  end type static_type

  type tend_type
    real, allocatable :: u_adv_lon_a(:,:)
    real, allocatable :: u_adv_lat_a(:,:)
    real, allocatable :: v_adv_lon_a(:,:)
    real, allocatable :: v_adv_lat_a(:,:)
    real, allocatable :: fu_a(:,:)
    real, allocatable :: fv_a(:,:)
    real, allocatable :: cu_a(:,:)
    real, allocatable :: cv_a(:,:)
    real, allocatable :: u_pgf_a(:,:)
    real, allocatable :: v_pgf_a(:,:)
    real, allocatable :: u_adv_lon_c(:,:)
    real, allocatable :: u_adv_lat_c(:,:)
    real, allocatable :: v_adv_lon_c(:,:)
    real, allocatable :: v_adv_lat_c(:,:)
    real, allocatable :: fu_c(:,:)
    real, allocatable :: fv_c(:,:)
    real, allocatable :: cu_c(:,:)
    real, allocatable :: cv_c(:,:)
    real, allocatable :: u_pgf_c(:,:)
    real, allocatable :: v_pgf_c(:,:)
    real, allocatable :: mass_div_lon(:,:)
    real, allocatable :: mass_div_lat(:,:)
    real, allocatable :: du_a(:,:)
    real, allocatable :: dv_a(:,:)
    real, allocatable :: du_c(:,:)
    real, allocatable :: dv_c(:,:)
    real, allocatable :: dgd(:,:)
  end type tend_type

  interface allocate_data
    module procedure allocate_coef_data
    module procedure allocate_static_data
    module procedure allocate_state_data
    module procedure allocate_tend_data
  end interface allocate_data

  interface deallocate_data
    module procedure deallocate_coef_data
    module procedure deallocate_static_data
    module procedure deallocate_state_data
    module procedure deallocate_tend_data
  end interface deallocate_data

contains

  subroutine allocate_coef_data(coef)

    type(coef_type), intent(out) :: coef

    allocate(coef%cori(mesh%num_full_lat))
    allocate(coef%curv(mesh%num_full_lat))
    allocate(coef%full_dlon(mesh%num_full_lat))
    allocate(coef%half_dlon(mesh%num_half_lat))
    allocate(coef%full_dlat(mesh%num_full_lat))
    allocate(coef%half_dlat(mesh%num_half_lat))

  end subroutine allocate_coef_data

  subroutine allocate_static_data(static)

    type(static_type), intent(out) :: static

    if (.not. allocated(static%ghs)) call parallel_allocate(static%ghs, extended_halo=.true.)

  end subroutine allocate_static_data

  subroutine allocate_state_data(state)

    type(state_type), intent(out) :: state

    if (.not. allocated(state%u_a))         call parallel_allocate(state%u_a,                            extended_halo=.true.)
    if (.not. allocated(state%v_a))         call parallel_allocate(state%v_a,                            extended_halo=.true.)
    if (.not. allocated(state%u_c))         call parallel_allocate(state%u_c,           half_lon=.true., extended_halo=.true.)
    if (.not. allocated(state%v_c))         call parallel_allocate(state%v_c,           half_lat=.true., extended_halo=.true.)
    if (.not. allocated(state%gd))          call parallel_allocate(state%gd,                             extended_halo=.true.)
    if (.not. allocated(state%iap%u_a))     call parallel_allocate(state%iap%u_a,                        extended_halo=.true.)
    if (.not. allocated(state%iap%v_a))     call parallel_allocate(state%iap%v_a,                        extended_halo=.true.)
    if (.not. allocated(state%iap%u_c))     call parallel_allocate(state%iap%u_c,       half_lon=.true., extended_halo=.true.)
    if (.not. allocated(state%iap%v_c))     call parallel_allocate(state%iap%v_c,       half_lat=.true., extended_halo=.true.)
    if (.not. allocated(state%iap%gd))      call parallel_allocate(state%iap%gd,                         extended_halo=.true.)

  end subroutine allocate_state_data

  subroutine allocate_tend_data(tend)

    type(tend_type), intent(out) :: tend

    if (.not. allocated(tend%u_adv_lon_a))  call parallel_allocate(tend%u_adv_lon_a,                     extended_halo=.true.)
    if (.not. allocated(tend%u_adv_lat_a))  call parallel_allocate(tend%u_adv_lat_a,                     extended_halo=.true.)
    if (.not. allocated(tend%fv_a))         call parallel_allocate(tend%fv_a,                            extended_halo=.true.)
    if (.not. allocated(tend%cv_a))         call parallel_allocate(tend%cv_a,                            extended_halo=.true.)
    if (.not. allocated(tend%u_pgf_a))      call parallel_allocate(tend%u_pgf_a,                         extended_halo=.true.)
    if (.not. allocated(tend%v_adv_lon_a))  call parallel_allocate(tend%v_adv_lon_a,                     extended_halo=.true.)
    if (.not. allocated(tend%v_adv_lat_a))  call parallel_allocate(tend%v_adv_lat_a,                     extended_halo=.true.)
    if (.not. allocated(tend%fu_a))         call parallel_allocate(tend%fu_a,                            extended_halo=.true.)
    if (.not. allocated(tend%cu_a))         call parallel_allocate(tend%cu_a,                            extended_halo=.true.)
    if (.not. allocated(tend%v_pgf_a))      call parallel_allocate(tend%v_pgf_a,                         extended_halo=.true.)
    if (.not. allocated(tend%u_adv_lon_c))  call parallel_allocate(tend%u_adv_lon_c,    half_lon=.true., extended_halo=.true.)
    if (.not. allocated(tend%u_adv_lat_c))  call parallel_allocate(tend%u_adv_lat_c,    half_lon=.true., extended_halo=.true.)
    if (.not. allocated(tend%fv_c))         call parallel_allocate(tend%fv_c,           half_lon=.true., extended_halo=.true.)
    if (.not. allocated(tend%cv_c))         call parallel_allocate(tend%cv_c,           half_lon=.true., extended_halo=.true.)
    if (.not. allocated(tend%u_pgf_c))      call parallel_allocate(tend%u_pgf_c,        half_lon=.true., extended_halo=.true.)
    if (.not. allocated(tend%v_adv_lon_c))  call parallel_allocate(tend%v_adv_lon_c,    half_lat=.true., extended_halo=.true.)
    if (.not. allocated(tend%v_adv_lat_c))  call parallel_allocate(tend%v_adv_lat_c,    half_lat=.true., extended_halo=.true.)
    if (.not. allocated(tend%fu_c))         call parallel_allocate(tend%fu_c,           half_lat=.true., extended_halo=.true.)
    if (.not. allocated(tend%cu_c))         call parallel_allocate(tend%cu_c,           half_lat=.true., extended_halo=.true.)
    if (.not. allocated(tend%v_pgf_c))      call parallel_allocate(tend%v_pgf_c,        half_lat=.true., extended_halo=.true.)
    if (.not. allocated(tend%mass_div_lon)) call parallel_allocate(tend%mass_div_lon,                    extended_halo=.true.)
    if (.not. allocated(tend%mass_div_lat)) call parallel_allocate(tend%mass_div_lat,                    extended_halo=.true.)
    if (.not. allocated(tend%du_a))         call parallel_allocate(tend%du_a,                            extended_halo=.true.)
    if (.not. allocated(tend%dv_a))         call parallel_allocate(tend%dv_a,                            extended_halo=.true.)
    if (.not. allocated(tend%du_c))         call parallel_allocate(tend%du_c,           half_lon=.true., extended_halo=.true.)
    if (.not. allocated(tend%dv_c))         call parallel_allocate(tend%dv_c,           half_lat=.true., extended_halo=.true.)
    if (.not. allocated(tend%dgd))          call parallel_allocate(tend%dgd,                             extended_halo=.true.)

  end subroutine allocate_tend_data

  subroutine deallocate_coef_data(coef)

    type(coef_type), intent(inout) :: coef

    if (allocated(coef%cori)) deallocate(coef%cori)
    if (allocated(coef%curv)) deallocate(coef%curv)
    if (allocated(coef%full_dlon)) deallocate(coef%full_dlon)
    if (allocated(coef%half_dlon)) deallocate(coef%half_dlon)
    if (allocated(coef%full_dlat)) deallocate(coef%full_dlat)
    if (allocated(coef%half_dlat)) deallocate(coef%half_dlat)

  end subroutine deallocate_coef_data

  subroutine deallocate_static_data(static)

    type(static_type), intent(inout) :: static

    if (allocated(static%ghs)) deallocate(static%ghs)

  end subroutine deallocate_static_data

  subroutine deallocate_state_data(state)

    type(state_type), intent(inout) :: state

    if (allocated(state%u_a)) deallocate(state%u_a)
    if (allocated(state%v_a)) deallocate(state%v_a)
    if (allocated(state%u_c)) deallocate(state%u_c)
    if (allocated(state%v_c)) deallocate(state%v_c)
    if (allocated(state%gd)) deallocate(state%gd)

    if (allocated(state%iap%u_a)) deallocate(state%iap%u_a)
    if (allocated(state%iap%v_a)) deallocate(state%iap%v_a)
    if (allocated(state%iap%u_c)) deallocate(state%iap%u_c)
    if (allocated(state%iap%v_c)) deallocate(state%iap%v_c)
    if (allocated(state%iap%gd)) deallocate(state%iap%gd)

  end subroutine deallocate_state_data

  subroutine deallocate_tend_data(tend)

    type(tend_type), intent(inout) :: tend

    if (allocated(tend%u_adv_lon_a)) deallocate(tend%u_adv_lon_a)
    if (allocated(tend%u_adv_lat_a)) deallocate(tend%u_adv_lat_a)
    if (allocated(tend%v_adv_lon_a)) deallocate(tend%v_adv_lon_a)
    if (allocated(tend%v_adv_lat_a)) deallocate(tend%v_adv_lat_a)
    if (allocated(tend%fu_a)) deallocate(tend%fu_a)
    if (allocated(tend%fv_a)) deallocate(tend%fv_a)
    if (allocated(tend%cu_a)) deallocate(tend%cu_a)
    if (allocated(tend%cv_a)) deallocate(tend%cv_a)
    if (allocated(tend%u_pgf_a)) deallocate(tend%u_pgf_a)
    if (allocated(tend%v_pgf_a)) deallocate(tend%v_pgf_a)
    if (allocated(tend%u_adv_lon_c)) deallocate(tend%u_adv_lon_c)
    if (allocated(tend%u_adv_lat_c)) deallocate(tend%u_adv_lat_c)
    if (allocated(tend%v_adv_lon_c)) deallocate(tend%v_adv_lon_c)
    if (allocated(tend%v_adv_lat_c)) deallocate(tend%v_adv_lat_c)
    if (allocated(tend%fu_c)) deallocate(tend%fu_c)
    if (allocated(tend%fv_c)) deallocate(tend%fv_c)
    if (allocated(tend%cu_c)) deallocate(tend%cu_c)
    if (allocated(tend%cv_c)) deallocate(tend%cv_c)
    if (allocated(tend%u_pgf_c)) deallocate(tend%u_pgf_c)
    if (allocated(tend%v_pgf_c)) deallocate(tend%v_pgf_c)
    if (allocated(tend%mass_div_lon)) deallocate(tend%mass_div_lon)
    if (allocated(tend%mass_div_lat)) deallocate(tend%mass_div_lat)
    if (allocated(tend%du_a)) deallocate(tend%du_a)
    if (allocated(tend%dv_a)) deallocate(tend%dv_a)
    if (allocated(tend%du_c)) deallocate(tend%du_c)
    if (allocated(tend%dv_c)) deallocate(tend%dv_c)
    if (allocated(tend%dgd)) deallocate(tend%dgd)

  end subroutine deallocate_tend_data

  subroutine iap_transform(state)

    type(state_type), intent(inout) :: state

    integer i, j

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%gd(i,j) = sqrt(state%gd(i,j))
      end do
    end do

    call parallel_fill_halo(state%iap%gd(:,:), all_halo=.true.)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%u_a(i,j) = state%iap%gd(i,j) * state%u_a(i,j)
      end do
    end do

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%v_a(i,j) = state%iap%gd(i,j) * state%v_a(i,j)
      end do
    end do

    call parallel_fill_halo(state%iap%u_a(:,:), all_halo=.true.)
    call parallel_fill_halo(state%iap%v_a(:,:), all_halo=.true.)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%iap%u_c(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i+1,j)) * state%u_c(i,j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%v_c(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%v_c(i,j)
      end do
    end do

    call parallel_fill_halo(state%iap%u_c(:,:), all_halo=.true.)
    call parallel_fill_halo(state%iap%v_c(:,:), all_halo=.true.)

  end subroutine iap_transform

  subroutine copy_state(state1, state2)

    type(state_type), intent(in) :: state1
    type(state_type), intent(inout) :: state2

    state2%u_a = state1%u_a
    state2%iap%u_a = state1%iap%u_a
    state2%v_a = state1%v_a
    state2%iap%v_a = state1%iap%v_a
    state2%u_c = state1%u_c
    state2%iap%u_c = state1%iap%u_c
    state2%v_c = state1%v_c
    state2%iap%v_c = state1%iap%v_c
    state2%gd = state1%gd
    state2%iap%gd = state1%iap%gd

  end subroutine copy_state

  subroutine average_state(state1, state2, state3)

    type(state_type), intent(in) :: state1
    type(state_type), intent(in) :: state2
    type(state_type), intent(inout) :: state3

    state3%u_a = (state1%u_a + state2%u_a) * 0.5
    state3%iap%u_a = (state1%iap%u_a + state2%iap%u_a) * 0.5
    state3%v_a = (state1%v_a + state2%v_a) * 0.5
    state3%iap%v_a = (state1%iap%v_a + state2%iap%v_a) * 0.5
    state3%u_c = (state1%u_c + state2%u_c) * 0.5
    state3%iap%u_c = (state1%iap%u_c + state2%iap%u_c) * 0.5
    state3%v_c = (state1%v_c + state2%v_c) * 0.5
    state3%iap%v_c = (state1%iap%v_c + state2%iap%v_c) * 0.5
    state3%gd = (state1%gd + state2%gd) * 0.5
    state3%iap%gd = (state1%iap%gd + state2%iap%gd) * 0.5

  end subroutine average_state

end module types_mod
