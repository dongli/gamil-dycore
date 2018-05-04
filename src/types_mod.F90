module types_mod

  use mesh_mod
  use parallel_mod

  implicit none

  private

  public allocate_data
  public deallocate_data
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
    real, allocatable :: u(:,:)
    real, allocatable :: v(:,:)
    real, allocatable :: gd(:,:)
    real, allocatable :: reduced_u(:,:)
    real, allocatable :: reduced_v(:,:,:) ! The third dimension is due to v is at half meridional grids.
    real, allocatable :: reduced_gd(:,:)
  end type iap_type

  type state_type
    real, allocatable :: u(:,:)
    real, allocatable :: v(:,:)
    real, allocatable :: gd(:,:) ! Geopotential depth
    type(iap_type) iap
    ! Zonal maximum CFL number
    real, allocatable :: max_cfl(:)
    integer, allocatable :: reduce_factor(:)
    real, allocatable :: reduced_u(:,:)
    real, allocatable :: reduced_gd(:,:)
  end type state_type

  type static_type
    real, allocatable :: ghs(:,:) ! Surface geopotential
    real, allocatable :: reduced_ghs(:,:)
  end type static_type

  type tend_type
    real, allocatable :: u_adv_lon(:,:)
    real, allocatable :: u_adv_lat(:,:)
    real, allocatable :: v_adv_lon(:,:)
    real, allocatable :: v_adv_lat(:,:)
    real, allocatable :: fu(:,:)
    real, allocatable :: fv(:,:)
    real, allocatable :: cu(:,:)
    real, allocatable :: cv(:,:)
    real, allocatable :: u_pgf(:,:)
    real, allocatable :: v_pgf(:,:)
    real, allocatable :: mass_div_lon(:,:)
    real, allocatable :: mass_div_lat(:,:)
    real, allocatable :: du(:,:)
    real, allocatable :: dv(:,:)
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

    if (.not. allocated(static%ghs)) call parallel_allocate(static%ghs)
    if (.not. allocated(static%reduced_ghs)) call parallel_allocate(static%reduced_ghs)

  end subroutine allocate_static_data

  subroutine allocate_state_data(state)

    type(state_type), intent(out) :: state

    if (.not. allocated(state%u)) call parallel_allocate(state%u, half_lon=.true.)
    if (.not. allocated(state%v)) call parallel_allocate(state%v, half_lat=.true.)
    if (.not. allocated(state%gd)) call parallel_allocate(state%gd)
    if (.not. allocated(state%max_cfl)) call parallel_allocate(state%max_cfl, full_lat=.true.)
    if (.not. allocated(state%reduce_factor)) call parallel_allocate(state%reduce_factor, full_lat=.true.)
    if (.not. allocated(state%reduced_u)) call parallel_allocate(state%reduced_u, half_lon=.true.)
    if (.not. allocated(state%reduced_gd)) call parallel_allocate(state%reduced_gd)

    if (.not. allocated(state%iap%u)) call parallel_allocate(state%iap%u, half_lon=.true.)
    if (.not. allocated(state%iap%v)) call parallel_allocate(state%iap%v, half_lat=.true.)
    if (.not. allocated(state%iap%gd)) call parallel_allocate(state%iap%gd)
    if (.not. allocated(state%iap%reduced_u)) call parallel_allocate(state%iap%reduced_u, half_lon=.true.)
    if (.not. allocated(state%iap%reduced_v)) call parallel_allocate(state%iap%reduced_v, dim=3, size=2)
    if (.not. allocated(state%iap%reduced_gd)) call parallel_allocate(state%iap%reduced_gd)

  end subroutine allocate_state_data

  subroutine allocate_tend_data(tend)

    type(tend_type), intent(out) :: tend

    if (.not. allocated(tend%u_adv_lon)) call parallel_allocate(tend%u_adv_lon, half_lon=.true.)
    if (.not. allocated(tend%u_adv_lat)) call parallel_allocate(tend%u_adv_lat, half_lon=.true.)
    if (.not. allocated(tend%fv)) call parallel_allocate(tend%fv, half_lon=.true.)
    if (.not. allocated(tend%cv)) call parallel_allocate(tend%cv, half_lon=.true.)
    if (.not. allocated(tend%u_pgf)) call parallel_allocate(tend%u_pgf, half_lon=.true.)
    if (.not. allocated(tend%v_adv_lon)) call parallel_allocate(tend%v_adv_lon, half_lat=.true.)
    if (.not. allocated(tend%v_adv_lat)) call parallel_allocate(tend%v_adv_lat, half_lat=.true.)
    if (.not. allocated(tend%fu)) call parallel_allocate(tend%fu, half_lat=.true.)
    if (.not. allocated(tend%cu)) call parallel_allocate(tend%cu, half_lat=.true.)
    if (.not. allocated(tend%v_pgf)) call parallel_allocate(tend%v_pgf, half_lat=.true.)
    if (.not. allocated(tend%mass_div_lon)) call parallel_allocate(tend%mass_div_lon)
    if (.not. allocated(tend%mass_div_lat)) call parallel_allocate(tend%mass_div_lat)
    if (.not. allocated(tend%du)) call parallel_allocate(tend%du, half_lon=.true.)
    if (.not. allocated(tend%dv)) call parallel_allocate(tend%dv, half_lat=.true.)
    if (.not. allocated(tend%dgd)) call parallel_allocate(tend%dgd)

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
    if (allocated(static%reduced_ghs)) deallocate(static%reduced_ghs)

  end subroutine deallocate_static_data

  subroutine deallocate_state_data(state)

    type(state_type), intent(inout) :: state

    if (allocated(state%u)) deallocate(state%u)
    if (allocated(state%v)) deallocate(state%v)
    if (allocated(state%gd)) deallocate(state%gd)
    if (allocated(state%max_cfl)) deallocate(state%max_cfl)
    if (allocated(state%reduce_factor)) deallocate(state%reduce_factor)
    if (allocated(state%reduced_u)) deallocate(state%reduced_u)
    if (allocated(state%reduced_gd)) deallocate(state%reduced_gd)

    if (allocated(state%iap%u)) deallocate(state%iap%u)
    if (allocated(state%iap%v)) deallocate(state%iap%v)
    if (allocated(state%iap%gd)) deallocate(state%iap%gd)
    if (allocated(state%iap%reduced_u)) deallocate(state%iap%reduced_u)
    if (allocated(state%iap%reduced_v)) deallocate(state%iap%reduced_v)
    if (allocated(state%iap%reduced_gd)) deallocate(state%iap%reduced_gd)

  end subroutine deallocate_state_data

  subroutine deallocate_tend_data(tend)

    type(tend_type), intent(inout) :: tend

    if (allocated(tend%u_adv_lon)) deallocate(tend%u_adv_lon)
    if (allocated(tend%u_adv_lat)) deallocate(tend%u_adv_lat)
    if (allocated(tend%v_adv_lon)) deallocate(tend%v_adv_lon)
    if (allocated(tend%v_adv_lat)) deallocate(tend%v_adv_lat)
    if (allocated(tend%fu)) deallocate(tend%fu)
    if (allocated(tend%fv)) deallocate(tend%fv)
    if (allocated(tend%cu)) deallocate(tend%cu)
    if (allocated(tend%cv)) deallocate(tend%cv)
    if (allocated(tend%u_pgf)) deallocate(tend%u_pgf)
    if (allocated(tend%v_pgf)) deallocate(tend%v_pgf)
    if (allocated(tend%mass_div_lon)) deallocate(tend%mass_div_lon)
    if (allocated(tend%mass_div_lat)) deallocate(tend%mass_div_lat)
    if (allocated(tend%du)) deallocate(tend%du)
    if (allocated(tend%dv)) deallocate(tend%dv)
    if (allocated(tend%dgd)) deallocate(tend%dgd)

  end subroutine deallocate_tend_data

  subroutine copy_state(state1, state2)

    type(state_type), intent(in) :: state1
    type(state_type), intent(inout) :: state2

    state2%u = state1%u
    state2%iap%u = state1%iap%u
    state2%v = state1%v
    state2%iap%v = state1%iap%v
    state2%gd = state1%gd
    state2%iap%gd = state1%iap%gd

  end subroutine copy_state

  subroutine average_state(state1, state2, state3)

    type(state_type), intent(in) :: state1
    type(state_type), intent(in) :: state2
    type(state_type), intent(inout) :: state3

    state3%u = (state1%u + state2%u) * 0.5
    state3%iap%u = (state1%iap%u + state2%iap%u) * 0.5
    state3%v = (state1%v + state2%v) * 0.5
    state3%iap%v = (state1%iap%v + state2%iap%v) * 0.5
    state3%gd = (state1%gd + state2%gd) * 0.5
    state3%iap%gd = (state1%iap%gd + state2%iap%gd) * 0.5

  end subroutine average_state

end module types_mod
