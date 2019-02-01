module types_mod

  use mesh_mod
  use parallel_mod

  implicit none

  private

  public allocate_data
  public deallocate_data
  public copy_data
  public zero_data
  public add_data
  public sub_data
  public scale_data
  public inner_product
  public iap_transform
  public coef_type
  public state_type
  public static_type
  public iap_type
  public tend_type

  type coef_type
    ! Coriolis coefficient at full/half meridional grids
    real, allocatable :: full_f(:)
    real, allocatable :: half_f(:)
    ! Curvature coefficient at full/half meridional grids 
    real, allocatable :: full_c(:)
    real, allocatable :: half_c(:)
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
  end type iap_type

  type state_type
    real, allocatable :: u(:,:)
    real, allocatable :: v(:,:)
    real, allocatable :: gd(:,:) ! Geopotential depth
    type(iap_type) iap
  end type state_type

  type static_type
    real, allocatable :: ghs(:,:) ! Surface geopotential
  end type static_type

  type tend_type
    real, allocatable :: u_adv_lon(:,:)
    real, allocatable :: u_adv_lat(:,:)
    real, allocatable :: v_adv_lon(:,:)
    real, allocatable :: v_adv_lat(:,:)
    real, allocatable :: fu(:,:)
    real, allocatable :: fv(:,:)
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

  interface copy_data
    module procedure copy_state_data
    module procedure copy_tend_data
  end interface copy_data

  interface zero_data
    module procedure zero_tend_data
  end interface zero_data

  interface add_data
    module procedure add_tend_data
  end interface add_data

  interface sub_data
    module procedure sub_tend_data
  end interface sub_data

  interface scale_data
    module procedure scale_tend_data
  end interface scale_data

  interface inner_product
    module procedure inner_product_tend_tend
    module procedure inner_product_tend_state
  end interface inner_product

contains

  subroutine allocate_coef_data(coef)

    type(coef_type), intent(out) :: coef

    allocate(coef%full_f(mesh%num_full_lat))
    allocate(coef%half_f(mesh%num_half_lat))
    allocate(coef%full_c(mesh%num_full_lat))
    allocate(coef%half_c(mesh%num_half_lat))
    allocate(coef%full_dlon(mesh%num_full_lat))
    allocate(coef%half_dlon(mesh%num_half_lat))
    allocate(coef%full_dlat(mesh%num_full_lat))
    allocate(coef%half_dlat(mesh%num_half_lat))

  end subroutine allocate_coef_data

  subroutine allocate_static_data(static)

    type(static_type), intent(out) :: static

    if (.not. allocated(static%ghs))        call parallel_allocate(static%ghs)

  end subroutine allocate_static_data

  subroutine allocate_state_data(state)

    type(state_type), intent(out) :: state

    if (.not. allocated(state%u))           call parallel_allocate(state%u,             half_lon=.true.)
    if (.not. allocated(state%v))           call parallel_allocate(state%v,             half_lat=.true.)
    if (.not. allocated(state%gd))          call parallel_allocate(state%gd)
    if (.not. allocated(state%iap%u))       call parallel_allocate(state%iap%u,         half_lon=.true.)
    if (.not. allocated(state%iap%v))       call parallel_allocate(state%iap%v,         half_lat=.true.)
    if (.not. allocated(state%iap%gd))      call parallel_allocate(state%iap%gd)

  end subroutine allocate_state_data

  subroutine allocate_tend_data(tend)

    type(tend_type), intent(out) :: tend

    if (.not. allocated(tend%u_adv_lon))    call parallel_allocate(tend%u_adv_lon,      half_lon=.true.)
    if (.not. allocated(tend%u_adv_lat))    call parallel_allocate(tend%u_adv_lat,      half_lon=.true.)
    if (.not. allocated(tend%fv))           call parallel_allocate(tend%fv,             half_lon=.true.)
    if (.not. allocated(tend%u_pgf))        call parallel_allocate(tend%u_pgf,          half_lon=.true.)
    if (.not. allocated(tend%du))           call parallel_allocate(tend%du,             half_lon=.true.)
    if (.not. allocated(tend%v_adv_lon))    call parallel_allocate(tend%v_adv_lon,      half_lat=.true.)
    if (.not. allocated(tend%v_adv_lat))    call parallel_allocate(tend%v_adv_lat,      half_lat=.true.)
    if (.not. allocated(tend%fu))           call parallel_allocate(tend%fu,             half_lat=.true.)
    if (.not. allocated(tend%v_pgf))        call parallel_allocate(tend%v_pgf,          half_lat=.true.)
    if (.not. allocated(tend%dv))           call parallel_allocate(tend%dv,             half_lat=.true.)
    if (.not. allocated(tend%mass_div_lon)) call parallel_allocate(tend%mass_div_lon)
    if (.not. allocated(tend%mass_div_lat)) call parallel_allocate(tend%mass_div_lat)
    if (.not. allocated(tend%dgd))          call parallel_allocate(tend%dgd)

  end subroutine allocate_tend_data

  subroutine deallocate_coef_data(coef)

    type(coef_type), intent(inout) :: coef

    if (allocated(coef%full_f))    deallocate(coef%full_f)
    if (allocated(coef%half_f))    deallocate(coef%half_f)
    if (allocated(coef%full_c))    deallocate(coef%full_c)
    if (allocated(coef%half_c))    deallocate(coef%half_c)
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

    if (allocated(state%u))  deallocate(state%u)
    if (allocated(state%v))  deallocate(state%v)
    if (allocated(state%gd)) deallocate(state%gd)

    if (allocated(state%iap%u))  deallocate(state%iap%u)
    if (allocated(state%iap%v))  deallocate(state%iap%v)
    if (allocated(state%iap%gd)) deallocate(state%iap%gd)

  end subroutine deallocate_state_data

  subroutine deallocate_tend_data(tend)

    type(tend_type), intent(inout) :: tend

    if (allocated(tend%u_adv_lon))    deallocate(tend%u_adv_lon)
    if (allocated(tend%u_adv_lat))    deallocate(tend%u_adv_lat)
    if (allocated(tend%v_adv_lon))    deallocate(tend%v_adv_lon)
    if (allocated(tend%v_adv_lat))    deallocate(tend%v_adv_lat)
    if (allocated(tend%fu))           deallocate(tend%fu)
    if (allocated(tend%fv))           deallocate(tend%fv)
    if (allocated(tend%u_pgf))        deallocate(tend%u_pgf)
    if (allocated(tend%v_pgf))        deallocate(tend%v_pgf)
    if (allocated(tend%mass_div_lon)) deallocate(tend%mass_div_lon)
    if (allocated(tend%mass_div_lat)) deallocate(tend%mass_div_lat)
    if (allocated(tend%du))           deallocate(tend%du)
    if (allocated(tend%dv))           deallocate(tend%dv)
    if (allocated(tend%dgd))          deallocate(tend%dgd)

  end subroutine deallocate_tend_data

  subroutine copy_state_data(state1, state2)

    type(state_type), intent(in) :: state1
    type(state_type), intent(inout) :: state2

    state2%u           = state1%u
    state2%v           = state1%v
    state2%gd          = state1%gd
    state2%iap%u       = state1%iap%u
    state2%iap%v       = state1%iap%v
    state2%iap%gd      = state1%iap%gd

  end subroutine copy_state_data

  subroutine copy_tend_data(tend1, tend2)

    type(tend_type), intent(in) :: tend1
    type(tend_type), intent(inout) :: tend2

    tend2%u_adv_lon    = tend1%u_adv_lon
    tend2%u_adv_lat    = tend1%u_adv_lat
    tend2%v_adv_lon    = tend1%v_adv_lon
    tend2%v_adv_lat    = tend1%v_adv_lat
    tend2%fu           = tend1%fu
    tend2%fv           = tend1%fv
    tend2%u_pgf        = tend1%u_pgf
    tend2%v_pgf        = tend1%v_pgf
    tend2%mass_div_lon = tend1%mass_div_lon
    tend2%mass_div_lat = tend1%mass_div_lat
    tend2%du           = tend1%du
    tend2%dv           = tend1%dv
    tend2%dgd          = tend1%dgd

  end subroutine copy_tend_data

  subroutine zero_tend_data(tend)

    type(tend_type), intent(inout) :: tend

    tend%u_adv_lon    = 0.0
    tend%u_adv_lat    = 0.0
    tend%v_adv_lon    = 0.0
    tend%v_adv_lat    = 0.0
    tend%fu           = 0.0
    tend%fv           = 0.0
    tend%u_pgf        = 0.0
    tend%v_pgf        = 0.0
    tend%mass_div_lon = 0.0
    tend%mass_div_lat = 0.0
    tend%du           = 0.0
    tend%dv           = 0.0
    tend%dgd          = 0.0

  end subroutine zero_tend_data

  subroutine add_tend_data(tend1, tend2)

    type(tend_type), intent(in) :: tend1
    type(tend_type), intent(inout) :: tend2

    tend2%u_adv_lon    = tend2%u_adv_lon    + tend1%u_adv_lon
    tend2%u_adv_lat    = tend2%u_adv_lat    + tend1%u_adv_lat
    tend2%v_adv_lon    = tend2%v_adv_lon    + tend1%v_adv_lon
    tend2%v_adv_lat    = tend2%v_adv_lat    + tend1%v_adv_lat
    tend2%fu           = tend2%fu           + tend1%fu
    tend2%fv           = tend2%fv           + tend1%fv
    tend2%u_pgf        = tend2%u_pgf        + tend1%u_pgf
    tend2%v_pgf        = tend2%v_pgf        + tend1%v_pgf
    tend2%mass_div_lon = tend2%mass_div_lon + tend1%mass_div_lon
    tend2%mass_div_lat = tend2%mass_div_lat + tend1%mass_div_lat
    tend2%du           = tend2%du           + tend1%du
    tend2%dv           = tend2%dv           + tend1%dv
    tend2%dgd          = tend2%dgd          + tend1%dgd

  end subroutine add_tend_data

  subroutine sub_tend_data(tend1, tend2)

    type(tend_type), intent(in) :: tend1
    type(tend_type), intent(inout) :: tend2

    tend2%u_adv_lon    = tend2%u_adv_lon    - tend1%u_adv_lon
    tend2%u_adv_lat    = tend2%u_adv_lat    - tend1%u_adv_lat
    tend2%v_adv_lon    = tend2%v_adv_lon    - tend1%v_adv_lon
    tend2%v_adv_lat    = tend2%v_adv_lat    - tend1%v_adv_lat
    tend2%fu           = tend2%fu           - tend1%fu
    tend2%fv           = tend2%fv           - tend1%fv
    tend2%u_pgf        = tend2%u_pgf        - tend1%u_pgf
    tend2%v_pgf        = tend2%v_pgf        - tend1%v_pgf
    tend2%mass_div_lon = tend2%mass_div_lon - tend1%mass_div_lon
    tend2%mass_div_lat = tend2%mass_div_lat - tend1%mass_div_lat
    tend2%du           = tend2%du           - tend1%du
    tend2%dv           = tend2%dv           - tend1%dv
    tend2%dgd          = tend2%dgd          - tend1%dgd

  end subroutine sub_tend_data

  subroutine scale_tend_data(scale, tend)

    real, intent(in) :: scale
    type(tend_type), intent(inout) :: tend

    tend%u_adv_lon    = tend%u_adv_lon    * scale
    tend%u_adv_lat    = tend%u_adv_lat    * scale
    tend%v_adv_lon    = tend%v_adv_lon    * scale
    tend%v_adv_lat    = tend%v_adv_lat    * scale
    tend%fu           = tend%fu           * scale
    tend%fv           = tend%fv           * scale
    tend%u_pgf        = tend%u_pgf        * scale
    tend%v_pgf        = tend%v_pgf        * scale
    tend%mass_div_lon = tend%mass_div_lon * scale
    tend%mass_div_lat = tend%mass_div_lat * scale
    tend%du           = tend%du           * scale
    tend%dv           = tend%dv           * scale
    tend%dgd          = tend%dgd          * scale

  end subroutine scale_tend_data

  real function inner_product_tend_tend(tend1, tend2) result(res)

    type(tend_type), intent(in) :: tend1
    type(tend_type), intent(in) :: tend2

    integer i, j

    res = 0.0
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        res = res + tend1%du(i,j) * tend2%du(i,j) * mesh%full_cos_lat(j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%dv(i,j) * tend2%dv(i,j) * mesh%half_cos_lat(j)
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend1%dgd(i,j) * tend2%dgd(i,j) * mesh%full_cos_lat(j)
      end do
    end do

  end function inner_product_tend_tend

  real function inner_product_tend_state(tend, state) result(res)

    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: state

    integer i, j

    res = 0.0
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        res = res + tend%du(i,j) * state%iap%u(i,j) * mesh%full_cos_lat(j)
      end do
    end do
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend%dv(i,j) * state%iap%v(i,j) * mesh%half_cos_lat(j)
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + tend%dgd(i,j) * state%gd(i,j) * mesh%full_cos_lat(j)
      end do
    end do

  end function inner_product_tend_state

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
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%iap%u(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i+1,j)) * state%u(i,j)
      end do
    end do
    call parallel_fill_halo(state%iap%u(:,:), all_halo=.true.)

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%v(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%v(i,j)
      end do
    end do
    call parallel_fill_halo(state%iap%v(:,:), all_halo=.true.)

  end subroutine iap_transform

end module types_mod
