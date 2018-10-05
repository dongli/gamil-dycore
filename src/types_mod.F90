module types_mod

  use mesh_mod
  use parallel_mod

  implicit none

  private

  public allocate_data
  public deallocate_data
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
    real, allocatable :: u_c(:,:)
    real, allocatable :: v_c(:,:)
    real, allocatable :: u_d(:,:)
    real, allocatable :: v_d(:,:)
    real, allocatable :: gd(:,:)
  end type iap_type

  type state_type
    real, allocatable :: u_c(:,:)
    real, allocatable :: v_c(:,:)
    real, allocatable :: u_d(:,:)
    real, allocatable :: v_d(:,:)
    real, allocatable :: gd(:,:) ! Geopotential depth
    type(iap_type) iap
  end type state_type

  type static_type
    real, allocatable :: ghs(:,:) ! Surface geopotential
  end type static_type

  type tend_type
    real, allocatable :: u_adv_lon_c(:,:)
    real, allocatable :: u_adv_lat_c(:,:)
    real, allocatable :: v_adv_lon_c(:,:)
    real, allocatable :: v_adv_lat_c(:,:)
    real, allocatable :: fu_c(:,:)
    real, allocatable :: fv_c(:,:)
    real, allocatable :: u_pgf_c(:,:)
    real, allocatable :: v_pgf_c(:,:)
    real, allocatable :: du_c(:,:)
    real, allocatable :: dv_c(:,:)
    real, allocatable :: u_adv_lon_d(:,:)
    real, allocatable :: u_adv_lat_d(:,:)
    real, allocatable :: v_adv_lon_d(:,:)
    real, allocatable :: v_adv_lat_d(:,:)
    real, allocatable :: fu_d(:,:)
    real, allocatable :: fv_d(:,:)
    real, allocatable :: u_pgf_d(:,:)
    real, allocatable :: v_pgf_d(:,:)
    real, allocatable :: du_d(:,:)
    real, allocatable :: dv_d(:,:)
    real, allocatable :: mass_div_lon(:,:)
    real, allocatable :: mass_div_lat(:,:)
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

    if (.not. allocated(static%ghs)) call parallel_allocate(static%ghs)

  end subroutine allocate_static_data

  subroutine allocate_state_data(state)

    type(state_type), intent(out) :: state

    if (.not. allocated(state%u_c))     call parallel_allocate(state%u_c,           half_lon=.true.)
    if (.not. allocated(state%v_c))     call parallel_allocate(state%v_c,           half_lat=.true.)
    if (.not. allocated(state%u_d))     call parallel_allocate(state%u_d,           half_lat=.true.)
    if (.not. allocated(state%v_d))     call parallel_allocate(state%v_d,           half_lon=.true.)
    if (.not. allocated(state%gd))      call parallel_allocate(state%gd)
    if (.not. allocated(state%iap%u_c)) call parallel_allocate(state%iap%u_c,       half_lon=.true.)
    if (.not. allocated(state%iap%v_c)) call parallel_allocate(state%iap%v_c,       half_lat=.true.)
    if (.not. allocated(state%iap%u_d)) call parallel_allocate(state%iap%u_d,       half_lat=.true.)
    if (.not. allocated(state%iap%v_d)) call parallel_allocate(state%iap%v_d,       half_lon=.true.)
    if (.not. allocated(state%iap%gd))  call parallel_allocate(state%iap%gd)

  end subroutine allocate_state_data

  subroutine allocate_tend_data(tend)

    type(tend_type), intent(out) :: tend

    if (.not. allocated(tend%u_adv_lon_c))  call parallel_allocate(tend%u_adv_lon_c,    half_lon=.true.)
    if (.not. allocated(tend%u_adv_lat_c))  call parallel_allocate(tend%u_adv_lat_c,    half_lon=.true.)
    if (.not. allocated(tend%v_adv_lon_c))  call parallel_allocate(tend%v_adv_lon_c,    half_lat=.true.)
    if (.not. allocated(tend%v_adv_lat_c))  call parallel_allocate(tend%v_adv_lat_c,    half_lat=.true.)
    if (.not. allocated(tend%fu_c))         call parallel_allocate(tend%fu_c,           half_lat=.true.)
    if (.not. allocated(tend%fv_c))         call parallel_allocate(tend%fv_c,           half_lon=.true.)
    if (.not. allocated(tend%u_pgf_c))      call parallel_allocate(tend%u_pgf_c,        half_lon=.true.)
    if (.not. allocated(tend%v_pgf_c))      call parallel_allocate(tend%v_pgf_c,        half_lat=.true.)
    if (.not. allocated(tend%du_c))         call parallel_allocate(tend%du_c,           half_lon=.true.)
    if (.not. allocated(tend%dv_c))         call parallel_allocate(tend%dv_c,           half_lat=.true.)
    if (.not. allocated(tend%u_adv_lon_d))  call parallel_allocate(tend%u_adv_lon_d,    half_lat=.true.)
    if (.not. allocated(tend%u_adv_lat_d))  call parallel_allocate(tend%u_adv_lat_d,    half_lat=.true.)
    if (.not. allocated(tend%v_adv_lon_d))  call parallel_allocate(tend%v_adv_lon_d,    half_lon=.true.)
    if (.not. allocated(tend%v_adv_lat_d))  call parallel_allocate(tend%v_adv_lat_d,    half_lon=.true.)
    if (.not. allocated(tend%fu_d))         call parallel_allocate(tend%fu_d,           half_lon=.true.)
    if (.not. allocated(tend%fv_d))         call parallel_allocate(tend%fv_d,           half_lat=.true.)
    if (.not. allocated(tend%u_pgf_d))      call parallel_allocate(tend%u_pgf_d,        half_lat=.true.)
    if (.not. allocated(tend%v_pgf_d))      call parallel_allocate(tend%v_pgf_d,        half_lon=.true.)
    if (.not. allocated(tend%du_d))         call parallel_allocate(tend%du_d,           half_lat=.true.)
    if (.not. allocated(tend%dv_d))         call parallel_allocate(tend%dv_d,           half_lon=.true.)
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

    if (allocated(state%u_c))     deallocate(state%u_c)
    if (allocated(state%v_c))     deallocate(state%v_c)
    if (allocated(state%u_d))     deallocate(state%u_d)
    if (allocated(state%v_d))     deallocate(state%v_d)
    if (allocated(state%gd))      deallocate(state%gd)
    if (allocated(state%iap%u_d)) deallocate(state%iap%u_d)
    if (allocated(state%iap%v_d)) deallocate(state%iap%v_d)
    if (allocated(state%iap%gd))  deallocate(state%iap%gd)

  end subroutine deallocate_state_data

  subroutine deallocate_tend_data(tend)

    type(tend_type), intent(inout) :: tend

    if (allocated(tend%u_adv_lon_c))  deallocate(tend%u_adv_lon_c)
    if (allocated(tend%u_adv_lat_c))  deallocate(tend%u_adv_lat_c)
    if (allocated(tend%v_adv_lon_c))  deallocate(tend%v_adv_lon_c)
    if (allocated(tend%v_adv_lat_c))  deallocate(tend%v_adv_lat_c)
    if (allocated(tend%fu_c))         deallocate(tend%fu_c)
    if (allocated(tend%fv_c))         deallocate(tend%fv_c)
    if (allocated(tend%u_pgf_c))      deallocate(tend%u_pgf_c)
    if (allocated(tend%v_pgf_c))      deallocate(tend%v_pgf_c)
    if (allocated(tend%du_c))         deallocate(tend%du_c)
    if (allocated(tend%dv_c))         deallocate(tend%dv_c)
    if (allocated(tend%u_adv_lon_d))  deallocate(tend%u_adv_lon_d)
    if (allocated(tend%u_adv_lat_d))  deallocate(tend%u_adv_lat_d)
    if (allocated(tend%v_adv_lon_d))  deallocate(tend%v_adv_lon_d)
    if (allocated(tend%v_adv_lat_d))  deallocate(tend%v_adv_lat_d)
    if (allocated(tend%fu_d))         deallocate(tend%fu_d)
    if (allocated(tend%fv_d))         deallocate(tend%fv_d)
    if (allocated(tend%u_pgf_d))      deallocate(tend%u_pgf_d)
    if (allocated(tend%v_pgf_d))      deallocate(tend%v_pgf_d)
    if (allocated(tend%du_d))         deallocate(tend%du_d)
    if (allocated(tend%dv_d))         deallocate(tend%dv_d)
    if (allocated(tend%mass_div_lon)) deallocate(tend%mass_div_lon)
    if (allocated(tend%mass_div_lat)) deallocate(tend%mass_div_lat)
    if (allocated(tend%dgd))          deallocate(tend%dgd)

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
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%iap%u_c(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i+1,j)) * state%u_c(i,j)
        state%iap%v_d(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i+1,j)) * state%v_d(i,j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%iap%u_d(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%u_d(i,j)
        state%iap%v_c(i,j) = 0.5 * (state%iap%gd(i,j) + state%iap%gd(i,j+1)) * state%v_c(i,j)
      end do
    end do

    call parallel_fill_halo(state%iap%u_c(:,:), all_halo=.true.)
    call parallel_fill_halo(state%iap%v_c(:,:), all_halo=.true.)
    call parallel_fill_halo(state%iap%u_d(:,:), all_halo=.true.)
    call parallel_fill_halo(state%iap%v_d(:,:), all_halo=.true.)

  end subroutine iap_transform

end module types_mod
