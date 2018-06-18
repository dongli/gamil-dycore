module parallel_mod

#ifndef NO_MPI
  use mpi
#endif
  use log_mod
  use mesh_mod
  use params_mod

  implicit none

  private

  public parallel
  public parallel_init
  public parallel_final
  public parallel_allocate
  public parallel_fill_halo
  public parallel_overlay_inner_halo
  public parallel_zonal_sum

  type parallel_params_type
    integer full_lon_start_idx
    integer full_lon_end_idx
    integer half_lon_start_idx
    integer half_lon_end_idx
    integer full_lat_start_idx
    integer full_lat_end_idx
    integer half_lat_start_idx
    integer half_lat_end_idx
    integer full_lat_start_idx_no_pole
    integer full_lat_end_idx_no_pole
    integer half_lat_start_idx_no_pole
    integer half_lat_end_idx_no_pole
    integer full_lat_south_pole_idx
    integer full_lat_north_pole_idx
    integer half_lat_south_pole_idx
    integer half_lat_north_pole_idx
    integer lon_halo_width
    integer lat_halo_width
    integer lon_halo_width_for_reduce
    integer lat_halo_width_for_reduce
    integer full_lon_lb
    integer half_lon_lb
    integer full_lon_ub
    integer half_lon_ub
    integer full_lat_lb
    integer half_lat_lb
    integer full_lat_ub
    integer half_lat_ub
    integer full_lon_lb_for_reduce
    integer half_lon_lb_for_reduce
    integer full_lon_ub_for_reduce
    integer half_lon_ub_for_reduce
    integer full_lat_lb_for_reduce
    integer half_lat_lb_for_reduce
    integer full_lat_ub_for_reduce
    integer half_lat_ub_for_reduce
    logical has_south_pole
    logical has_north_pole
  end type parallel_params_type

  type(parallel_params_type) parallel

  interface parallel_allocate
    module procedure parallel_allocate_1
    module procedure parallel_allocate_1_add_dim_2
    module procedure parallel_allocate_1_i
    module procedure parallel_allocate_2
    module procedure parallel_allocate_2_add_dim_1
    module procedure parallel_allocate_2_add_dim_2
  end interface parallel_allocate

  interface parallel_fill_halo
    module procedure parallel_fill_halo_1
    module procedure parallel_fill_halo_2
  end interface parallel_fill_halo

  interface parallel_overlay_inner_halo
    module procedure parallel_overlay_inner_halo_1
  end interface parallel_overlay_inner_halo

contains

  subroutine parallel_init()

    parallel%full_lon_start_idx = 1
    parallel%full_lon_end_idx = mesh%num_full_lon
    parallel%half_lon_start_idx = 1
    parallel%half_lon_end_idx = mesh%num_half_lon
    parallel%full_lat_start_idx = 1
    parallel%full_lat_end_idx = mesh%num_full_lat
    parallel%half_lat_start_idx = 1
    parallel%half_lat_end_idx = mesh%num_half_lat

    parallel%full_lat_start_idx_no_pole = parallel%full_lat_start_idx + 1
    parallel%full_lat_end_idx_no_pole = parallel%full_lat_end_idx - 1
    parallel%half_lat_start_idx_no_pole = parallel%half_lat_start_idx + 1
    parallel%half_lat_end_idx_no_pole = parallel%half_lat_end_idx - 1

    parallel%full_lat_south_pole_idx = 1
    parallel%full_lat_north_pole_idx = mesh%num_full_lat
    parallel%half_lat_south_pole_idx = 1
    parallel%half_lat_north_pole_idx = mesh%num_half_lat

    parallel%has_south_pole = parallel%full_lat_start_idx == parallel%full_lat_south_pole_idx
    parallel%has_north_pole = parallel%full_lat_end_idx == parallel%full_lat_north_pole_idx

    parallel%lon_halo_width = 1
    parallel%lat_halo_width = 1

    parallel%full_lon_lb = parallel%full_lon_start_idx - parallel%lon_halo_width
    parallel%full_lon_ub = parallel%full_lon_end_idx + parallel%lon_halo_width
    parallel%half_lon_lb = parallel%half_lon_start_idx - parallel%lon_halo_width
    parallel%half_lon_ub = parallel%half_lon_end_idx + parallel%lon_halo_width
    parallel%full_lat_lb = parallel%full_lat_start_idx - parallel%lat_halo_width
    parallel%full_lat_ub = parallel%full_lat_end_idx + parallel%lat_halo_width
    parallel%half_lat_lb = parallel%half_lat_start_idx - parallel%lat_halo_width
    parallel%half_lat_ub = parallel%half_lat_end_idx + parallel%lat_halo_width

    parallel%lon_halo_width_for_reduce = max(1, maxval(zonal_reduce_factors))
    parallel%lat_halo_width_for_reduce = 1

    parallel%full_lon_lb_for_reduce = parallel%full_lon_start_idx - parallel%lon_halo_width_for_reduce
    parallel%full_lon_ub_for_reduce = parallel%full_lon_end_idx + parallel%lon_halo_width_for_reduce
    parallel%half_lon_lb_for_reduce = parallel%half_lon_start_idx - parallel%lon_halo_width_for_reduce
    parallel%half_lon_ub_for_reduce = parallel%half_lon_end_idx + parallel%lon_halo_width_for_reduce
    parallel%full_lat_lb_for_reduce = parallel%full_lat_start_idx - parallel%lat_halo_width_for_reduce
    parallel%full_lat_ub_for_reduce = parallel%full_lat_end_idx + parallel%lat_halo_width_for_reduce
    parallel%half_lat_lb_for_reduce = parallel%half_lat_start_idx - parallel%lat_halo_width_for_reduce
    parallel%half_lat_ub_for_reduce = parallel%half_lat_end_idx + parallel%lat_halo_width_for_reduce

    call log_notice('Parallel module is initialized.')

  end subroutine parallel_init

  subroutine parallel_final()

    call log_notice('Parallel module is finalized.')

  end subroutine parallel_final

  subroutine parallel_allocate_1(field, full_lon, half_lon, full_lat, half_lat, extended_halo)

    real, intent(out), allocatable :: field(:)
    logical, intent(in), optional :: full_lon
    logical, intent(in), optional :: half_lon
    logical, intent(in), optional :: full_lat
    logical, intent(in), optional :: half_lat
    logical, intent(in), optional :: extended_halo

    if (present(full_lon) .and. full_lon) then
      if (present(extended_halo) .and. extended_halo) then
        allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce))
      else
        allocate(field(parallel%full_lon_lb:parallel%full_lon_ub))
      end if
    else if (present(half_lon) .and. half_lon) then
      if (present(extended_halo) .and. extended_halo) then
        allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce))
      else
        allocate(field(parallel%half_lon_lb:parallel%half_lon_ub))
      end if
    else if (present(full_lat) .and. full_lat) then
      if (present(extended_halo) .and. extended_halo) then
        allocate(field(parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
      else
        allocate(field(parallel%full_lat_lb:parallel%full_lat_ub))
      end if
    else if (present(half_lat) .and. half_lat) then
      if (present(extended_halo) .and. extended_halo) then
        allocate(field(parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
      else
        allocate(field(parallel%half_lat_lb:parallel%half_lat_ub))
      end if
    end if

    field(:) = 0.0

  end subroutine parallel_allocate_1

  subroutine parallel_allocate_1_add_dim_2(field, dim, size, full_lon, half_lon, full_lat, half_lat, extended_halo)

    real, intent(out), allocatable ::  field(:,:,:)
    integer, intent(in) :: dim(2)
    integer, intent(in) :: size(2)
    logical, intent(in), optional :: full_lon
    logical, intent(in), optional :: half_lon
    logical, intent(in), optional :: full_lat
    logical, intent(in), optional :: half_lat
    logical, intent(in), optional :: extended_halo

    if (dim(1) == 1 .and. dim(2) == 2) then
      if (present(half_lon) .and. half_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),size(2),parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce))
        else
          allocate(field(size(1),size(2),parallel%half_lon_lb:parallel%half_lon_ub))
        end if
      else if (present(full_lon) .and. full_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),size(2),parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce))
        else
          allocate(field(size(1),size(2),parallel%full_lon_lb:parallel%full_lon_ub))
        end if
      else if (present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),size(2),parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
        else
          allocate(field(size(1),size(2),parallel%half_lat_lb:parallel%half_lat_ub))
        end if
      else if (present(full_lat) .and. full_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),size(2),parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
        else
          allocate(field(size(1),size(2),parallel%full_lat_lb:parallel%full_lat_ub))
        end if
      end if
    else if (dim(1) == 1 .and. dim(2) == 3) then
      if (present(half_lon) .and. half_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,size(2)))
        else
          allocate(field(size(1),parallel%half_lon_lb:parallel%half_lon_ub,size(2)))
        end if
      else if (present(full_lon) .and. full_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,size(2)))
        else
          allocate(field(size(1),parallel%full_lon_lb:parallel%full_lon_ub,size(2)))
        end if
      else if (present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce,size(2)))
        else
          allocate(field(size(1),parallel%half_lat_lb:parallel%half_lat_ub,size(2)))
        end if
      else if (present(full_lat) .and. full_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce,size(2)))
        else
          allocate(field(size(1),parallel%full_lat_lb:parallel%full_lat_ub,size(2)))
        end if
      end if
    else if (dim(1) == 2 .and. dim(2) == 3) then
      if (present(half_lon) .and. half_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,size(1),size(2)))
        else
          allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,size(1),size(2)))
        end if
      else if (present(full_lon) .and. full_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,size(1),size(2)))
        else
          allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,size(1),size(2)))
        end if
      else if (present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce,size(1),size(2)))
        else
          allocate(field(parallel%half_lat_lb:parallel%half_lat_ub,size(1),size(2)))
        end if
      else if (present(full_lat) .and. full_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce,size(1),size(2)))
        else
          allocate(field(parallel%full_lat_lb:parallel%full_lat_ub,size(1),size(2)))
        end if
      end if
    else
      call log_error('Inconsistent dimension arguments!')
    end if

    ! Initialize field with zeros.
    field(:,:,:) = 0.0

  end subroutine parallel_allocate_1_add_dim_2

  subroutine parallel_allocate_1_i(field, full_lon, half_lon, full_lat, half_lat, extended_halo)

    integer, intent(out), allocatable :: field(:)
    logical, intent(in), optional :: full_lon
    logical, intent(in), optional :: half_lon
    logical, intent(in), optional :: full_lat
    logical, intent(in), optional :: half_lat
    logical, intent(in), optional :: extended_halo

    if (present(full_lon) .and. full_lon) then
      if (present(extended_halo) .and. extended_halo) then
        allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce))
      else
        allocate(field(parallel%full_lon_lb:parallel%full_lon_ub))
      end if
    else if (present(half_lon) .and. half_lon) then
      if (present(extended_halo) .and. extended_halo) then
        allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce))
      else
        allocate(field(parallel%half_lon_lb:parallel%half_lon_ub))
      end if
    else if (present(full_lat) .and. full_lat) then
      if (present(extended_halo) .and. extended_halo) then
        allocate(field(parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
      else
        allocate(field(parallel%full_lat_lb:parallel%full_lat_ub))
      end if
    else if (present(half_lat) .and. half_lat) then
      if (present(extended_halo) .and. extended_halo) then
        allocate(field(parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
      else
        allocate(field(parallel%half_lat_lb:parallel%half_lat_ub))
      end if
    end if

    field(:) = 0

  end subroutine parallel_allocate_1_i

  subroutine parallel_allocate_2(field, half_lon, half_lat, extended_halo)

    real, intent(out), allocatable ::  field(:,:)
    logical, intent(in), optional :: half_lon
    logical, intent(in), optional :: half_lat
    logical, intent(in), optional :: extended_halo

    if (present(half_lon) .and. half_lon .and. present(half_lat) .and. half_lat) then
      if (present(extended_halo) .and. extended_halo) then
        allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
      else
        allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub))
      end if
    else if (present(half_lon) .and. half_lon) then
      if (present(extended_halo) .and. extended_halo) then
        allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
      else
        allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub))
      end if
    else if (present(half_lat) .and. half_lat) then
      if (present(extended_halo) .and. extended_halo) then
        allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
      else
        allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub))
      end if
    else
      if (present(extended_halo) .and. extended_halo) then
        allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
      else
        allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub))
      end if
    end if

    ! Initialize field with zeros.
    field(:,:) = 0.0

  end subroutine parallel_allocate_2

  subroutine parallel_allocate_2_add_dim_1(field, dim, size, half_lon, half_lat, extended_halo)

    real, intent(out), allocatable ::  field(:,:,:)
    integer, intent(in) :: dim
    integer, intent(in) :: size
    logical, intent(in), optional :: half_lon
    logical, intent(in), optional :: half_lat
    logical, intent(in), optional :: extended_halo

    select case (dim)
    case (1)
      if (present(half_lon) .and. half_lon .and. present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size,parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
        else
          allocate(field(size,parallel%half_lon_lb:parallel%half_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub))
        end if
      else if (present(half_lon) .and. half_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size,parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
        else
          allocate(field(size,parallel%half_lon_lb:parallel%half_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub))
        end if
      else if (present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size,parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
        else
          allocate(field(size,parallel%full_lon_lb:parallel%full_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub))
        end if
      else
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size,parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
        else
          allocate(field(size,parallel%full_lon_lb:parallel%full_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub))
        end if
      end if
    case (2)
      if (present(half_lon) .and. half_lon .and. present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,size,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
        else
          allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,size,parallel%half_lat_lb:parallel%half_lat_ub))
        end if
      else if (present(half_lon) .and. half_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,size,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
        else
          allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,size,parallel%full_lat_lb:parallel%full_lat_ub))
        end if
      else if (present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,size,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
        else
          allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,size,parallel%half_lat_lb:parallel%half_lat_ub))
        end if
      else
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,size,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
        else
          allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,size,parallel%full_lat_lb:parallel%full_lat_ub))
        end if
      end if
    case (3)
      if (present(half_lon) .and. half_lon .and. present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce,size))
        else
          allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub,size))
        end if
      else if (present(half_lon) .and. half_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce,size))
        else
          allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub,size))
        end if
      else if (present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce,size))
        else
          allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub,size))
        end if
      else
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce,size))
        else
          allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub,size))
        end if
      end if
    end select

    ! Initialize field with zeros.
    field(:,:,:) = 0.0

  end subroutine parallel_allocate_2_add_dim_1

  subroutine parallel_allocate_2_add_dim_2(field, dim, size, half_lon, half_lat, extended_halo)

    real, intent(out), allocatable ::  field(:,:,:,:)
    integer, intent(in) :: dim(2)
    integer, intent(in) :: size(2)
    logical, intent(in), optional :: half_lon
    logical, intent(in), optional :: half_lat
    logical, intent(in), optional :: extended_halo

    if (dim(1) == 1 .and. dim(2) == 2) then
      if (present(half_lon) .and. half_lon .and. present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),size(2),parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
        else
          allocate(field(size(1),size(2),parallel%half_lon_lb:parallel%half_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub))
        end if
      else if (present(half_lon) .and. half_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),size(2),parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
        else
          allocate(field(size(1),size(2),parallel%half_lon_lb:parallel%half_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub))
        end if
      else if (present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),size(2),parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
        else
          allocate(field(size(1),size(2),parallel%full_lon_lb:parallel%full_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub))
        end if
      else
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),size(2),parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
        else
          allocate(field(size(1),size(2),parallel%full_lon_lb:parallel%full_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub))
        end if
      end if
    else if (dim(1) == 1 .and. dim(2) == 3) then
      if (present(half_lon) .and. half_lon .and. present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,size(2),parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
        else
          allocate(field(size(1),parallel%half_lon_lb:parallel%half_lon_ub,size(2),parallel%half_lat_lb:parallel%half_lat_ub))
        end if
      else if (present(half_lon) .and. half_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,size(2),parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
        else
          allocate(field(size(1),parallel%half_lon_lb:parallel%half_lon_ub,size(2),parallel%full_lat_lb:parallel%full_lat_ub))
        end if
      else if (present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,size(2),parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
        else
          allocate(field(size(1),parallel%full_lon_lb:parallel%full_lon_ub,size(2),parallel%half_lat_lb:parallel%half_lat_ub))
        end if
      else
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,size(2),parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
        else
          allocate(field(size(1),parallel%full_lon_lb:parallel%full_lon_ub,size(2),parallel%full_lat_lb:parallel%full_lat_ub))
        end if
      end if
    else if (dim(1) == 1 .and. dim(2) == 4) then
      if (present(half_lon) .and. half_lon .and. present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce,size(2)))
        else
          allocate(field(size(1),parallel%half_lon_lb:parallel%half_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub,size(2)))
        end if
      else if (present(half_lon) .and. half_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce,size(2)))
        else
          allocate(field(size(1),parallel%half_lon_lb:parallel%half_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub,size(2)))
        end if
      else if (present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce,size(2)))
        else
          allocate(field(size(1),parallel%full_lon_lb:parallel%full_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub,size(2)))
        end if
      else
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(size(1),parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce,size(2)))
        else
          allocate(field(size(1),parallel%full_lon_lb:parallel%full_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub,size(2)))
        end if
      end if
    else if (dim(1) == 2 .and. dim(2) == 3) then
      if (present(half_lon) .and. half_lon .and. present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,size(1),size(2),parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
        else
          allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,size(1),size(2),parallel%half_lat_lb:parallel%half_lat_ub))
        end if
      else if (present(half_lon) .and. half_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,size(1),size(2),parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
        else
          allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,size(1),size(2),parallel%full_lat_lb:parallel%full_lat_ub))
        end if
      else if (present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,size(1),size(2),parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce))
        else
          allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,size(1),size(2),parallel%half_lat_lb:parallel%half_lat_ub))
        end if
      else
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,size(1),size(2),parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce))
        else
          allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,size(1),size(2),parallel%full_lat_lb:parallel%full_lat_ub))
        end if
      end if
    else if (dim(1) == 2 .and. dim(2) == 4) then
      if (present(half_lon) .and. half_lon .and. present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,size(1),parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce,size(2)))
        else
          allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,size(1),parallel%half_lat_lb:parallel%half_lat_ub,size(2)))
        end if
      else if (present(half_lon) .and. half_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,size(1),parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce,size(2)))
        else
          allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,size(1),parallel%full_lat_lb:parallel%full_lat_ub,size(2)))
        end if
      else if (present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,size(1),parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce,size(2)))
        else
          allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,size(1),parallel%half_lat_lb:parallel%half_lat_ub,size(2)))
        end if
      else
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,size(1),parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce,size(2)))
        else
          allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,size(1),parallel%full_lat_lb:parallel%full_lat_ub,size(2)))
        end if
      end if
    else if (dim(1) == 3 .and. dim(2) == 4) then
      if (present(half_lon) .and. half_lon .and. present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce,size(1),size(2)))
        else
          allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub,size(1),size(2)))
        end if
      else if (present(half_lon) .and. half_lon) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%half_lon_lb_for_reduce:parallel%half_lon_ub_for_reduce,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce,size(1),size(2)))
        else
          allocate(field(parallel%half_lon_lb:parallel%half_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub,size(1),size(2)))
        end if
      else if (present(half_lat) .and. half_lat) then
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,parallel%half_lat_lb_for_reduce:parallel%half_lat_ub_for_reduce,size(1),size(2)))
        else
          allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,parallel%half_lat_lb:parallel%half_lat_ub,size(1),size(2)))
        end if
      else
        if (present(extended_halo) .and. extended_halo) then
          allocate(field(parallel%full_lon_lb_for_reduce:parallel%full_lon_ub_for_reduce,parallel%full_lat_lb_for_reduce:parallel%full_lat_ub_for_reduce,size(1),size(2)))
        else
          allocate(field(parallel%full_lon_lb:parallel%full_lon_ub,parallel%full_lat_lb:parallel%full_lat_ub,size(1),size(2)))
        end if
      end if
    else
      call log_error('Inconsistent dimension arguments!')
    end if

    ! Initialize field with zeros.
    field(:,:,:,:) = 0.0

  end subroutine parallel_allocate_2_add_dim_2

  subroutine parallel_fill_halo_1(field, all_halo, left_halo, right_halo, top_halo, bottom_halo)

    real, intent(inout) :: field(:,:)
    logical, intent(in), optional :: all_halo
    logical, intent(in), optional :: left_halo
    logical, intent(in), optional :: right_halo
    logical, intent(in), optional :: top_halo
    logical, intent(in), optional :: bottom_halo

    integer i, j, m, n

    if ((present(all_halo) .and. all_halo) .or. (present(left_halo) .and. left_halo)) then
      m = lbound(field, 1) - 1
      n = ubound(field, 1) - 2 * parallel%lon_halo_width_for_reduce
      do j = lbound(field, 2), ubound(field, 2)
        do i = 1, parallel%lon_halo_width_for_reduce
          field(m+i,j) = field(n+i,j)
        end do
      end do
    end if

    ! |             |                             |              |              |
    ! lb            lb + w                        ub - 2w        ub - w         ub
    if ((present(all_halo) .and. all_halo) .or. (present(right_halo) .and. right_halo)) then
      m = ubound(field, 1) - parallel%lon_halo_width_for_reduce
      n = lbound(field, 1) + parallel%lon_halo_width_for_reduce - 1
      do j = lbound(field, 2), ubound(field, 2)
        do i = 1, parallel%lon_halo_width_for_reduce
          field(m+i,j) = field(n+i,j)
        end do
      end do
    end if

  end subroutine parallel_fill_halo_1

  subroutine parallel_fill_halo_2(field, left_halo, right_halo)

    real, intent(inout) :: field(:)
    logical, intent(in), optional :: left_halo
    logical, intent(in), optional :: right_halo

    integer i, m, n

    if (present(left_halo) .and. left_halo) then
      m = lbound(field, 1) - 1
      n = ubound(field, 1) - 2 * parallel%lon_halo_width_for_reduce
      do i = 1, parallel%lon_halo_width_for_reduce
        field(m+i) = field(n+i)
      end do
    end if

    ! |             |                             |              |              |
    ! lb            lb + w                        ub - 2w        ub - w         ub
    if (present(right_halo) .and. right_halo) then
      m = ubound(field, 1) - parallel%lon_halo_width_for_reduce
      n = lbound(field, 1) + parallel%lon_halo_width_for_reduce - 1
      do i = 1, parallel%lon_halo_width_for_reduce
        field(m+i) = field(n+i)
      end do
    end if

  end subroutine parallel_fill_halo_2

  subroutine parallel_overlay_inner_halo_1(array, left_halo, right_halo)

    real, intent(inout) :: array(:)
    logical, intent(in), optional :: left_halo
    logical, intent(in), optional :: right_halo

    integer i, m, n

    if (present(left_halo) .and. left_halo) then
      m = lbound(array, 1) + parallel%lon_halo_width_for_reduce - 1
      n = ubound(array, 1) - parallel%lon_halo_width_for_reduce
      do i = 1, parallel%lon_halo_width_for_reduce - 1
        array(m+i) = array(m+i) + array(n+i)
      end do
    end if

    if (present(right_halo) .and. right_halo) then
      m = ubound(array, 1) - 2 * parallel%lon_halo_width_for_reduce
      n = lbound(array, 1)
      do i = 1, parallel%lon_halo_width_for_reduce - 1
        array(m+i) = array(m+i) + array(n+i)
      end do
    end if

  end subroutine parallel_overlay_inner_halo_1

  subroutine parallel_zonal_sum(send_buf, recv_buf)

    real, intent(in) :: send_buf
    real, intent(out) :: recv_buf

    recv_buf = send_buf

  end subroutine parallel_zonal_sum

end module parallel_mod
