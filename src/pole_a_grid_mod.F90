module pole_a_grid_mod

  use mesh_mod
  use parallel_mod
  use log_mod
  use params_mod
  use types_mod

  implicit none

  private

  public pole_a_grid_init
  public pole_a_grid_init_state
  public full_lat_start_idx_a_grid
  public full_lat_end_idx_a_grid
  public full_lat_start_idx_c_grid
  public full_lat_end_idx_c_grid
  public half_lat_start_idx_c_grid
  public half_lat_end_idx_c_grid
  public full_lat_start_idx_no_pole_a_grid
  public full_lat_end_idx_no_pole_a_grid
  public full_lat_start_idx_no_pole_c_grid
  public full_lat_end_idx_no_pole_c_grid
  public half_lat_start_idx_no_pole_c_grid
  public half_lat_end_idx_no_pole_c_grid

  integer full_lat_start_idx_no_pole_a_grid(2), full_lat_end_idx_no_pole_a_grid(2)
  integer full_lat_start_idx_no_pole_c_grid,    full_lat_end_idx_no_pole_c_grid
  integer half_lat_start_idx_no_pole_c_grid,    half_lat_end_idx_no_pole_c_grid

  integer full_lat_start_idx_a_grid(2), full_lat_end_idx_a_grid(2)
  integer full_lat_start_idx_c_grid,    full_lat_end_idx_c_grid
  integer half_lat_start_idx_c_grid,    half_lat_end_idx_c_grid

contains

  subroutine pole_a_grid_init()

    integer j

    full_lat_start_idx_a_grid(:) = 0
    full_lat_end_idx_a_grid(:)   = -1
    full_lat_start_idx_c_grid = parallel%full_lat_start_idx
    full_lat_end_idx_c_grid   = parallel%full_lat_end_idx
    half_lat_start_idx_c_grid = parallel%half_lat_start_idx
    half_lat_end_idx_c_grid   = parallel%half_lat_end_idx

    full_lat_start_idx_no_pole_a_grid = 0
    full_lat_end_idx_no_pole_a_grid   = -1
    full_lat_start_idx_no_pole_c_grid = parallel%full_lat_start_idx_no_pole
    full_lat_end_idx_no_pole_c_grid   = parallel%full_lat_end_idx_no_pole
    half_lat_start_idx_no_pole_c_grid = parallel%half_lat_start_idx_no_pole
    half_lat_end_idx_no_pole_c_grid   = parallel%half_lat_end_idx_no_pole

    ! Handle unreasonable input parameters.
    if (use_pole_a_grid .and. pole_a_grid_tag(1) == 0) then
      use_pole_a_grid = .false.
    end if
    if (use_pole_a_grid .and. pole_a_grid_tag(1) == 2) then
      ! Use A grid globally.
      use_all_a_grid = .true.
      full_lat_start_idx_a_grid(:) = [parallel%full_lat_start_idx, 45]
      full_lat_end_idx_a_grid(:)   = [44, parallel%full_lat_end_idx]
      full_lat_start_idx_c_grid = 1
      full_lat_end_idx_c_grid   = -1
      half_lat_start_idx_c_grid = 1
      half_lat_end_idx_c_grid   = -1

      full_lat_start_idx_no_pole_a_grid(:) = [parallel%full_lat_start_idx_no_pole, 45]
      full_lat_end_idx_no_pole_a_grid(:)   = [44, parallel%full_lat_end_idx_no_pole]
      full_lat_start_idx_no_pole_c_grid = 1
      full_lat_end_idx_no_pole_c_grid   = -1
      half_lat_start_idx_no_pole_c_grid = 1
      half_lat_end_idx_no_pole_c_grid   = -1
    else
      if (parallel%has_south_pole .and. use_pole_a_grid) then
        ! Assume pole_a_grid_tag(1) will be nonzero.
        do j = 2, size(pole_a_grid_tag)
          if (pole_a_grid_tag(j) == 0) then
            full_lat_start_idx_a_grid(1) = parallel%full_lat_start_idx
            full_lat_end_idx_a_grid(1) = parallel%full_lat_start_idx + j - 2
            full_lat_start_idx_c_grid = parallel%full_lat_start_idx + j - 1
            half_lat_start_idx_c_grid = parallel%half_lat_start_idx + j - 2

            full_lat_start_idx_no_pole_a_grid(1) = full_lat_start_idx_a_grid(1) + 1
            full_lat_end_idx_no_pole_a_grid(1) = full_lat_end_idx_a_grid(1)
            full_lat_start_idx_no_pole_c_grid = full_lat_start_idx_c_grid
            half_lat_start_idx_no_pole_c_grid = half_lat_start_idx_c_grid
            exit
          end if
        end do
      end if
      if (parallel%has_north_pole .and. use_pole_a_grid) then
        do j = 2, size(pole_a_grid_tag)
          if (pole_a_grid_tag(j) == 0) then
            full_lat_start_idx_a_grid(2) = parallel%full_lat_end_idx - j + 2
            full_lat_end_idx_a_grid(2) = parallel%full_lat_end_idx
            full_lat_end_idx_c_grid = parallel%full_lat_end_idx - j + 1
            half_lat_end_idx_c_grid = parallel%half_lat_end_idx - j + 2

            full_lat_start_idx_no_pole_a_grid(2) = full_lat_start_idx_a_grid(2)
            full_lat_end_idx_no_pole_a_grid(2) = full_lat_end_idx_a_grid(2) - 1
            full_lat_end_idx_no_pole_c_grid = full_lat_end_idx_c_grid
            half_lat_end_idx_no_pole_c_grid = half_lat_end_idx_c_grid
            exit
          end if
        end do
      end if
    end if

  end subroutine pole_a_grid_init

  subroutine pole_a_grid_init_state(state)

    type(state_type), intent(inout) :: state

    integer i, j

    ! NOTE: We assume state is set on A grid.
    do j = full_lat_start_idx_c_grid, full_lat_end_idx_c_grid
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state%u_c(i,j) = (state%u_a(i,j) + state%u_a(i+1,j)) * 0.5
      end do
    end do
    do j = half_lat_start_idx_c_grid, half_lat_end_idx_c_grid
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state%v_c(i,j) = (state%v_a(i,j) + state%v_a(i,j+1)) * 0.5
      end do
    end do

    call parallel_fill_halo(state%u_c, all_halo=.true.)
    call parallel_fill_halo(state%v_c, all_halo=.true.)

    ! TODO: Do we need to zero out the wind speeds in the void area for both of the grids?
    if (use_pole_a_grid) then
      do j = full_lat_start_idx_c_grid, full_lat_end_idx_c_grid
        state%u_a(:,j) = 0.0
        state%v_a(:,j) = 0.0
      end do
    else
      state%u_a(:,:) = 0.0
      state%v_a(:,:) = 0.0
    end if

    do j = full_lat_start_idx_a_grid(south_pole), full_lat_end_idx_a_grid(south_pole)
      state%u_c(:,j) = 0.0
      state%v_c(:,j) = 0.0
    end do

    do j = full_lat_start_idx_a_grid(north_pole), full_lat_end_idx_a_grid(north_pole)
      state%u_c(:,j) = 0.0
      state%v_c(:,j) = 0.0
    end do

  end subroutine pole_a_grid_init_state

end module pole_a_grid_mod