module pole_a_grid_mod

  use mesh_mod
  use parallel_mod
  use log_mod
  use params_mod

  implicit none

  private

  public pole_a_grid_init
  public south_ac_lat_idx
  public north_ac_lat_idx

  integer south_ac_lat_idx
  integer north_ac_lat_idx

contains

  subroutine pole_a_grid_init()

    integer j

    if (parallel%has_south_pole) then
      do j = 1, size(pole_a_grid_tag)
        if (pole_a_grid_tag(j) == 0) then
          south_ac_lat_idx = parallel%full_lat_start_idx + j - 1
          exit
        end if
      end do
    end if
    if (parallel%has_north_pole) then
      do j = 1, size(pole_a_grid_tag)
        if (pole_a_grid_tag(j) == 0) then
          north_ac_lat_idx = parallel%full_lat_end_idx - j + 1
          exit
        end if
      end do
    end if

  end subroutine pole_a_grid_init

end module pole_a_grid_mod