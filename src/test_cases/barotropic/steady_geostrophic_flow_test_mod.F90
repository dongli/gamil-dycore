module steady_geostrophic_flow_test_mod

  use mesh_mod
  use parallel_mod
  use params_mod
  use data_mod

  implicit none

  private

  public steady_geostrophic_flow_test_set_initial_condition

  real, parameter :: alpha = 0.0
  real, parameter :: u0 = 2 * pi * radius / (12 * 86400)
  real, parameter :: gd0 = 2.94e4 ! m2 s-2

contains

  subroutine steady_geostrophic_flow_test_set_initial_condition()

    real cos_lat, sin_lat, cos_lon, sin_lon, cos_alpha, sin_alpha
    integer i, j

    write(6, *) '[Notice]: Use steady geostrophic flow initial condition.'

    cos_alpha = cos(alpha)
    sin_alpha = sin(alpha)
    static%ghs(:,:) = 0.0

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        cos_lon = mesh%full_cos_lon(i)
        state(1)%u_a(i,j) = u0 * (cos_lat * cos_alpha + cos_lon * sin_lat * sin_alpha)
      end do
    end do

    call parallel_fill_halo(state(1)%u_a, all_halo=.true.)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sin_lon = mesh%full_cos_lon(i)
        state(1)%v_a(i,j) = - u0 * sin_lon * sin_alpha
      end do
    end do

    call parallel_fill_halo(state(1)%v_a, all_halo=.true.)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        cos_lon = mesh%full_cos_lon(i)
        state(1)%gd(i,j) = gd0 - (radius * omega * u0 + u0**2 * 0.5) * (sin_lat * cos_alpha - cos_lon * cos_lat * sin_alpha)**2
      end do
    end do

    call parallel_fill_halo(state(1)%gd, all_halo=.true.)

  end subroutine steady_geostrophic_flow_test_set_initial_condition

end module steady_geostrophic_flow_test_mod
