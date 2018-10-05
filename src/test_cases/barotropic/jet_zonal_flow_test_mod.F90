module jet_zonal_flow_test_mod

  use mesh_mod
  use parallel_mod
  use params_mod
  use data_mod
  use log_mod
  use string_mod

  implicit none

  private

  public jet_zonal_flow_test_set_initial_condition

  real, parameter :: u_max = 80.0
  real, parameter :: lat0 = pi / 7.0
  real, parameter :: lat1 = pi / 2.0 - lat0
  real, parameter :: en = exp(-4.0 / (lat1 - lat0)**2)
  real, parameter :: gh0 = g * 1.0e4
  real, parameter :: ghd = g * 120
  real, parameter :: lat2 = pi / 4.0
  real, parameter :: alpha = 1.0 / 3.0
  real, parameter :: beta = 1.0 / 15.0

contains

  subroutine jet_zonal_flow_test_set_initial_condition()

    integer i, j, neval, ierr
    real abserr

    write(6, *) '[Notice]: Use jet zonal flow initial condition.'

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        static%ghs(i,j) = 0.0
      end do
    end do

    call parallel_fill_halo(static%ghs, all_halo=.true.)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state(1)%u_c(i,j) = u_function(mesh%full_lat(j))
      end do
    end do

    call parallel_fill_halo(state(1)%u_c, all_halo=.true.)

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        state(1)%u_d(i,j) = u_function(mesh%half_lat(j))
      end do
    end do

    call parallel_fill_halo(state(1)%u_d, all_halo=.true.)

    state(1)%v_c(:,:) = 0.0
    state(1)%v_d(:,:) = 0.0

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      if (j == parallel%full_lat_start_idx) then
        state(1)%gd(0,j) = gh0
      else
        call qags(gh_integrand, -0.5*pi, mesh%full_lat(j), 1.0e-10, 1.0e-3, state(1)%gd(0,j), abserr, neval, ierr)
        if (ierr /= 0) then
          call log_error('Failed to calculate integration at (' // to_string(i) // ',' // to_string(j) // ')!')
        end if
        state(1)%gd(0,j) = gh0 - state(1)%gd(0,j)
      end if
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        state(1)%gd(i,j) = state(1)%gd(0,j)
        ! Add perturbation.
        state(1)%gd(i,j) = state(1)%gd(i,j) + ghd * &
          cos(mesh%full_lat(j)) * &
          exp(-((mesh%full_lon(i) - pi)  / alpha)**2) * &
          exp(-((lat2 - mesh%full_lat(j)) / beta)**2)
      end do
    end do

    call parallel_fill_halo(state(1)%gd, all_halo=.true.)

  end subroutine jet_zonal_flow_test_set_initial_condition

  real function gh_integrand(lat) result(res)

    real, intent(in) :: lat

    real u_c, f

    u_c = u_function(lat)
    f = 2 * omega * sin(lat)
    res = radius * u_c * (f + tan(lat) / radius * u_c)

  end function gh_integrand

  real function u_function(lat) result(res)

    real, intent(in) :: lat

    if (lat <= lat0 .or. lat >= lat1) then
      res = 0.0
    else
      res = u_max / en * exp(1 / (lat - lat0) / (lat - lat1))
    end if

  end function u_function

end module jet_zonal_flow_test_mod
