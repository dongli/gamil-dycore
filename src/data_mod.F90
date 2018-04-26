module data_mod

  use log_mod
  use mesh_mod
  use params_mod
  use types_mod

  implicit none

  private

  public coef
  public state
  public static
  public tend
  public data_init
  public data_final

  type(coef_type) coef
  type(state_type) state(-1:2)
  type(static_type) static
  type(tend_type) tend(0:2)

contains

  subroutine data_init()

    integer i, j, time_idx

    call allocate_data(coef)
    do j = 1, mesh%num_full_lat
      coef%cori(j) = 2.0 * omega * mesh%full_sin_lat(j)
      if (j == 1 .or. j == mesh%num_full_lat) then
        coef%curv(j) = 0.0
      else
        coef%curv(j) = mesh%full_sin_lat(j) / mesh%full_cos_lat(j) / radius
      end if
      coef%full_dlon(j) = 2.0 * radius * mesh%dlon * mesh%full_cos_lat(j)
      coef%full_dlat(j) = 2.0 * radius * mesh%dlat * mesh%full_cos_lat(j)
    end do

    do j = 1, mesh%num_half_lat
      coef%half_dlon(j) = 2.0 * radius * mesh%dlon * mesh%half_cos_lat(j)
      coef%half_dlat(j) = 2.0 * radius * mesh%dlat * mesh%half_cos_lat(j)
    end do

    do time_idx = 0, 2
      call allocate_data(state(time_idx))
    end do
    do time_idx = 0, 2
      call allocate_data(tend(time_idx))
    end do
    call allocate_data(static)

    call log_notice('Data module is initialized.')

  end subroutine data_init

  subroutine data_final()

    integer time_idx

    call deallocate_data(coef)
    do time_idx = lbound(state, 1), ubound(state, 1)
      call deallocate_data(state(time_idx))
    end do
    do time_idx = lbound(tend, 1), ubound(tend, 1)
      call deallocate_data(tend(time_idx))
    end do
    call deallocate_data(static)

    call log_notice('Data module is finalized.')

  end subroutine data_final

end module data_mod