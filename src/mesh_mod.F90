module mesh_mod

  use log_mod
  use params_mod

  implicit none

  private

  public mesh
  public mesh_init
  public mesh_final

  type mesh_type
    integer num_full_lon
    integer num_half_lon
    integer num_full_lat
    integer num_half_lat
    real dlon
    real dlat
    real, allocatable :: full_lon(:)
    real, allocatable :: half_lon(:)
    real, allocatable :: full_lat(:)
    real, allocatable :: half_lat(:)
    real, allocatable :: full_cos_lat(:)
    real, allocatable :: half_cos_lat(:)
    real, allocatable :: full_sin_lat(:)
    real, allocatable :: half_sin_lat(:)
    ! For output
    real, allocatable :: full_lon_deg(:)
    real, allocatable :: half_lon_deg(:)
    real, allocatable :: full_lat_deg(:)
    real, allocatable :: half_lat_deg(:)
  end type mesh_type

  type(mesh_type) mesh

contains

  subroutine mesh_init()

    integer i, j

    mesh%num_full_lon = num_lon
    mesh%num_half_lon = num_lon
    mesh%num_full_lat = num_lat
    mesh%num_half_lat = num_lat - 1

    allocate(mesh%full_lon(mesh%num_full_lon))
    allocate(mesh%half_lon(mesh%num_half_lon))
    allocate(mesh%full_lat(mesh%num_full_lat))
    allocate(mesh%half_lat(mesh%num_half_lat))
    allocate(mesh%full_cos_lat(mesh%num_full_lat))
    allocate(mesh%half_cos_lat(mesh%num_half_lat))
    allocate(mesh%full_sin_lat(mesh%num_full_lat))
    allocate(mesh%half_sin_lat(mesh%num_half_lat))
    allocate(mesh%full_lon_deg(mesh%num_full_lon))
    allocate(mesh%half_lon_deg(mesh%num_half_lon))
    allocate(mesh%full_lat_deg(mesh%num_full_lat))
    allocate(mesh%half_lat_deg(mesh%num_half_lat))

    mesh%dlon = 2 * pi / mesh%num_full_lon
    do i = 1, mesh%num_full_lon
      mesh%full_lon(i) = (i - 1) * mesh%dlon
      mesh%half_lon(i) = mesh%full_lon(i) + 0.5 * mesh%dlon
      mesh%full_lon_deg(i) = mesh%full_lon(i) * rad_to_deg
      mesh%half_lon_deg(i) = mesh%half_lon(i) * rad_to_deg
    end do

    mesh%dlat = pi / mesh%num_half_lat
    do j = 1, mesh%num_half_lat
      mesh%full_lat(j) = - 0.5 * pi + (j - 1) * mesh%dlat
      mesh%half_lat(j) = mesh%full_lat(j) + 0.5 * mesh%dlat
      mesh%full_lat_deg(j) = mesh%full_lat(j) * rad_to_deg
      mesh%half_lat_deg(j) = mesh%half_lat(j) * rad_to_deg
    end do
    mesh%full_lat(num_lat) = 0.5 * pi
    mesh%full_lat_deg(num_lat) = 90.0

    do j = 1, mesh%num_half_lat
      mesh%half_cos_lat(j) = cos(mesh%half_lat(j))
      mesh%half_sin_lat(j) = sin(mesh%half_lat(j))
    end do

    do j = 1, mesh%num_full_lat
      mesh%full_cos_lat(j) = cos(mesh%full_lat(j))
      mesh%full_sin_lat(j) = sin(mesh%full_lat(j))
    end do
    mesh%full_cos_lat(1) = 0.0
    mesh%full_cos_lat(mesh%num_full_lat) = 0.0
    mesh%full_sin_lat(1) = -1.0
    mesh%full_sin_lat(mesh%num_full_lat) = 1.0

    call log_notice('Mesh module is initialized.')

  end subroutine mesh_init

  subroutine mesh_final()

    if (allocated(mesh%full_lon)) deallocate(mesh%full_lon)
    if (allocated(mesh%full_lat)) deallocate(mesh%full_lat)
    if (allocated(mesh%half_lon)) deallocate(mesh%half_lon)
    if (allocated(mesh%half_lat)) deallocate(mesh%half_lat)
    if (allocated(mesh%full_cos_lat)) deallocate(mesh%full_cos_lat)
    if (allocated(mesh%half_cos_lat)) deallocate(mesh%half_cos_lat)
    if (allocated(mesh%full_sin_lat)) deallocate(mesh%full_sin_lat)
    if (allocated(mesh%half_sin_lat)) deallocate(mesh%half_sin_lat)
    if (allocated(mesh%full_lon_deg)) deallocate(mesh%full_lon_deg)
    if (allocated(mesh%half_lon_deg)) deallocate(mesh%half_lon_deg)
    if (allocated(mesh%full_lat_deg)) deallocate(mesh%full_lat_deg)
    if (allocated(mesh%half_lat_deg)) deallocate(mesh%half_lat_deg)

    call log_notice('Mesh module is finalized.')

  end subroutine mesh_final

end module mesh_mod
