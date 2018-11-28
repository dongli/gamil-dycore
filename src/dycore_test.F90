program dycore_test

  use params_mod
  use io_mod
  use time_mod
  use dycore_mod
  use rossby_haurwitz_wave_test_mod
  use steady_geostrophic_flow_test_mod
  use mountain_zonal_flow_test_mod
  use jet_zonal_flow_test_mod
  use shallow_water_waves_test_mod

  character(256) namelist_file_path

  if (command_argument_count() /= 1) then
    write(6, *) 'Usage: ./qconswm <namelist_file_path>'
    stop 1
  end if

  call get_command_argument(1, namelist_file_path)

  call params_read(namelist_file_path)

  call dycore_init()

  if (is_restart_run) then
    call dycore_restart()
  else
    select case (test_case)
    case ('steady_geostrophic_flow')
      call steady_geostrophic_flow_test_set_initial_condition()
    case ('rossby_haurwitz_wave')
      call rossby_haurwitz_wave_test_set_initial_condition()
    case ('mountain_zonal_flow')
      call mountain_zonal_flow_test_set_initial_condition()
    case ('jet_zonal_flow')
      call jet_zonal_flow_test_set_initial_condition()
    case ('shallow_water_waves')
      call shallow_water_waves_test_set_initial_condition()
    case default
      write(6, *) '[Error]: Unknown test case ' // trim(test_case) // '!'
    end select
  end if

  call dycore_run()

  call dycore_final()

end program dycore_test
