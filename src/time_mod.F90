module time_mod

  use datetime_mod
  use timedelta_mod
  use map_mod
  use params_mod, &
    start_time_in => start_time, &
    end_time_in => end_time, &
    time_step_size_in => time_step_size

  implicit none

  private

  public time_init
  public time_reset_start_time
  public time_swap_indices
  public time_advance
  public time_elapsed_seconds
  public time_is_finished
  public time_add_alert
  public time_is_alerted

  public curr_time
  public start_time_format
  public curr_time_format
  public time_step
  public old_time_idx
  public new_time_idx

  type alert_type
    type(timedelta_type) period
    type(datetime_type) last_time
  end type alert_type

  type(datetime_type) start_time
  type(datetime_type) end_time
  type(datetime_type) curr_time
  type(timedelta_type) time_step_size
  real(8) elapsed_seconds
  type(map_type) alerts
  integer time_step
  integer old_time_idx
  integer new_time_idx
  character(30) start_time_format
  character(30) curr_time_format

contains

  subroutine time_init()

    if (days > 0 .or. hours > 0 .or. minutes > 0) then
      start_time = datetime(days=0, hours=0, minutes=0)
    else if (sum(start_time_in) > 0) then
      start_time = datetime(year=start_time_in(1), month=start_time_in(2), day=start_time_in(3), &
        hour=start_time_in(4), minute=start_time_in(5))
    end if
    if (days > 0 .or. hours > 0 .or. minutes > 0) then
      end_time = datetime(days=days, hours=hours, minutes=minutes)
    else if (sum(end_time_in) > 0) then
      end_time = datetime(year=end_time_in(1), month=end_time_in(2), day=end_time_in(3), &
        hour=end_time_in(4), minute=end_time_in(5))
    end if

    time_step = 0
    elapsed_seconds = 0
    old_time_idx = 1
    new_time_idx = 2
    time_step_size = timedelta(seconds=time_step_size_in)

    curr_time = start_time

    start_time_format = start_time%isoformat()
    curr_time_format = curr_time%isoformat()

  end subroutine time_init

  subroutine time_reset_start_time(time)

    type(datetime_type), intent(in) :: time

    start_time = time
    curr_time = start_time

    start_time_format = start_time%isoformat()
    curr_time_format = curr_time%isoformat()

  end subroutine time_reset_start_time

  subroutine time_swap_indices(i, j)

    integer, intent(inout) :: i
    integer, intent(inout) :: j

    integer tmp

    tmp = i
    i = j
    j = tmp

  end subroutine time_swap_indices

  subroutine time_advance()

    call time_swap_indices(old_time_idx, new_time_idx)

    time_step = time_step + 1
    elapsed_seconds = elapsed_seconds + time_step_size%total_seconds()
    curr_time = curr_time + time_step_size
    curr_time_format = curr_time%isoformat()

  end subroutine time_advance

  real function time_elapsed_seconds() result(res)

    res = elapsed_seconds

  end function time_elapsed_seconds

  logical function time_is_finished() result(res)

    res = curr_time >= end_time

  end function time_is_finished

  subroutine time_add_alert(name, days, hours, minutes, seconds)

    character(*), intent(in) :: name
    class(*), intent(in), optional :: days
    class(*), intent(in), optional :: hours
    class(*), intent(in), optional :: minutes
    class(*), intent(in), optional :: seconds

    type(alert_type) alert

    alert%period = timedelta(days, hours, minutes, seconds)
    alert%last_time = start_time - alert%period
    call alerts%insert(trim(name), alert)

  end subroutine time_add_alert

  function time_is_alerted(name) result(res)

    character(*), intent(in) :: name
    logical res

    type(alert_type), pointer :: alert => null()
    type(datetime_type) time

    alert => get_alert(name)
    if (associated(alert)) then
      time = alert%last_time + alert%period
      if (time <= curr_time) then
        alert%last_time = curr_time
        res = .true.
      else
        res = .false.
      end if
    else
      res = .false.
    end if

  end function time_is_alerted

  function get_alert(name) result(res)

    character(*), intent(in) :: name
    type(alert_type), pointer :: res

    class(*), pointer :: value

    if (alerts%mapped(name)) then
      value => alerts%value(name)
      select type (value)
      type is (alert_type)
        res => value
      end select
    else
      nullify(res)
    end if

  end function get_alert

end module time_mod
