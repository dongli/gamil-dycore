module log_mod

  use hash_table_mod
  use time_mod
  use string_mod

  implicit none

  type(hash_table_type) diags

  interface log_add_diag
    module procedure log_add_diag_1
    module procedure log_add_diag_2
  end interface log_add_diag

contains

  subroutine log_add_diag_1(name, value)

    character(*), intent(in) :: name
    integer, intent(in) :: value

    call diags%insert(name, value)

  end subroutine log_add_diag_1

  subroutine log_add_diag_2(name, value)

    character(*), intent(in) :: name
    real, intent(in) :: value

    call diags%insert(name, value)

  end subroutine log_add_diag_2

  subroutine log_notice(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    if (present(file) .and. present(line)) then
      write(6, *) '[Notice]: ' // trim(file) // ': ' // to_string(line) // ': ' // trim(message)
    else
      write(6, *) '[Notice]: ' // trim(message)
    end if

  end subroutine log_notice

  subroutine log_warning(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    if (present(file) .and. present(line)) then
      write(6, *) '[Warning]: ' // trim(file) // ': ' // to_string(line) // ': ' // trim(message)
    else
      write(6, *) '[Warning]: ' // trim(message)
    end if

  end subroutine log_warning

  subroutine log_error(message, file, line)

    character(*), intent(in) :: message
    character(*), intent(in), optional :: file
    integer, intent(in), optional :: line

    if (present(file) .and. present(line)) then
      write(6, *) '[Error]: ' // trim(file) // ': ' // to_string(line) // ': ' // trim(message)
    else
      write(6, *) '[Error]: ' // trim(message)
    end if
    stop 1

  end subroutine log_error

  subroutine log_step()

    type(hash_table_iterator_type) iter

    write(6, '(" => ", A)', advance='no') trim(curr_time_format)

    iter = hash_table_iterator(diags)
    do while (.not. iter%ended())
      select type(value => iter%value)
      type is (integer)
        write(6, '(X, A)', advance='no') trim(to_string(value))
      type is (real)
        write(6, '(X, A)', advance='no') trim(to_string(value, 20))
      end select
      call iter%next()
    end do
    write(6, *)

  end subroutine log_step

end module log_mod