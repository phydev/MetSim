program main

  use run_cells_m
  use global_m

  implicit none

  call get_command_argument(1,sim_id)
  call get_command_argument(2,arg_iseed)
  call get_command_argument(3,arg_porosity)

  read (arg_porosity,'(F10.4)') porosity
  read (arg_iseed,'(I10)') iseed

  call run_cells(sim_id, porosity, iseed)



end program main
