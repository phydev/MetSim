module sim_init_m

  use global_m
  use misc_m

  implicit none

  private

  public :: single_cell_init, cell_pos_init, space_init, parameters_init, print_header

  contains
   subroutine  single_cell_init(cell, tcell, ncell, lxyz_part, lxyz_inv_part, sphere, np_sphere, np_part, first)

      implicit none

      ! input variables
      type(mesh_t), allocatable, intent(inout) :: cell(:,:)
      integer, allocatable, intent(in) :: lxyz_part(:,:), lxyz_inv_part(:,:)
      integer, intent(in) ::  tcell, sphere(:,:), np_sphere, np_part, ncell(:)
      logical, intent(in) :: first

      ! internal only
      integer :: ip, i, j, l, icell, itype

      itype = 1
      i = 1
      if(first) then
         cell(0:np_part,1)%phi = -1.d0
         cell(0:np_part,1)%itype = 1.d0
         do l = 1, np_sphere
            cell(lxyz_inv_part(sphere(l,1),sphere(l,2)),1)%phi = 1.d0
            cell(lxyz_inv_part(sphere(l,1),sphere(l,2)),1)%itype = 1.d0
         end do
      else
         do icell=2,tcell
            i = i +1
            cell(0:np_part,icell)%phi = cell(0:np_part,1)%phi
            if(i.gt.ncell(itype)) then
               itype = itype + 1
               i = 1
            end if
            cell(0:np_part,icell)%itype = itype

         end do

      end if

    end subroutine single_cell_init



    subroutine cell_pos_init(r, icell, nleap,  dr, dri, drf, ncell, cell_radius, lxyz, lxyz_inv,&
         Lsize, np, sphere, np_sphere, iseed)

      implicit none

      ! input variables
      integer, allocatable, intent(inout) :: r(:)
      integer, intent(inout) ::  dr(2), icell
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:)
      integer, intent(in) ::  np, iseed, Lsize(2),  ncell, sphere(:,:),&
           np_sphere, drf(2), dri(2), nleap
      real, intent(in) :: cell_radius

      ! internal only
      integer :: ip, i, j, l


      i=-Lsize(1)+dri(1)
      do while  (i<= Lsize(1)-drf(1))
          j= -Lsize(2)+dri(2)
         do while (j<= Lsize(2)-drf(2))

               ip = lxyz_inv(i,j)

               r(icell) = ip

               icell = icell + 1
               j = j + dr(2)

               !if(icell.ge.nleap) dr(2) = 2*cell_radius
            if (icell.gt.ncell) EXIT

         end do
         i = i + dr(1)
         if (icell.gt.ncell) EXIT
      end do



    end subroutine cell_pos_init

    subroutine space_init(Lsize, lxyz, lxyz_inv, np_bndry, np, periodic)

      implicit none

      ! input/output variables
      integer, intent(in) ::  Lsize(1:2), np_bndry
      integer, allocatable, intent(inout) :: lxyz(:,:), lxyz_inv(:,:)
      ! internal variables
      integer :: iint, jint, l, m, n, ip, ip_part, np, ki, kf
      real :: hs(1:2), i, j
      logical :: boundary, periodic

      ip = 0
      ip_part = np

      ! bulk points
      ! allocated from 1 to np

      i=-Lsize(1)
      do while (i<=Lsize(1)-1)
         j=-Lsize(2)
         do while(j<= Lsize(2)-1)


               ip = ip + 1

               lxyz(ip,1) = iint
               lxyz(ip,2) = jint
               lxyzf(ip,1) = i
               lxyzf(ip,2) = j

               lxyz_inv(iint,jint) = ip
               jint = jint + 1
               j = j + dr(2)
         end do
         iint = iint + 1
         i = i + dr(1)
      end do

      ! boundary points
      ! allocated from np to np_part

      do i=-Lsize(1)-np_bndry, Lsize(1)-1+np_bndry
         do j=-Lsize(2)-np_bndry, Lsize(2)-1+np_bndry

               l = i
               m = j

               boundary = .false.

               hs(1) = heaviside(real(i))
               hs(2) = heaviside(real(j))

               if( abs(i)>Lsize(1)-hs(1)) then
                  boundary = .true.
                  if(periodic) then
                     l = i - SIGN(1,i)*(2*Lsize(1))
                  else
                     l = i - SIGN(1,i)
                  end if
               end if

               if( abs(j)>Lsize(2)-hs(2)) then
                  boundary = .true.

                  if(periodic) then
                     m = j  - SIGN(1,j)*(2*Lsize(2))
                  else
                     m = j - SIGN(1,j)
                  end if
               end if



               if(boundary) then
                  if(periodic) then
                     ip_part = ip_part + 1

                     lxyz(ip_part,1) = l
                     lxyz(ip_part,2) = m

                     lxyz_inv(i,j) = lxyz_inv(l, m)
                  else
                     !ip_part = ip_part + 1
                     !lxyz(ip_part,1) = l
                     !lxyz(ip_part,2) = m

                     lxyz_inv(i,j) =  0 !lxyz_inv(l, m)

                  end if
               end if



               ! the updating of boundaries should be
               ! cell(ip_part) = cell( lxyz_inv( lxyz(np_part,1), lxyz(np_part,2), lxyz(np_part,3) ) )
               ! because:
               ! lxyz(np_part,1), lxyz(np_part,2), lxyz(np_part,3) ->(l, m, n) in the bulk
               ! then lxyz_inv(m,n,l) will give the updated value in the bulk
               ! for the respective boundary position cell(ip_part)

         end do
      end do

    end subroutine space_init


    subroutine parameters_init(cell_radius, ntypes, density,  interface_width, tstep, dt, Lsize, dr, dir_name, iseed,&
        np_bndry, depletion_weight,adh1, adh2, output_period, periodic)

      implicit none

      real, allocatable, intent(inout) :: density(:)
      real, intent(inout) :: cell_radius,  interface_width, dt, depletion_weight, adh1, adh2
      integer, intent(inout) :: tstep, Lsize(2), iseed, np_bndry, dr(2), output_period, ntypes
      character(len=3), intent(inout) :: dir_name
      character(len=255) :: temp
      logical :: periodic

      OPEN (UNIT=1,FILE='input_file')
      read(1,*) cell_radius, temp ! R_c - Cell Radius
      read(1,*) ntypes, temp ! cell types

      ALLOCATE(density(ntypes)) ! allocating the

      read(1,*) density(1:ntypes), temp ! density
      if(sum(density(1:ntypes)).gt.1.0) then
         write(*,'(A)') "Warnning:"
         write(*,'(A)') "The sum of each cell type density must be less or equal 1.0!"
         write(*,'(A)') "The program will be initiate with the default density value."
         density(1:ntypes) = 1.0/real(ntypes)
      end if
      read(1,*) interface_width, temp ! Eps - Interface Witdh
      read(1,*) tstep, temp ! tstep - Total time step
      read(1,*) dt, temp ! dt - Time increment
      read(1,*) Lsize(1:2), dr(1:2), temp !Box Length - x,y,z, dr
      read(1,*) dir_name, temp ! Simulation name
      read(1,*) iseed, temp ! Initial Seed for RAN2
      read(1,*) np_bndry, temp ! boundary points
      read(1,*) depletion_weight, temp ! repulsive force coefficient
      read(1,*) adh1, temp ! adhesion const 1
      read(1,*) adh2, temp ! adhesion const 2
      read(1,*) output_period, temp ! the period which the data will be written to the output
      read(1,*) periodic, temp ! boundary conditions
      CLOSE(1)
      call system('mkdir '//dir_name)

      OPEN (UNIT=2,FILE=dir_name//'/parameters'//dir_name//'.init')
      write(2,'(F10.2,A)') cell_radius, " cell_radius" ! R_c - Cell Radius
      write(2,'(I10,A)') ntypes, "ntypes" ! ntypes - number of cell types
      write(2,*) density(1:ntypes), " density" ! density
      write(2,'(F10.2,A)') interface_width, " interface_width" ! Eps - Interface Witdh
      write(2,'(I10,A)')   tstep, " tstep" ! tstep - Total time step
      write(2,'(F10.5,A)') dt, " dt" ! dt - Time increment
      write(2,'(I3, I3,  I3, I3, A)') Lsize(1:2), dr(1:2), " box_length_xyz_dr" !Box Length - x,y,z, dr
      write(2,'(A,A)') dir_name, " dir_name" ! Simulation name
      write(2,'(I10,A)') iseed, " iseed" ! Initial Seed for RAN2
      write(2,'(I10,A)') np_bndry, " bounary_points"! boundary points
      write(2,'(F10.2,A)') depletion_weight, " depletion_weight" ! repulsive force coefficient
      write(2,'(F10.2,A)') adh1, " adh1" ! constant adhesion 1
      write(2,'(F10.2,A)') adh2, " adh2" ! constant adhesion 2
      write(2,'(I10,A)') output_period, " output_period" ! the period which the data will be written to the output
      write(2,'(L1,A) ') periodic, " periodic" ! boundary conditions
      CLOSE(2)
    end subroutine parameters_init



    subroutine print_header(Lsize, tcell, ntypes, ncell, dir_name, periodic)

      implicit none
      integer, intent(in) :: Lsize(1:2), ntypes, ncell(:), tcell
      character(len=255) :: cwd, hostname
      character(len=32) :: username
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      character(3)  :: dir_name
      integer, dimension(8) :: values
      logical :: periodic

      call date_and_time(date,time,zone,values)
      call date_and_time(DATE=date,ZONE=zone)
      call date_and_time(TIME=time)

      call hostnm(hostname)
      call getcwd(cwd)
      call getlog(username)
      write(*,'(A)') "                                Running Cell"
      write(*,'(A)') "       "
      write(*,'(A)') "Version        :       0.0.1"
      write(*,'(A,A)') "Locate         :       ", trim(cwd)
      write(*,'(A,A)') "User           :       ", trim(username)
      write(*,'(A)') "Developer      :       Moreira, M."
      write(*,'(A)') "       "
      write(*,'(A,A)') "                      The code is running in ", trim(hostname)
      write(*,'(A)') "       "
      write(*,'(A)') "       "
      write(*,'(A,2X,A,2X,A)') "             Calculation started on", date(7:8)//"/"//date(5:6)//"/"//date(1:4),&
           time(1:2)//":"//time(3:4)//":"//time(5:6)
      write(*,'(A)') "       "
      write(*,'(A)') "************************************ Grid *************************************"
      write(*,'(A)') "Simulation Box:"
      write(*,'(A,I3,A,I3,A,I3,A)') "  Lengths = (",Lsize(1),",", Lsize(2), ")"
      write(*,'(A)') "  the code will run in 3 dimension(s)."
      if(periodic) then
         write(*,'(A)') "  the code will treat the system as periodic in 3 dimension(s)."
      else
         write(*,'(A)') "  the code will treat the system as periodic in 2 dimension(s)."
      end if
      write(*,'(A)') "*******************************************************************************"
      write(*,'(A,I3)') "  Number of cells:", tcell
      write(*,'(A,I3)') "  Cell types:", ntypes
      write(*,*) " Cells per type:", ncell(1:ntypes)
      write(*,'(A,A)') "  Simulation ID: ", dir_name
    end subroutine print_header

end module sim_init_m
