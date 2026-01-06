module barnes_hut
    use geometry
    use particle
    !$ use omp_lib
    implicit none

    ! Barnes–Hut opening angle
    real(dp), parameter :: theta = 1.0_dp

    ! Number of particles
    integer :: n

    ! Particle data (positions, velocities, masses)
    type(particle3d), allocatable :: particles(:)

    ! Accelerations
    type(vector3d), allocatable :: acc(:)

    !=========================================================
    ! Definition of spatial range (min/max corners of a cell)
    !=========================================================
    type :: range
        type(point3d) :: min, max
    end type range

    type :: cptr
        class(cell), pointer :: ptr => null()
    end type cptr

    type :: cell
        type(range) :: range
        type(point3d) :: part
        integer :: pos
        integer :: type      ! 0 = empty, 1 = particle, 2 = conglomerate
        real(dp) :: mass
        type(point3d) :: c_o_m
        type(cptr), dimension(2,2,2) :: subcell
    end type cell

contains

    !=========================================================
    ! calculate_ranges
    !
    ! Computes the bounding box containing all particles and
    ! assigns it to the root cell.
    !=========================================================
    subroutine calculate_ranges(root)
        type(cell), pointer :: root
        real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax
        real(dp) :: span
        integer :: i

        xmin = particles(1)%p%x
        xmax = particles(1)%p%x
        ymin = particles(1)%p%y
        ymax = particles(1)%p%y
        zmin = particles(1)%p%z
        zmax = particles(1)%p%z

        do i = 2, n
            xmin = min(xmin, particles(i)%p%x)
            xmax = max(xmax, particles(i)%p%x)
            ymin = min(ymin, particles(i)%p%y)
            ymax = max(ymax, particles(i)%p%y)
            zmin = min(zmin, particles(i)%p%z)
            zmax = max(zmax, particles(i)%p%z)
        end do

        span = max(xmax-xmin, max(ymax-ymin, zmax-zmin)) * 1.1_dp

        root%range%min = point3d( (xmin+xmax)/2 - span/2, &
                                  (ymin+ymax)/2 - span/2, &
                                  (zmin+zmax)/2 - span/2 )

        root%range%max = point3d( (xmin+xmax)/2 + span/2, &
                                  (ymin+ymax)/2 + span/2, &
                                  (zmin+zmax)/2 + span/2 )
    end subroutine calculate_ranges

    !=========================================================
    ! nullify_pointers
    !
    ! Sets all subcell pointers of a cell to null().
    !=========================================================
    subroutine nullify_pointers(c)
        type(cell), pointer :: c
        integer :: i, j, k

        do i = 1,2
            do j = 1,2
                do k = 1,2
                    nullify(c%subcell(i,j,k)%ptr)
                end do
            end do
        end do
    end subroutine nullify_pointers

    !=========================================================
    ! belongs
    !
    ! Returns .true. if point p lies inside cell c.
    !=========================================================
    logical function belongs(p, c)
        type(point3d), intent(in) :: p
        type(cell), pointer, intent(in) :: c

        belongs = (p%x >= c%range%min%x .and. p%x <= c%range%max%x .and. &
                   p%y >= c%range%min%y .and. p%y <= c%range%max%y .and. &
                   p%z >= c%range%min%z .and. p%z <= c%range%max%z)
    end function belongs

    !=========================================================
    ! find_cell
    !
    ! Finds the leaf cell where a particle should be placed.
    !=========================================================
    recursive subroutine find_cell(root, goal, p)
        type(cell), pointer :: root, goal, temp
        type(point3d), intent(in) :: p
        integer :: i, j, k

        select case (root%type)
        case (2)
            do i = 1,2
                do j = 1,2
                    do k = 1,2
                        if (associated(root%subcell(i,j,k)%ptr)) then
                            if (belongs(p, root%subcell(i,j,k)%ptr)) then
                                call find_cell(root%subcell(i,j,k)%ptr, temp, p)
                                goal => temp
                                return
                            end if
                        end if
                    end do
                end do
            end do
        case default
            goal => root
        end select
    end subroutine find_cell

    !=========================================================
    ! place_cell
    !
    ! Places a particle into a leaf cell. If the cell already
    ! contains a particle, it is subdivided and both particles
    ! are placed into the appropriate subcells.
    !=========================================================
    recursive subroutine place_cell(c, p, idx)
        type(cell), pointer :: c, temp
        type(point3d), intent(in) :: p
        integer, intent(in) :: idx

        select case (c%type)
        case (0)
            c%type = 1
            c%part = p
            c%pos  = idx

        case (1)
            call create_subcells(c)
            call find_cell(c, temp, p)
            call place_cell(temp, p, idx)
        end select
    end subroutine place_cell

    !=========================================================
    ! create_subcells
    !
    ! Subdivides a leaf cell into 8 octants.
    !=========================================================
    subroutine create_subcells(c)
        type(cell), pointer :: c
        type(point3d) :: p_old, mid
        integer :: i, j, k

        p_old = c%part
        mid = point3d( (c%range%min%x + c%range%max%x)/2, &
                    (c%range%min%y + c%range%max%y)/2, &
                    (c%range%min%z + c%range%max%z)/2 )

        c%type = 2

        do i = 1,2
            do j = 1,2
                do k = 1,2
                    allocate(c%subcell(i,j,k)%ptr)
                    call nullify_pointers(c%subcell(i,j,k)%ptr)

                    ! X range
                    if (i == 1) then
                        c%subcell(i,j,k)%ptr%range%min%x = c%range%min%x
                        c%subcell(i,j,k)%ptr%range%max%x = mid%x
                    else
                        c%subcell(i,j,k)%ptr%range%min%x = mid%x
                        c%subcell(i,j,k)%ptr%range%max%x = c%range%max%x
                    end if

                    ! Y range
                    if (j == 1) then
                        c%subcell(i,j,k)%ptr%range%min%y = c%range%min%y
                        c%subcell(i,j,k)%ptr%range%max%y = mid%y
                    else
                        c%subcell(i,j,k)%ptr%range%min%y = mid%y
                        c%subcell(i,j,k)%ptr%range%max%y = c%range%max%y
                    end if

                    ! Z range
                    if (k == 1) then
                        c%subcell(i,j,k)%ptr%range%min%z = c%range%min%z
                        c%subcell(i,j,k)%ptr%range%max%z = mid%z
                    else
                        c%subcell(i,j,k)%ptr%range%min%z = mid%z
                        c%subcell(i,j,k)%ptr%range%max%z = c%range%max%z
                    end if

                    if (belongs(p_old, c%subcell(i,j,k)%ptr)) then
                        c%subcell(i,j,k)%ptr%type = 1
                        c%subcell(i,j,k)%ptr%part = p_old
                        c%subcell(i,j,k)%ptr%pos  = c%pos
                    else
                        c%subcell(i,j,k)%ptr%type = 0
                    end if
                end do
            end do
        end do
    end subroutine create_subcells

    !=========================================================
    ! calculate_masses
    !
    ! Recursively computes mass and center of mass of cells.
    !=========================================================
    recursive subroutine calculate_masses(c)
        type(cell), pointer :: c
        integer :: i, j, k
        real(dp) :: m_old

        c%mass = 0.0_dp
        c%c_o_m = point3d(0.0_dp, 0.0_dp, 0.0_dp)

        select case (c%type)
        case (1)
            c%mass = particles(c%pos)%m
            c%c_o_m = particles(c%pos)%p
        case (2)
            do i = 1,2
                do j = 1,2
                    do k = 1,2
                        if (associated(c%subcell(i,j,k)%ptr)) then
                            call calculate_masses(c%subcell(i,j,k)%ptr)
                            m_old = c%mass
                            c%mass = c%mass + c%subcell(i,j,k)%ptr%mass
                            if (c%mass > 0.0_dp) then
                                c%c_o_m%x = (m_old*c%c_o_m%x + &
                                    c%subcell(i,j,k)%ptr%mass * &
                                    c%subcell(i,j,k)%ptr%c_o_m%x) / c%mass

                                c%c_o_m%y = (m_old*c%c_o_m%y + &
                                    c%subcell(i,j,k)%ptr%mass * &
                                    c%subcell(i,j,k)%ptr%c_o_m%y) / c%mass

                                c%c_o_m%z = (m_old*c%c_o_m%z + &
                                    c%subcell(i,j,k)%ptr%mass * &
                                    c%subcell(i,j,k)%ptr%c_o_m%z) / c%mass
                            end if
                        end if
                    end do
                end do
            end do
        end select
    end subroutine calculate_masses

    !=========================================================
    ! calculate_forces_aux
    !
    ! Recursively computes force contribution of a cell on
    ! particle idx.
    !=========================================================
    recursive subroutine calculate_forces_aux(idx, c)
        integer, intent(in) :: idx
        type(cell), pointer :: c
        integer :: i, j, k
        type(vector3d) :: rji
        real(dp) :: r2, r3, l, d

        select case (c%type)
        case (1)
            if (idx /= c%pos) then
                rji = c%c_o_m - particles(idx)%p
                r2 = rji%x**2 + rji%y**2 + rji%z**2
                r3 = r2 * sqrt(r2)
                acc(idx) = acc(idx) + (particles(c%pos)%m / r3) * rji
            end if

        case (2)
            l = c%range%max%x - c%range%min%x
            rji = c%c_o_m - particles(idx)%p
            r2 = rji%x**2 + rji%y**2 + rji%z**2
            d  = sqrt(r2)

            if (l/d < theta) then
                r3 = r2 * d
                acc(idx) = acc(idx) + (c%mass / r3) * rji
            else
                do i = 1,2
                    do j = 1,2
                        do k = 1,2
                            if (associated(c%subcell(i,j,k)%ptr)) then
                                call calculate_forces_aux(idx, c%subcell(i,j,k)%ptr)
                            end if
                        end do
                    end do
                end do
            end if
        end select
    end subroutine calculate_forces_aux

    !=========================================================
    ! calculate_forces
    !
    ! Computes accelerations for all particles.
    ! This loop is parallelised with OpenMP.
    !=========================================================
    subroutine calculate_forces(root)
        type(cell), pointer :: root
        integer :: i

        acc = vector3d(0.0_dp, 0.0_dp, 0.0_dp)

        !$omp parallel do default(shared) private(i)
        do i = 1, n
            call calculate_forces_aux(i, root)
        end do
        !$omp end parallel do
    end subroutine calculate_forces

end module barnes_hut
