!!$     This program is free software: you can redistribute it and/or modify it under the terms of the
!!$     GNU General Public License as published by the Free Software Foundation, either version 3 of
!!$     the License, or (at your option) any later version.
!!$
!!$     This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!!$     without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!!$     See the GNU General Public License for more details.
!!$
!!$     You should have received a copy of the GNU General Public License along with this program.
!!$     If not, see <https://www.gnu.org/licenses/>.
module remeshsimplex
!> @brief type definition for the original mesh structure.
!> this type is used to represent the original mesh in the simulation.
!> it contains all the necessary information about the mesh elements,
!> nodes, and connectivity required for the remeshing process.
   type :: originalmeshtype
!*** nr. nodes:
      integer::nno
!*** nr. edges
      integer::nar
!*** conectivities
      integer::nel
      integer, dimension(:), allocatable::elno, elni
!*** edge data
      integer, dimension(:), allocatable::elir, elar, arni, arno, noir, noar
!*** list of marked edges defined from nodes
      integer::nmarkednodepairs
      integer, dimension(:, :), allocatable::markednodepairs
   end type originalmeshtype
!> @brief type definition for a renewed mesh structure.
!>
!> this type is used to represent a mesh that has been renewed or updated.
!> it contains all the necessary information and data structures to describe
!> the new state of the mesh after a remeshing operation.
   type :: renewedmeshtype
      integer::nno, npoint, nbar2, ntria3, ntetra4
      integer::nel
      integer, dimension(:), allocatable::elno, elni
      integer, dimension(:, :), allocatable::parentnodes
   end type renewedmeshtype
contains
!> @brief splits the given mesh into a new mesh.
!>
!> this subroutine takes an existing mesh and splits it into a new mesh.
!> the details of how the mesh is split are determined by the implementation
!> within this subroutine.
!>
!> @param[in] mesh the original mesh to be split.
!> @param[out] newmesh the resulting new mesh after splitting.
   subroutine splitmesh(mesh, newmesh)
      USE OMP_LIB
      real(8)::starttime, endtime
      type(originalmeshtype)::mesh
      type(renewedmeshtype)::newmesh
      integer, dimension(:), allocatable::mark
!*** determine edge relation indices for the original mesh
      starttime = omp_get_wtime()
      call edgerelations(mesh)
      call edgerelationsmarked(mesh, mark)
!*** create edge nodes
      call createedgenodes(mesh, newmesh, mark)
!*** nummber of elements by type
      newmesh%npoint = 0
      newmesh%nbar2 = 0
      newmesh%ntria3 = 0
      newmesh%ntetra4 = 0
!*** dimensions the point elements
      call splitpoint0(mesh, newmesh%npoint)
!*** dimensions the split of 2-node elements (bar2)
      call splitbar20(mark, mesh, newmesh%nbar2)
!*** dimensions the split of 3-node elements (tria3)
      call splittria30(mark, mesh, newmesh%ntria3)
!*** dimensions the split of 4-node elements (tetra4)
      call splittetra40(mark, mesh, newmesh%ntetra4)
!*** now dimensions the essentials of the new mesh
      call dimensionsnewmesh(newmesh)
!*** point
      call splitpoint1(mesh, newmesh)
!*** bar2
      call splitbar21(mark, mesh, newmesh)
!*** tria3
      call splittria31(mark, mesh, newmesh)
!*** tetra4
      call splittetra41(mark, mesh, newmesh)
      endtime = omp_get_wtime()
      print *, "splitting took", endtime - starttime, "seconds"
   end subroutine splitmesh
!> @brief splits a point in the mesh based on the given mark.
!>
!> this subroutine takes a mark and a mesh, and generates a new mesh by splitting a point.
!>
!> @param mark an integer or logical value indicating the point to be split.
!> @param mesh the original mesh data structure.
!> @param newmesh the new mesh data structure after the point has been split.
!> this subroutine is responsible for splitting a point in the mesh
!> according to the specified mark. it updates the mesh structure
!> and the number of points accordingly.
!>
!> @param[in] mark an integer value representing the mark used to determine the split.
!> @param[inout] mesh a derived type or array representing the mesh structure that will be updated.
!> @param[inout] npoint an integer representing the number of points in the mesh, which will be updated.
   subroutine splitpoint0(mesh, npoint)
      type(originalmeshtype)::mesh
      npoint = 0
      do iel = 1, mesh%nel
         if (mesh%elni(iel + 1) - mesh%elni(iel) .eq. 1) then
            npoint = npoint + 1
         end if
      end do
      write (*, *) "npoint before=", npoint
      write (*, *) "npoint after=", npoint
   end subroutine splitpoint0
!> @brief splits a point in the mesh based on the provided mark.
!>
!> this subroutine takes a mark and the current mesh, and generates a new mesh
!> by splitting a point according to the mark. it is used in mesh refinement
!> processes.
!>
!> @param[in] mark an integer or logical value indicating the point to be split.
!> @param[in] mesh the original mesh data structure.
!> @param[out] newmesh the new mesh data structure after the point has been split.
   subroutine splitpoint1(mesh, newmesh)
      type(originalmeshtype)::mesh
      type(renewedmeshtype)::newmesh
      ipoint = 0
      last = 0
      do iel = 1, mesh%nel
         if (mesh%elni(iel + 1) - mesh%elni(iel) .eq. 1) then
            ipoint = ipoint + 1
            ing = mesh%elno(mesh%elni(iel))! node 1 of edge
            newmesh%elno(newmesh%elni(last + ipoint):newmesh%elni(last + 1 + ipoint) - 1) = [ing, iel]
         end if
      end do
   end subroutine splitpoint1
!> @brief subroutine to split a bar element in a mesh.
!>
!> this subroutine is responsible for splitting a bar element in the mesh
!> based on the provided mark and mesh data. the number of bar elements
!> after the split is updated in `nbar2`.
!>
!> @param mark an integer array indicating the elements to be split.
!> @param mesh a derived type or array containing the mesh data.
!> @param nbar2 an integer representing the number of bar elements after the split.
   subroutine splitbar20(mark, mesh, nbar2)
      type(originalmeshtype)::mesh
      integer, dimension(*)::mark
      nbar2 = 0
      do iel = 1, mesh%nel
         if (mesh%elni(iel + 1) - mesh%elni(iel) .eq. 2) then
            nbar2 = nbar2 + 1
         end if
      end do
      write (*, *) "nbar2 before=", nbar2
      nbar2 = 0
      do iel = 1, mesh%nel
         if (mesh%elni(iel + 1) - mesh%elni(iel) .eq. 2) then
            ia = mesh%elar(mesh%elir(iel))
            if (mark(ia) .ne. 0) then
               nbar2 = nbar2 + 2
            else
               nbar2 = nbar2 + 1
            end if
         end if
      end do
      write (*, *) "nbar2 after=", nbar2
   end subroutine splitbar20
!> @brief splits a bar element into two new elements.
!>
!> this subroutine takes a marked element from the mesh and splits it into two new elements,
!> updating the mesh structure accordingly.
!>
!> @param mark an integer array indicating which elements are marked for splitting.
!> @param mesh the original mesh structure containing the elements to be split.
!> @param newmesh the updated mesh structure after the elements have been split.
   subroutine splitbar21(mark, mesh, newmesh)
      type(originalmeshtype)::mesh
      type(renewedmeshtype)::newmesh
      integer, dimension(*)::mark
      integer::ibar2, ing1, ing2, ing3, ia, iel, last
      ibar2 = 0
      last = newmesh%npoint
      do iel = 1, mesh%nel
         if (mesh%elni(iel + 1) - mesh%elni(iel) .eq. 2) then
            ia = mesh%elar(mesh%elir(iel))
            ing2 = mark(ia)
            select case (ing2)
            case (0)
               ing1 = mesh%elno(mesh%elni(iel))! node 1 of edge
               ing3 = mesh%elno(mesh%elni(iel) + 1)! node 2 of edge
               ibar2 = ibar2 + 1
               newmesh%elno(newmesh%elni(last + ibar2):newmesh%elni(last + 1 + ibar2) - 1) = [ing1, ing3, iel]
            case default
               ing1 = mesh%elno(mesh%elni(iel))! node 1 of edge
               ing3 = mesh%elno(mesh%elni(iel) + 1)! node 2 of edge
               ibar2 = ibar2 + 1
               newmesh%elno(newmesh%elni(last + ibar2):newmesh%elni(last + 1 + ibar2) - 1) = [ing1, ing2, iel]
               ibar2 = ibar2 + 1
               newmesh%elno(newmesh%elni(last + ibar2):newmesh%elni(last + 1 + ibar2) - 1) = [ing2, ing3, iel]
            end select
         end if
      end do
   end subroutine splitbar21
!> @brief splits a triangular element into smaller elements.
!>
!> this subroutine is responsible for splitting a triangular element
!> into smaller elements based on certain criteria. it modifies the
!> mesh structure and updates the number of triangular elements.
!>
!> @param mark an integer array that marks the elements to be split.
!> @param mesh a derived type or array representing the mesh structure.
!> @param ntria3 an integer representing the number of triangular elements.
!>
!> @note ensure that the mesh structure is properly initialized before
!>       calling this subroutine.
   subroutine splittria30(mark, mesh, ntria3)
      type(originalmeshtype)::mesh
      integer, dimension(*) ::mark
      ntria3 = 0
      do iel = 1, mesh%nel
         if (mesh%elni(iel + 1) - mesh%elni(iel) .eq. 3) then
            ntria3 = ntria3 + 1
         end if
      end do
      write (*, *) "ntria3 before=", ntria3
      ntria3 = 0
      do iel = 1, mesh%nel
         if (mesh%elni(iel + 1) - mesh%elni(iel) .eq. 3) then
            kk = 1
            do ik = mesh%elir(iel), mesh%elir(iel + 1) - 1
               if (mark(mesh%elar(ik)) .ne. 0) kk = kk + 1
            end do
            ntria3 = ntria3 + kk
         end if
      end do
      write (*, *) "ntria3 after=", ntria3
   end subroutine splittria30
!> @brief splits a triangular element into three new elements.
!>
!> this subroutine takes a marked triangular element and splits it into
!> new triangular elements. the new elements are added to the mesh.
!>
!> @param mark an integer indicating the marked triangular element to be split.
!> @param mesh the original mesh containing the triangular elements.
!> @param newmesh the new mesh where the split triangular elements will be added.
   subroutine splittria31(mark, mesh, newmesh)
      implicit none
      type(originalmeshtype), intent(in) :: mesh
      type(renewedmeshtype), intent(inout) :: newmesh
      integer, dimension(:), intent(in) :: mark
      integer :: elem, markededgecount, edge, localedge, node1, node2, node3
      integer :: newnode1, newnode2, newnode3, newnode4, newnode5, newnode6, tricounter
      integer :: edgenode1, edgenode2, prevnode1, prevnode2, prevnode3, mark1, mark2, newmeshstart
      tricounter = 0  ! initialize triangle counter
      newmeshstart = newmesh%npoint + newmesh%nbar2
! loop through each element in the mesh
      do elem = 1, mesh%nel
! check if the element has 3 nodes (triangular)
         if (mesh%elni(elem + 1) - mesh%elni(elem) == 3) then
            markededgecount = 0
! count the number of marked edges
            do edge = mesh%elir(elem), mesh%elir(elem + 1) - 1
               if (mark(mesh%elar(edge)) /= 0) markededgecount = markededgecount + 1
            end do
            select case (markededgecount)
            case (0)
! no marked edges, retain original triangle
               newnode1 = mesh%elno(mesh%elni(elem))
               newnode2 = mesh%elno(mesh%elni(elem) + 1)
               newnode3 = mesh%elno(mesh%elni(elem) + 2)
               tricounter = tricounter + 1
             newmesh%elno(newmesh%elni(tricounter + newmeshstart):newmesh%elni(tricounter + newmeshstart + 1) - 1) = [newnode1, newnode2, newnode3,elem]
            case (1)
! one marked edge, split triangle
               do edge = mesh%elir(elem), mesh%elir(elem + 1) - 1
                  if (mark(mesh%elar(edge)) /= 0) then
                     localedge = localedgeinelement(mesh, elem, mesh%elar(edge))
                     node1 = localedge
                     node2 = iclock(3, node1 + 1)
                     node3 = iclock(3, node1 + 2)
                     newnode1 = mesh%elno(mesh%elni(elem) - 1 + node1)
                     newnode2 = mark(mesh%elar(edge))
                     newnode3 = mesh%elno(mesh%elni(elem) - 1 + node2)
                     newnode4 = mesh%elno(mesh%elni(elem) - 1 + node3)
! create two new triangles
                     tricounter = tricounter + 1
                   newmesh%elno(newmesh%elni(tricounter + newmeshstart):newmesh%elni(tricounter + newmeshstart + 1) - 1) = [newnode1, newnode2, newnode4,elem]
                     tricounter = tricounter + 1
                   newmesh%elno(newmesh%elni(tricounter + newmeshstart):newmesh%elni(tricounter + newmeshstart + 1) - 1) = [newnode2, newnode3, newnode4,elem]
                  end if
               end do
            case (2)
               do prevnode1 = 1, 3
                  prevnode2 = iclock(3, prevnode1 + 1)
                  prevnode3 = iclock(3, prevnode1 + 2)
                  edgenode1 = mesh%elar(mesh%elir(elem) - 1 + prevnode1)
                  edgenode2 = mesh%elar(mesh%elir(elem) - 1 + prevnode2)
                  mark1 = mark(edgenode1)
                  mark2 = mark(edgenode2)
                  if (mark1 /= 0 .and. mark2 /= 0) then
                     newnode1 = mesh%elno(mesh%elni(elem) - 1 + prevnode1)!node1
                     newnode2 = mark(edgenode1)
                     newnode3 = mesh%elno(mesh%elni(elem) - 1 + prevnode2)
                     newnode4 = mark(edgenode2)
                     newnode5 = mesh%elno(mesh%elni(elem) - 1 + prevnode3)
                     if (newnode2 > newnode4) then
                        tricounter = tricounter + 1
                      newmesh%elno(newmesh%elni(tricounter + newmeshstart):newmesh%elni(tricounter + newmeshstart + 1) - 1) = [newnode1, newnode2, newnode5,elem]
                        tricounter = tricounter + 1
                      newmesh%elno(newmesh%elni(tricounter + newmeshstart):newmesh%elni(tricounter + newmeshstart + 1) - 1) = [newnode2, newnode4, newnode5,elem]
                     else
                        tricounter = tricounter + 1
                      newmesh%elno(newmesh%elni(tricounter + newmeshstart):newmesh%elni(tricounter + newmeshstart + 1) - 1) = [newnode1, newnode2, newnode4,elem]
                        tricounter = tricounter + 1
                      newmesh%elno(newmesh%elni(tricounter + newmeshstart):newmesh%elni(tricounter + newmeshstart + 1) - 1) = [newnode1, newnode4, newnode5,elem]
                     end if
                     tricounter = tricounter + 1
                   newmesh%elno(newmesh%elni(tricounter + newmeshstart):newmesh%elni(tricounter + newmeshstart + 1) - 1) = [newnode2, newnode3, newnode4,elem]                 
                  end if
               end do
            case (3)
! three marked edges, fully subsplit triangle
               newnode1 = mesh%elno(mesh%elni(elem))
               newnode2 = mark(mesh%elar(mesh%elir(elem)))
               newnode3 = mesh%elno(mesh%elni(elem) + 1)
               newnode4 = mark(mesh%elar(mesh%elir(elem) + 1))
               newnode5 = mesh%elno(mesh%elni(elem) + 2)
               newnode6 = mark(mesh%elar(mesh%elir(elem) + 2))
! create four new triangles
               tricounter = tricounter + 1
             newmesh%elno(newmesh%elni(tricounter + newmeshstart):newmesh%elni(tricounter + newmeshstart + 1) - 1) = [newnode1, newnode2, newnode6,elem]
               tricounter = tricounter + 1
             newmesh%elno(newmesh%elni(tricounter + newmeshstart):newmesh%elni(tricounter + newmeshstart + 1) - 1) = [newnode2, newnode3, newnode4,elem]
               tricounter = tricounter + 1
             newmesh%elno(newmesh%elni(tricounter + newmeshstart):newmesh%elni(tricounter + newmeshstart + 1) - 1) = [newnode4, newnode5, newnode6,elem]
               tricounter = tricounter + 1
             newmesh%elno(newmesh%elni(tricounter + newmeshstart):newmesh%elni(tricounter + newmeshstart + 1) - 1) = [newnode2, newnode4, newnode6,elem]
            end select
         end if
      end do
   end subroutine splittria31
!> @brief splits tetrahedral elements in a mesh.
!>
!> this subroutine is responsible for splitting tetrahedral elements in a mesh
!> based on certain criteria. it modifies the mesh structure by dividing the
!> tetrahedral elements.
!>
!> @param mark an array indicating which tetrahedral elements need to be split.
!> @param mesh the mesh data structure containing the tetrahedral elements.
!> @param ntetra4 the number of tetrahedral elements in the mesh.
   subroutine splittetra40(mark, mesh, ntetra4)
      type(originalmeshtype)::mesh
      integer, dimension(*)::mark
      integer, dimension(6)::midnodes
!*** determines the new number of elements
      ntetra4 = 0
      do iel = 1, mesh%nel
         if (mesh%elni(iel + 1) - mesh%elni(iel) .eq. 4) then
            ntetra4 = ntetra4 + 1
         end if
      end do
      write (*, *) "ntetra4 before=", ntetra4
      ntetra4 = 0
      do iel = 1, mesh%nel
         if (mesh%elni(iel + 1) - mesh%elni(iel) .eq. 4) then
            midnodes = 0
            do ik = mesh%elir(iel), mesh%elir(iel + 1) - 1
               iar = mesh%elar(ik)
               midnodes(ik - mesh%elir(iel) + 1) = mark(iar)
            end do
            call splittetwork0(midnodes, ntets)
            ntetra4 = ntetra4 + ntets
         end if
      end do
      write (*, *) "ntetra4 after=", ntetra4
   end subroutine splittetra40
!> splits a tetrahedron into smaller tetrahedra based on the provided
!> marking and mesh information.
!>
!> this subroutine takes an existing mesh and a marking array to
!> generate a new refined mesh by splitting the tetrahedra.
!>
!> @param mark an array indicating which tetrahedra to split.
!> @param mesh the original mesh data structure.
!> @param newmesh the new mesh data structure after splitting.
   subroutine splittetra41(mark, mesh, newmesh)
      type(originalmeshtype)::mesh
      type(renewedmeshtype)::newmesh
      integer, dimension(:) :: mark
      integer, dimension(4, 8)::tets
      integer, dimension(6)::midnodes
      itetra4 = 0
      last = newmesh%npoint + newmesh%nbar2 + newmesh%ntria3
      do iel = 1, mesh%nel
         if (mesh%elni(iel + 1) - mesh%elni(iel) .eq. 4) then
            midnodes = 0
            do ik = mesh%elir(iel), mesh%elir(iel + 1) - 1
               iar = mesh%elar(ik)
               midnodes(ik - mesh%elir(iel) + 1) = mark(iar)
            end do
            call splittetwork1(mesh%elno(mesh%elni(iel):mesh%elni(iel + 1) - 1), midnodes, ntets, tets)
            do i = 1, ntets
               itetra4 = itetra4 + 1
        newmesh%elno(newmesh%elni(itetra4+last):newmesh%elni(itetra4+last+1)-1) = [tets(1,i),tets(2,i),tets(3,i),tets(4,i),iel]!iel]
            end do
         end if
      end do
   end subroutine splittetra41
!> @brief splits tetrahedral elements in the mesh.
!>
!> this subroutine is responsible for splitting the work tetrahedral elements
!> of the mesh. it takes the current nodes, mid-nodes, and the number of tetrahedral
!> elements as input parameters and performs the necessary operations to split them.
!>
!> @param nodes an array containing the coordinates of the nodes in the mesh.
!> @param midnodes an array containing the coordinates of the mid-nodes in the mesh.
!> @param ntets the number of tetrahedral elements in the mesh.
   subroutine splittetwork0(midnodes, ntets)
      integer, dimension(6)::midnodes, marked
      integer, dimension(3, 4)::nodeedges
      integer, dimension(5)::imarks
!*** edge description (constants):
      nodeedges(1, 1) = 1!
      nodeedges(1, 2) = 1!
!--------------------
      nodeedges(2, 2) = 2!
      nodeedges(1, 3) = 2!
!--------------------
      nodeedges(2, 1) = 3!
      nodeedges(2, 3) = 3!
!--------------------
      nodeedges(3, 3) = 4!
      nodeedges(1, 4) = 4!
!--------------------
      nodeedges(2, 4) = 5!
      nodeedges(3, 1) = 5!
!--------------------
      nodeedges(3, 2) = 6!
      nodeedges(3, 4) = 6!
      nmarked = 0
      do i = 1, 6
         if (midnodes(i) .gt. 0) then
            nmarked = nmarked + 1
            marked(nmarked) = i
         end if
      end do
!*** number of marks per node
!*** and local edges
      imarks = 1
      do i = 1, 4
         imarks(i + 1) = imarks(i)
         do j = 1, 3
            ie = nodeedges(j, i)! all elementwise edges of each node
            if (midnodes(ie) .gt. 0) then! if this edge is marked
               imarks(i + 1) = imarks(i + 1) + 1! updates the pointer
            end if
         end do
      end do
!*** now the division
!*** according to the cases
      select case (nmarked)
      case (0)
         ntets = 1
      case (1)
         ntets = 2
      case (2)
         i1 = 0
         i2 = 0
         mm = 0
         do i = 1, 4
            mm = max(mm, imarks(i + 1) - imarks(i))
         end do
!*** one node with two marked edges:
         if (mm .eq. 2) then
            ntets = 3
         else
            ntets = 4
         end if
      case (3)
!*** now with 3 marked edges
!*** find case
         n4 = 0
         n3 = 0
         do i = 1, 4
            if (imarks(i + 1) - imarks(i) .eq. 3) then
               n4 = i
            else if (imarks(i + 1) - imarks(i) .eq. 0) then
               n3 = i
            end if
         end do
         if (n3 .ne. 0) then
            ntets = 4
         else if (n4 .ne. 0) then
            ntets = 4
         else
            ntets = 5
         end if
      case (4)
         ntets = 6
      case (5)
         ntets = 7
      case (6)
         ntets = 8
      end select
   end subroutine splittetwork0
!> splits the network of tetrahedral elements.
!>
!> this subroutine takes the existing nodes and midnodes of the tetrahedral
!> elements and splits them into smaller tetrahedra.
!>
!> @param nodes an array containing the coordinates of the original nodes.
!> @param midnodes an array containing the coordinates of the midpoints of the edges.
!> @param ntets the number of tetrahedral elements.
!> @param tets an array containing the indices of the nodes that form each tetrahedral element.
   subroutine splittetwork1(nodes, midnodes, ntets, tets)
      implicit none
! inputs
      integer, dimension(4) :: nodes          ! node indices of the tetrahedron's vertices
      integer, dimension(6) :: midnodes       ! midpoint node indices for each edge (0 if not present)
! outputs
      integer :: ntets                        ! number of tetrahedra after division
      integer, dimension(4, *) :: tets        ! node indices of the resulting tetrahedra
! internal variables
      integer, dimension(6) :: marked         ! indices of marked edges (edges with midnodes)
      integer, dimension(6) :: post, prev     ! next and previous nodes for each edge
      integer, dimension(2, 6) :: edgenodes    ! node indices for each edge
      integer, dimension(3, 4) :: nodeedges    ! edge indices connected to each node
      integer, dimension(5) :: imarks         ! start indices in jmarks for each node
      integer, dimension(12) :: jmarks        ! marked edges per node
      integer :: nmarked                      ! number of marked edges
      integer :: i, j, ie, mm, ie1, itemp
      integer :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10
! define the edges of the tetrahedron with node pairs
! edgenodes(1:2, 1:6): local node indices for each edge
! nodeedges(1:3, 1:4): edge indices connected to each node
! prev(1:6): node 1 of opposing edge
! post(1:6): node 2 of opposing edge
! edge 1: nodes 1 and 2
      edgenodes(1, 1) = 1
      edgenodes(2, 1) = 2
      nodeedges(1, 1) = 1   ! node 1 connected to edge 1
      nodeedges(1, 2) = 1   ! node 2 connected to edge 1
      prev(1) = 3          ! opposing edge node 1
      post(1) = 4          ! opposing edge node 2
! edge 2: nodes 2 and 3
      edgenodes(1, 2) = 2
      edgenodes(2, 2) = 3
      nodeedges(2, 2) = 2   ! node 2 connected to edge 2
      nodeedges(1, 3) = 2   ! node 3 connected to edge 2
      prev(2) = 1          ! opposing edge node 1
      post(2) = 4          ! opposing edge node 2
! edge 3: nodes 1 and 3
      edgenodes(1, 3) = 1
      edgenodes(2, 3) = 3
      nodeedges(2, 1) = 3   ! node 1 connected to edge 3
      nodeedges(2, 3) = 3   ! node 3 connected to edge 3
      prev(3) = 4          ! opposing edge node 1
      post(3) = 2          ! opposing edge node 2
! edge 4: nodes 3 and 4
      edgenodes(1, 4) = 3
      edgenodes(2, 4) = 4
      nodeedges(3, 3) = 4   ! node 3 connected to edge 4
      nodeedges(1, 4) = 4   ! node 4 connected to edge 4
      prev(4) = 1          ! opposing edge node 1
      post(4) = 2          ! opposing edge node 2
! edge 5: nodes 4 and 1
      edgenodes(1, 5) = 4
      edgenodes(2, 5) = 1
      nodeedges(2, 4) = 5   ! node 4 connected to edge 5
      nodeedges(3, 1) = 5   ! node 1 connected to edge 5
      prev(5) = 3          ! opposing edge node 1
      post(5) = 2          ! opposing edge node 2
! edge 6: nodes 2 and 4
      edgenodes(1, 6) = 2
      edgenodes(2, 6) = 4
      nodeedges(3, 2) = 6   ! node 2 connected to edge 6
      nodeedges(3, 4) = 6   ! node 4 connected to edge 6
      prev(6) = 3          ! opposing edge node 1
      post(6) = 1          ! opposing edge node 2
! identify the marked edges (edges with midnodes)
      nmarked = 0
      do i = 1, 6
         if (midnodes(i) .gt. 0) then
            nmarked = nmarked + 1
            marked(nmarked) = i
         end if
      end do
! determine the number of marked edges per node and record the edges
      imarks(1) = 1
      jmarks = 0
      do i = 1, 4
         imarks(i + 1) = imarks(i)
         do j = 1, 3
            ie = nodeedges(j, i)            ! edge connected to node i
            if (midnodes(ie) .gt. 0) then   ! if edge is marked
               jmarks(imarks(i + 1)) = ie    ! record edge index
               imarks(i + 1) = imarks(i + 1) + 1 ! update pointer
            end if
         end do
      end do
! split the tetrahedron based on the number of marked edges
      select case (nmarked)
      case (0)
! case #1: no marked edges (original tetrahedron)
         ntets = 1
         tets(1:4, 1) = nodes(1:4)
      case (1)
! case #2: one marked edge
! identify the nodes connected by the marked edge
         n1 = edgenodes(1, marked(1))
         n2 = edgenodes(2, marked(1))
! find the other two nodes (n3 and n4)
! of opposing edge
         call nodesfrom2nodes(edgenodes, post, prev, n3, n1, n2, n4)
! convert local node indices to global node indices
         n1 = nodes(n1)
         n2 = nodes(n2)
         n3 = nodes(n3)
         n4 = nodes(n4)
         n5 = midnodes(marked(1)) ! midpoint node of the marked edge
! build two new tetrahedra
         ntets = 2
         call splitetintoitself(n3, n4, n5, n2, tets(1:4, 1))
         call splitetintoitself(n4, n3, n5, n1, tets(1:4, 2))
      case (2)
! cases #3 and #4: two marked edges
! find the maximum number of marked edges connected to any node
         mm = 0
         do i = 1, 4
            mm = max(mm, imarks(i + 1) - imarks(i))
         end do
         if (mm == 2) then
! case #3: one node with two marked edges
! identify the node with 0 marked edges (n3) and the node with 2 marked edges (n4)
            n3 = 0
            n4 = 0
            do i = 1, 4
               if (imarks(i + 1) - imarks(i) == 0) then
                  n3 = i
               else if (imarks(i + 1) - imarks(i) == 2) then
                  n4 = i
               end if
            end do
            if (n3 == 0 .or. n4 == 0) then
               stop "error: either n3 or n4 is zero in case #3"
            end if
! get the remaining nodes (n1 and n2)
            call nodesfrom2nodes(edgenodes, post, prev, n1, n3, n4, n2)
! convert local node indices to global node indices
! midpoint nodes of the marked edges
            n5 = splitnode(n2, n4, midnodes, edgenodes)
            n6 = splitnode(n1, n4, midnodes, edgenodes)
            n1 = nodes(n1)
            n2 = nodes(n2)
            n3 = nodes(n3)
            n4 = nodes(n4)
! build three new tetrahedra
            ntets = 3
            call splitpyramidintotets(n1, n2, n5, n6, n3, tets(1:4, 1), tets(1:4, 2))
            call splitetintoitself(n6, n5, n3, n4, tets(1:4, 3))
         else
! case #4: each node has one marked edge
! identify two nodes connected by a marked edge
            ie1 = marked(1)
            n1 = edgenodes(1, ie1)
            n2 = edgenodes(2, ie1)
! find the other two nodes (n3 and n4)
            call nodesfrom2nodes(edgenodes, post, prev, n3, n1, n2, n4)
! midpoint nodes of the marked edges
            n5 = splitnode(n1, n2, midnodes, edgenodes)
            n6 = splitnode(n3, n4, midnodes, edgenodes)
! convert local node indices to global node indices
            n1 = nodes(n1)
            n2 = nodes(n2)
            n3 = nodes(n3)
            n4 = nodes(n4)
! build four new tetrahedra
            ntets = 4
            call splitetintoitself(n2, n3, n5, n6, tets(1:4, 1))
            call splitetintoitself(n2, n6, n5, n4, tets(1:4, 2))
            call splitetintoitself(n1, n5, n6, n4, tets(1:4, 3))
            call splitetintoitself(n1, n5, n3, n6, tets(1:4, 4))
         end if
      case (3)
! cases #5, #6, and #7: three marked edges
! identify nodes based on the number of marked edges connected
         n4 = 0
         n3 = 0
         do i = 1, 4
            if (imarks(i + 1) - imarks(i) == 3) then
               n4 = i
            else if (imarks(i + 1) - imarks(i) == 0) then
               n3 = i
            end if
         end do
         if (n4 .ne. 0) then
! case #6: one node with three marked edges
            n1 = iclock(4, n4 + 1)   ! next node in cyclic order
            call nodesfrom2nodes(edgenodes, post, prev, n2, n1, n4, n3)
! midpoint nodes of the marked edges
            n5 = splitnode(n3, n4, midnodes, edgenodes)
            n6 = splitnode(n1, n4, midnodes, edgenodes)
            n7 = splitnode(n2, n4, midnodes, edgenodes)
! convert local node indices to global node indices
            n1 = nodes(n1)
            n2 = nodes(n2)
            n3 = nodes(n3)
            n4 = nodes(n4)
! build four new tetrahedra
            ntets = 4
            call splitetintoitself(n4, n7, n6, n5, tets(1:4, 1))
            call splitprismintotets(n2, n7, n6, n1, n3, n5, tets(1:4, 2), tets(1:4, 3), tets(1:4, 4))
         else if (n3 .ne. 0) then
! case #5: one node with zero marked edges
            n4 = iclock(4, n3 + 1)
            call nodesfrom2nodes(edgenodes, post, prev, n1, n3, n4, n2)
! midpoint nodes of the marked edges
            n5 = splitnode(n1, n4, midnodes, edgenodes)
            n6 = splitnode(n1, n2, midnodes, edgenodes)
            n7 = splitnode(n2, n4, midnodes, edgenodes)
! convert local node indices to global node indices
            n1 = nodes(n1)
            n2 = nodes(n2)
            n3 = nodes(n3)
            n4 = nodes(n4)
! build four new tetrahedra
            ntets = 4
            call splitetintoitself(n1, n6, n3, n5, tets(1:4, 1))
            call splitetintoitself(n2, n3, n6, n7, tets(1:4, 2))
            call splitetintoitself(n6, n3, n5, n7, tets(1:4, 3))
            call splitetintoitself(n5, n7, n3, n4, tets(1:4, 4))
         else
! case #7: three nodes with one marked edge each
            n1 = 0
            n2 = 0
            do i = 1, 4
               if (imarks(i + 1) - imarks(i) == 1) then
                  ie = jmarks(imarks(i))
                  if (n1 == 0) then
                     n1 = i
                     if (edgenodes(1, ie) .eq. n1) then
                        n4 = edgenodes(2, ie)
                     else
                        n4 = edgenodes(1, ie)
                     end if
                  else
                     n2 = i
                     if (edgenodes(1, ie) .eq. n2) then
                        n3 = edgenodes(2, ie)
                     else
                        n3 = edgenodes(1, ie)
                     end if
                  end if
               end if
            end do
! ensure correct orientation
            if (tetra4symmetry(n1, n2, n3, n4) == -1) then
               itemp = n3
               n3 = n4
               n4 = itemp
               itemp = n1
               n1 = n2
               n2 = itemp
            end if
! midpoint nodes of the marked edges
            n5 = splitnode(n1, n4, midnodes, edgenodes)
            n6 = splitnode(n2, n3, midnodes, edgenodes)
            n7 = splitnode(n3, n4, midnodes, edgenodes)
! convert local node indices to global node indices
            n1 = nodes(n1)
            n2 = nodes(n2)
            n3 = nodes(n3)
            n4 = nodes(n4)
! build five new tetrahedra
            ntets = 5
            call splitpyramidintotets(n6, n7, n4, n2, n5, tets(1:4, 1), tets(1:4, 2))
            call splitpyramidintotets(n3, n1, n5, n7, n6, tets(1:4, 3), tets(1:4, 4))
            call splitetintoitself(n1, n2, n6, n5, tets(1:4, 5))
         end if
      case (4)
! cases #8 and #9: four marked edges
! find the maximum number of marked edges connected to any node
         mm = 0
         do i = 1, 4
            mm = max(mm, imarks(i + 1) - imarks(i))
         end do
         if (mm == 3) then
! case #8: one node with three marked edges
            n3 = 0
            n4 = 0
            do i = 1, 4
               if (imarks(i + 1) - imarks(i) == 1) then
                  n3 = i
               else if (imarks(i + 1) - imarks(i) == 3) then
                  n4 = i
               end if
            end do
            call nodesfrom2nodes(edgenodes, post, prev, n1, n3, n4, n2)
! midpoint nodes
            n5 = splitnode(n1, n2, midnodes, edgenodes)
            n6 = splitnode(n2, n4, midnodes, edgenodes)
            n7 = splitnode(n3, n4, midnodes, edgenodes)
            n8 = splitnode(n1, n4, midnodes, edgenodes)
! convert local node indices to global node indices
            n1 = nodes(n1)
            n2 = nodes(n2)
            n3 = nodes(n3)
            n4 = nodes(n4)
! build six new tetrahedra
            ntets = 6
            call splitetintoitself(n8, n6, n7, n4, tets(1:4, 1))
            call splitetintoitself(n7, n8, n5, n6, tets(1:4, 2))
            call splitpyramidintotets(n2, n3, n7, n6, n5, tets(1:4, 3), tets(1:4, 4))
            call splitpyramidintotets(n1, n8, n7, n3, n5, tets(1:4, 5), tets(1:4, 6))
         else
! case #9: edges form a quadrilateral
            ie = 0
            do i = 1, 6
               if (midnodes(i) == 0) then
                  ie = i
                  exit
               end if
            end do
            if (ie == 0) then
               stop "error: no unmarked edge found in case #9"
            end if
            n1 = edgenodes(1, ie)
            n2 = edgenodes(2, ie)
            call nodesfrom2nodes(edgenodes, post, prev, n3, n1, n2, n4)
! midpoint nodes
            n5 = splitnode(n1, n4, midnodes, edgenodes)
            n6 = splitnode(n1, n3, midnodes, edgenodes)
            n7 = splitnode(n2, n3, midnodes, edgenodes)
            n8 = splitnode(n2, n4, midnodes, edgenodes)
! convert local node indices to global node indices
            n1 = nodes(n1)
            n2 = nodes(n2)
            n3 = nodes(n3)
            n4 = nodes(n4)
! build six new tetrahedra
            ntets = 6
            call splitprismintotets(n8, n7, n3, n4, n5, n6, tets(1:4, 1), tets(1:4, 2), tets(1:4, 3))
            call splitprismintotets(n1, n2, n8, n5, n6, n7, tets(1:4, 4), tets(1:4, 5), tets(1:4, 6))
         end if
      case (5)
! case #10: five marked edges
         ie = 0
         do i = 1, 6
            if (midnodes(i) == 0) then
               ie = i
               exit
            end if
         end do
         if (ie == 0) then
            stop "error: no unmarked edge found in case #10"
         end if
         n1 = edgenodes(1, ie)
         n2 = edgenodes(2, ie)
         call nodesfrom2nodes(edgenodes, post, prev, n3, n1, n2, n4)
! midpoint nodes
         n5 = splitnode(n1, n4, midnodes, edgenodes)
         n6 = splitnode(n1, n3, midnodes, edgenodes)
         n7 = splitnode(n2, n3, midnodes, edgenodes)
         n8 = splitnode(n3, n4, midnodes, edgenodes)
         n9 = splitnode(n2, n4, midnodes, edgenodes)
! convert local node indices to global node indices
         n1 = nodes(n1)
         n2 = nodes(n2)
         n3 = nodes(n3)
         n4 = nodes(n4)
! build seven new tetrahedra
         ntets = 7
         call splitetintoitself(n9, n8, n5, n4, tets(1:4, 1))
         call splitetintoitself(n3, n6, n7, n8, tets(1:4, 2))
         call splitpyramidintotets(n6, n7, n9, n5, n8, tets(1:4, 3), tets(1:4, 4))
         call splitprismintotets(n1, n2, n9, n5, n6, n7, tets(1:4, 5), tets(1:4, 6), tets(1:4, 7))
      case (6)
! case #11: all edges are marked
! midpoint nodes of all edges
         n5 = splitnode(1, 4, midnodes, edgenodes)
         n6 = splitnode(2, 4, midnodes, edgenodes)
         n7 = splitnode(3, 4, midnodes, edgenodes)
         n8 = splitnode(1, 3, midnodes, edgenodes)
         n9 = splitnode(1, 2, midnodes, edgenodes)
         n10 = splitnode(2, 3, midnodes, edgenodes)
! convert local node indices to global node indices
         n1 = nodes(1)
         n2 = nodes(2)
         n3 = nodes(3)
         n4 = nodes(4)
! build eight new tetrahedra
         ntets = 8
         call splitetintoitself(n1, n9, n8, n5, tets(1:4, 1))
         call splitetintoitself(n2, n10, n9, n6, tets(1:4, 2))
         call splitetintoitself(n3, n8, n10, n7, tets(1:4, 3))
         call splitetintoitself(n5, n6, n7, n4, tets(1:4, 4))
         call splitoctointotets(n8, n6, n7, n5, n9, n10, tets(1:4, 5), tets(1:4, 6), tets(1:4, 7), tets(1:4, 8))
      end select
   end subroutine splittetwork1
!*** split a triangular prism in three tets
!*** chk0
!*** chk1
!*** 236 ok
!*** 1256 internal quad - always ok
!*** 3456 internal quad - always ok
!> @brief converts a prism element into three tetrahedral elements.
!>
!> this subroutine takes the indices of the vertices of a prism and
!> generates three tetrahedral elements from it.
!>
!> @param i1 index of the first vertex of the prism.
!> @param i2 index of the second vertex of the prism.
!> @param i3 index of the third vertex of the prism.
!> @param i4 index of the fourth vertex of the prism.
!> @param i5 index of the fifth vertex of the prism.
!> @param i6 index of the sixth vertex of the prism.
!> @param tet1 array to store the indices of the vertices of the first tetrahedron.
!> @param tet2 array to store the indices of the vertices of the second tetrahedron.
!> @param tet3 array to store the indices of the vertices of the third tetrahedron.
   subroutine splitprismintotets(i1, i2, i3, i4, i5, i6, tet1, tet2, tet3)
      integer, dimension(4) :: tet1, tet2, tet3
      logical :: swap46, swap13, swap16, swap25
! determine swaps
      swap46 = max(i4, i6) > max(i3, i5)
      swap13 = max(i1, i3) > max(i2, i4)
      swap16 = max(i1, i6) > max(i2, i5)
      swap25 = .not. swap16
      if (swap46) then
         if (swap13) then
            tet1 = [i1, i2, i6, i3]
            tet2 = [i1, i6, i4, i3]
            tet3 = [i1, i6, i5, i4]
         else if (swap16) then
            tet1 = [i1, i2, i6, i4]
            tet2 = [i1, i6, i5, i4]
            tet3 = [i2, i6, i4, i3]
         else
            tet1 = [i1, i2, i5, i4]
            tet2 = [i2, i3, i6, i4]
            tet3 = [i2, i6, i5, i4]
         end if
      else
         if (swap13) then
            if (swap16) then
               tet1 = [i1, i2, i6, i3]
               tet2 = [i1, i3, i5, i4]
               tet3 = [i1, i6, i5, i3]
            else
               tet1 = [i1, i3, i5, i4]
               tet2 = [i2, i5, i1, i3]
               tet3 = [i2, i6, i5, i3]
            end if
         else
            tet1 = [i1, i2, i5, i4]
            tet2 = [i2, i3, i5, i4]
            tet3 = [i2, i6, i5, i3]
         end if
      end if
   end subroutine splitprismintotets
!> \brief subroutine to perform octahedral to tetrahedral decomposition.
!>
!> this subroutine takes six input indices and decomposes an octahedron into four tetrahedra.
!>
!> \param[in] i1 index of the first vertex of the octahedron.
!> \param[in] i2 index of the second vertex of the octahedron.
!> \param[in] i3 index of the third vertex of the octahedron.
!> \param[in] i4 index of the fourth vertex of the octahedron.
!> \param[in] i5 index of the fifth vertex of the octahedron.
!> \param[in] i6 index of the sixth vertex of the octahedron.
!> \param[out] tet1 array containing the indices of the vertices of the first tetrahedron.
!> \param[out] tet2 array containing the indices of the vertices of the second tetrahedron.
!> \param[out] tet3 array containing the indices of the vertices of the third tetrahedron.
!> \param[out] tet4 array containing the indices of the vertices of the fourth tetrahedron.
   subroutine splitoctointotets(i1, i2, i3, i4, i5, i6, tet1, tet2, tet3, tet4)
      integer, dimension(4) :: tet1, tet2, tet3, tet4
      logical :: swap3546
! determine if diagonal swap between (35) and (46) is needed
      swap3546 = max(i3, i5) > max(i4, i6)
      if (swap3546) then
         tet1 = [i1, i3, i4, i5]
         tet2 = [i1, i3, i5, i6]
         tet3 = [i3, i4, i5, i2]
         tet4 = [i3, i5, i6, i2]
      else
         tet1 = [i1, i4, i6, i3]
         tet2 = [i1, i6, i4, i5]
         tet3 = [i2, i4, i6, i5]
         tet4 = [i2, i6, i4, i3]
      end if
   end subroutine splitoctointotets
!---------------------------------
!*** 34 ok
!*** 125 ok
!*** 345 ok
!*** 235 ok
!*** 12345 ok
!---------------------------------
!> @brief converts a pyramid element into two tetrahedral elements.
!>
!> this subroutine takes the indices of the vertices of a pyramid element and
!> generates two tetrahedral elements from it. the indices of the vertices of
!> the pyramid are provided as inputs, and the indices of the vertices of the
!> resulting tetrahedral elements are returned as outputs.
!>
!> @param i1 index of the first vertex of the pyramid.
!> @param i2 index of the second vertex of the pyramid.
!> @param i3 index of the third vertex of the pyramid.
!> @param i4 index of the fourth vertex of the pyramid.
!> @param i5 index of the fifth vertex of the pyramid (apex).
!> @param tet1 array containing the indices of the vertices of the first tetrahedral element.
!> @param tet2 array containing the indices of the vertices of the second tetrahedral element.
   subroutine splitpyramidintotets(i1, i2, i3, i4, i5, tet1, tet2)
      integer, dimension(4) :: tet1, tet2
      logical :: swap13
! determine if i1-i3 swap is needed
      swap13 = max(i1, i3) > max(i2, i4)
      if (swap13) then
         tet1 = [i1, i2, i5, i3]
         tet2 = [i1, i5, i4, i3]
      else
         tet1 = [i1, i2, i5, i4]
         tet2 = [i2, i3, i5, i4]
      end if
   end subroutine splitpyramidintotets
!> \brief subroutine to perform operations on a tetrahedron.
!>
!> this subroutine takes four integer indices representing the vertices of a tetrahedron
!> and performs operations on the tetrahedron.
!>
!> \param[in] i1 integer index of the first vertex.
!> \param[in] i2 integer index of the second vertex.
!> \param[in] i3 integer index of the third vertex.
!> \param[in] i4 integer index of the fourth vertex.
!> \param[out] tet resulting tetrahedron after the operation.
   subroutine splitetintoitself(i1, i2, i3, i4, tet)
      integer, dimension(4)::tet
      tet = [i1, i2, i3, i4]
   end subroutine splitetintoitself
!* @brief     determines the symmetry of a tetrahedron defined by four nodes.
!*
!* @param[in]  n1  integer representing the first node of the tetrahedron.
!* @param[in]  n2  integer representing the second node of the tetrahedron.
!* @param[in]  n3  integer representing the third node of the tetrahedron.
!* @param[in]  n4  integer representing the fourth node of the tetrahedron.
!*
!* @return     integer value indicating the symmetry of the tetrahedron.
   integer function tetra4symmetry(n1, n2, n3, n4)
      integer::num
      integer, dimension(12)::list
      num = n1*1000 + n2*100 + n3*10 + n4*1
      list(1) = 1234
      list(2) = 2314
      list(3) = 3124
      list(4) = 1243
      list(5) = 2413
      list(6) = 4123
      list(7) = 1342
      list(8) = 3412
      list(9) = 4132
      list(10) = 2341
      list(11) = 3421
      list(12) = 4231
      do i = 1, 12
         if (num .eq. list(i)) then
            tetra4symmetry = 1
            return
         end if
      end do
      tetra4symmetry = -1
   end function tetra4symmetry
!*** number of edges
!*** chk0
   integer function nedges(ity)
      ne = 0
      select case (ity)
      case (1)
         ne = 0
      case (2)
         ne = 1
      case (3)
         ne = 3
      case (4)
         ne = 6
      case default
         stop "unavailable element topology"
      end select
      nedges = ne
   end function nedges
!*** list of edges according to topology
!*** chk0
   subroutine listedges(ity, ll)
      integer, dimension(2, 12)::ll
      ll = 0
      select case (ity)
      case (1)
         return
      case (2)
         ll(1, 1) = 1
         ll(2, 1) = 2
      case (3)
         ll(1, 1) = 1
         ll(2, 1) = 2
         ll(1, 2) = 2
         ll(2, 2) = 3
         ll(1, 3) = 3
         ll(2, 3) = 1
      case (4)
         ll(1, 1) = 1
         ll(2, 1) = 2
         ll(1, 2) = 2
         ll(2, 2) = 3
         ll(1, 3) = 1
         ll(2, 3) = 3
         ll(1, 4) = 3
         ll(2, 4) = 4
         ll(1, 5) = 1
         ll(2, 5) = 4
         ll(1, 6) = 2
         ll(2, 6) = 4
      case default
         stop "wrong call to listedges"
      end select
   end subroutine listedges
!*** other node given an edge
!*** chk0
   integer function othernodeofedge(mesh, ino, iar)
      type(originalmeshtype)::mesh
      othernodeofedge = 0
      if (mesh%arno(mesh%arni(iar)) .eq. ino) then
         othernodeofedge = mesh%arno(mesh%arni(iar) + 1)
      else if (mesh%arno(mesh%arni(iar) + 1) .eq. ino) then
         othernodeofedge = mesh%arno(mesh%arni(iar))
      end if
   end function othernodeofedge
!> @brief subroutine to define the dimensions of a new mesh.
!>
!> this subroutine is responsible for setting up the dimensions
!> of a new mesh structure. it initializes and configures the
!> necessary parameters for the new mesh.
!>
!> @param newmesh the new mesh structure to be initialized.
   subroutine dimensionsnewmesh(newmesh)
      integer::last
      type(renewedmeshtype)::newmesh
      newmesh%nel = newmesh%npoint + newmesh%nbar2 + newmesh%ntria3 + newmesh%ntetra4
      allocate (newmesh%elni(newmesh%nel + 1))
      last = 0
      newmesh%elni(1) = 1
      do ipoint = 1, newmesh%npoint
         newmesh%elni(ipoint + 1 + last) = newmesh%elni(ipoint + last) + 2
      end do
      last = newmesh%npoint
      do ibar2 = 1, newmesh%nbar2
         newmesh%elni(ibar2 + 1 + last) = newmesh%elni(ibar2 + last) + 3
      end do
      last = last + newmesh%nbar2
      do itria3 = 1, newmesh%ntria3
         newmesh%elni(itria3 + 1 + last) = newmesh%elni(itria3 + last) + 4
      end do
      last = last + newmesh%ntria3
      do itetra4 = 1, newmesh%ntetra4
         newmesh%elni(itetra4 + 1 + last) = newmesh%elni(itetra4 + last) + 5
      end do
      allocate (newmesh%elno(newmesh%elni(newmesh%nel + 1) - 1))
   end subroutine dimensionsnewmesh
!* @brief     determines the local edge of a given element in the mesh.
!*
!* @param[in] mesh  the mesh data structure containing all elements.
!* @param[in] iel   the index of the element in the mesh.
!* @param[in] ing   the index of the node within the element.
!*
!* @return    the local edge index corresponding to the given node.
   integer function localedgeinelement(mesh, iel, ing)
      type(originalmeshtype)::mesh
      inl = 0
      isucc = 0
      do ik = mesh%elir(iel), mesh%elir(iel + 1) - 1
         inl = inl + 1
         if (mesh%elar(ik) .eq. ing) then
            isucc = 1
            exit
         end if
      end do
      if (isucc .eq. 1) then
         localedgeinelement = inl
      else
         localedgeinelement = 0
      end if
   end function localedgeinelement
!* @brief     determines if a node should be split between two given nodes.
!*
!* @param[in]  i1         the index of the first node.
!* @param[in]  i2         the index of the second node.
!* @param[in]  midnodes   an array containing the midpoints of edges.
!* @param[in]  edgenodes  an array containing the nodes on the edges.
!*
!* @return     an integer indicating whether the node should be split.
   integer function splitnode(i1, i2, midnodes, edgenodes)
      integer, dimension(2, 6)::edgenodes
      integer, dimension(6)::midnodes
      ie = edgefromnodeslocal(i1, i2, edgenodes)
      splitnode = midnodes(abs(ie))
   end function splitnode
!* @brief     determines the local edge index from given node indices.
!*
!* @param[in]  i1         the index of the first node.
!* @param[in]  i2         the index of the second node.
!* @param[in]  edgenodes  an array containing the node indices for each edge.
!*
!* @return     the local edge index corresponding to the given node indices.
   integer function edgefromnodeslocal(i1, i2, edgenodes)
      integer, dimension(2, 6)::edgenodes
      iedge = 0
      do i = 1, 6
         if (i1 .eq. edgenodes(1, i) .and. i2 .eq. edgenodes(2, i)) then
            iedge = i
         else if (i1 .eq. edgenodes(2, i) .and. i2 .eq. edgenodes(1, i)) then
            iedge = -i
         end if
      end do
      edgefromnodeslocal = iedge
   end function edgefromnodeslocal
!> \brief this subroutine determines nodes from two given nodes.
!>
!> \param[in] edgenodes an array containing the edge nodes.
!> \param[in] post an array containing the post nodes.
!> \param[in] prev an array containing the previous nodes.
!> \param[in] i1 an integer representing the first node index.
!> \param[in] i2 an integer representing the second node index.
!> \param[in] i3 an integer representing the third node index.
!> \param[in] i4 an integer representing the fourth node index.
!>
!> this subroutine processes the given edge nodes and determines the
!> corresponding nodes based on the provided indices.
   subroutine nodesfrom2nodes(edgenodes, post, prev, i1, i2, i3, i4)
      integer, dimension(2, 6)::edgenodes
      integer, dimension(6)::post, prev
      i1 = 0
      i4 = 0
      iedge = edgefromnodeslocal(i2, i3, edgenodes)
      if (iedge .gt. 0) then
         i1 = prev(iedge)
         i4 = post(iedge)
      else
         i4 = prev(-iedge)
         i1 = post(-iedge)
      end if
   end subroutine nodesfrom2nodes
!*** an edge given nodes
!*** chk0
   integer function edgefromnodes(mesh, ino1, ino2)
      type(originalmeshtype)::mesh
      jar = 0
      do ik = mesh%noir(ino1), mesh%noir(ino1 + 1) - 1
         iar = mesh%noar(ik)
         if (othernodeofedge(mesh, ino1, iar) .eq. ino2) then
            jar = iar
            exit
         end if
      end do
      edgefromnodes = jar
   end function edgefromnodes
!*** modular arithmetic
   integer function iclock(n, i)
      itmp = mod(i, n)
      if (itmp .eq. 0) then
         iclock = n
      else
         iclock = itmp
      end if
      if (iclock .le. 0) iclock = iclock + ceiling(-1.0d00*iclock/n)*n
   end function iclock
!> this subroutine establishes the relationships between edges in
!> the given mesh. it processes the mesh data and marks the edges
!> accordingly.
!>
!> @param mesh the mesh data structure containing the elements and
!>             nodes of the mesh.
!> @param mark an array or data structure used to mark or label the
!>             edges based on certain criteria.
   subroutine edgerelations(mesh)
      type(originalmeshtype) :: mesh
      integer :: edgenodes(2, 12)
      integer, allocatable :: itemp(:), originalpairfromunique(:), uniquepairfromoriginal(:)
      integer, allocatable :: edgepairs(:, :)
      integer :: iel, ik, n1, n2, totaledges, uniqueedges, iar
      integer :: numelements, numedgeselem, numnodeselem
      numelements = mesh%nel
! calculate edge indices for each element
      if (.not. allocated(mesh%elir)) allocate (mesh%elir(numelements + 1))
      mesh%elir(1) = 1
      do iel = 1, numelements
         mesh%elir(iel + 1) = mesh%elir(iel) + nedges(mesh%elni(iel + 1) - mesh%elni(iel))
      end do
      totaledges = mesh%elir(numelements + 1) - 1
      if (allocated(mesh%elar)) DEALLOCATE (MESH%ELAR)
      allocate (mesh%elar(totaledges))
      if (allocated(edgepairs)) deallocate (edgepairs)
      allocate (edgepairs(2, totaledges))
      if (allocated(uniquepairfromoriginal)) deallocate (uniquepairfromoriginal)
      allocate (uniquepairfromoriginal(totaledges))
!*** trims it
      totaledges = 0
      do iel = 1, numelements
         numnodeselem = mesh%elni(iel + 1) - mesh%elni(iel)
         numedgeselem = nedges(numnodeselem)
         call listedges(numnodeselem, edgenodes)
         do ik = 1, numedgeselem
            n1 = mesh%elno(mesh%elni(iel) + edgenodes(1, ik) - 1)
            n2 = mesh%elno(mesh%elni(iel) + edgenodes(2, ik) - 1)
            if (n1 == n2) stop "same nodes"
            totaledges = totaledges + 1
            edgepairs(1:2, totaledges) = [n1, n2]
            mesh%elar(mesh%elir(iel) + ik - 1) = totaledges
         end do
      end do
      if (allocated(originalpairfromunique)) deallocate (originalpairfromunique)
      allocate (originalpairfromunique(totaledges))
      call getuniqueindicesforpairs(totaledges, edgepairs, originalpairfromunique, uniquepairfromoriginal, uniqueedges)
      if (allocated(mesh%arno)) deallocate (mesh%arno)
      allocate (mesh%arno(2*uniqueedges))
      if (allocated(mesh%arni)) deallocate (mesh%arni)
      allocate (mesh%arni(uniqueedges + 1))
      mesh%arni(1) = 1
      do iar = 1, uniqueedges
         mesh%arni(iar + 1) = mesh%arni(iar) + 2
      end do
      do iar = 1, uniqueedges
         mesh%arno(2*iar - 1:2*iar) = edgepairs(1:2, originalpairfromunique(iar))
      end do
      mesh%nar = uniqueedges
      do iel = 1, numelements
         do ik = mesh%elir(iel), mesh%elir(iel + 1) - 1
            mesh%elar(ik) = uniquepairfromoriginal(mesh%elar(ik))
         end do
      end do
! count edge connections per node and build row pointer array
      if (allocated(mesh%noir)) deallocate (mesh%noir)
      allocate (mesh%noir(mesh%nno + 1))
      mesh%noir = 0
      do iar = 1, mesh%nar
         n1 = mesh%arno(2*iar - 1)
         n2 = mesh%arno(2*iar)
         mesh%noir(n1) = mesh%noir(n1) + 1
         mesh%noir(n2) = mesh%noir(n2) + 1
      end do
      lol = mesh%noir(1)
      mesh%noir(1) = 1
      do in = 1, mesh%nno
         newv = mesh%noir(in) + lol
         in1 = in + 1
         lol = mesh%noir(in1)
         mesh%noir(in1) = newv
      end do
! build the edge index array
      if (.not. allocated(mesh%noar)) allocate (mesh%noar(mesh%noir(mesh%nno + 1) - 1))
      allocate (itemp(mesh%nno))
      itemp = 0
      do iar = 1, mesh%nar
         n1 = mesh%arno(2*iar - 1)
         n2 = mesh%arno(2*iar)
         mesh%noar(mesh%noir(n1) + itemp(n1)) = iar
         itemp(n1) = itemp(n1) + 1
         mesh%noar(mesh%noir(n2) + itemp(n2)) = iar
         itemp(n2) = itemp(n2) + 1
      end do
! deallocate temporary arrays
      deallocate (edgepairs, itemp)
   end subroutine edgerelations

!> this subroutine establishes the relationships between edges in
!> the given mesh. it processes the mesh data and marks the edges
!> accordingly.
!>
!> @param mesh the mesh data structure containing the elements and
!>             nodes of the mesh.
!> @param mark an array or data structure used to mark or label the
!>             edges based on certain criteria.
   subroutine edgerelationsmarked(mesh, mark)
      type(originalmeshtype) :: mesh
      integer ::  n1, n2, iar
      integer, allocatable ::  mark(:)
! set the marked edges
      allocate (mark(mesh%nar))
      mark = 0
      do i = 1, mesh%nmarkednodepairs
         n1 = mesh%markednodepairs(1, i)
         n2 = mesh%markednodepairs(2, i)
         iar = edgefromnodes(mesh, n1, n2)
         mark(iar) = abs(iar)
      end do
   end subroutine edgerelationsmarked
! Subroutine to perform iterative parallel merge sort on pairs
   subroutine sortpair2(arr, size)
      integer, intent(inout) :: arr(3, *)
      integer, intent(in) :: size
      integer :: curr_size, left_start, mid, right_end
      integer, allocatable :: temp(:, :)
! Start with subarrays of size 1, then 2, 4, 8, and so on
      curr_size = 1
      do while (curr_size <= size - 1)
!$omp parallel private(left_start, mid, right_end, temp)
!$omp do schedule(dynamic)
         do left_start = 1, size, 2*curr_size
! Calculate mid and right_end for current subarray
            mid = min(left_start + curr_size - 1, size)
            right_end = min(left_start + 2*curr_size - 1, size)
            if (mid < right_end) then
! Allocate temp array based on the size of the current subarray
               allocate (temp(3, right_end - left_start + 1))
! Call the merge subroutine
               call merge(arr, temp, left_start, mid, right_end)
! Deallocate temp after use
               deallocate (temp)
            end if
         end do
!$omp end do
!$omp end parallel
         curr_size = 2*curr_size
      end do
   end subroutine sortpair2
! Subroutine to merge two sorted halves
   subroutine merge(arr, temp, left, mid, right)
      integer, intent(inout) :: arr(3, *), temp(3, *)
      integer, intent(in) :: left, mid, right
      integer :: i, j, k
      i = left
      j = mid + 1
      k = 1
! Merge the two halves into temp array
      do while (i <= mid .and. j <= right)
         if (arr(1, i) < arr(1, j) .or. (arr(1, i) == arr(1, j) .and. arr(2, i) <= arr(2, j))) then
            temp(:, k) = arr(:, i)
            i = i + 1
         else
            temp(:, k) = arr(:, j)
            j = j + 1
         end if
         k = k + 1
      end do
! Copy the remaining elements from the left half (if any)
      do while (i <= mid)
         temp(:, k) = arr(:, i)
         i = i + 1
         k = k + 1
      end do
! Copy the remaining elements from the right half (if any)
      do while (j <= right)
         temp(:, k) = arr(:, j)
         j = j + 1
         k = k + 1
      end do
! Copy the sorted temp array back into the original array
      do i = 1, right - left + 1
         arr(:, left + i - 1) = temp(:, i)
      end do
   end subroutine merge
!*** get unique pair indices from
!*** a list of pairs
!> @brief this subroutine identifies unique pairs from a given list of pairs.
!>
!> @param numberofpairs the total number of pairs provided.
!> @param pairs an array containing the pairs to be processed.
!> @param originalpairfromunique an array that maps each unique pair to its original pair index.
!> @param uniquepairfromoriginal an array that maps each original pair to its unique pair index.
!> @param numberofuniquepairs the total number of unique pairs identified.
   subroutine getuniqueindicesforpairs(numberofpairs, pairs, originalpairfromunique, uniquepairfromoriginal, numberofuniquepairs)
      implicit none
      integer, intent(in) :: numberofpairs
      integer, intent(in) :: pairs(2, numberofpairs)
      integer, intent(out), allocatable :: originalpairfromunique(:)  ! originalpairfromunique(unique) = original index of the unique pair
      integer, intent(out) :: uniquepairfromoriginal(numberofpairs)   ! uniquepairfromoriginal(original) = unique index of the original index
      integer, intent(out) :: numberofuniquepairs
      integer :: i
      integer, allocatable :: ab(:, :)  ! ab(1, :) = normalized first elements
! ab(2, :) = normalized second elements
! ab(3, :) = original indices
      if (numberofpairs > 0) then
         allocate (ab(3, numberofpairs), originalpairfromunique(numberofpairs))   ! allocate arrays
! normalize pairs and store original indices
         do i = 1, numberofpairs
            if (pairs(1, i) <= pairs(2, i)) then  ! corrected variable name
               ab(1, i) = pairs(1, i)
               ab(2, i) = pairs(2, i)
            else
               ab(1, i) = pairs(2, i)
               ab(2, i) = pairs(1, i)
            end if
            ab(3, i) = i  ! store original index
         end do
! sort the ab array (including original indices)
         call sortpair2(ab, numberofpairs)
! initialize mappings
         numberofuniquepairs = 1
         originalpairfromunique(numberofuniquepairs) = ab(3, 1)  ! map unique index to original index
         uniquepairfromoriginal(ab(3, 1)) = numberofuniquepairs  ! map original index to unique index
         do i = 2, numberofpairs
            if ((ab(1, i) /= ab(1, i - 1) .or. ab(2, i) /= ab(2, i - 1))) then
! unique pair found
               numberofuniquepairs = numberofuniquepairs + 1
               originalpairfromunique(numberofuniquepairs) = ab(3, i)  ! map unique index to original index
               uniquepairfromoriginal(ab(3, i)) = numberofuniquepairs  ! map original index to unique index
            else
! duplicate found
               uniquepairfromoriginal(ab(3, i)) = uniquepairfromoriginal(ab(3, i - 1))
            end if
         end do
         deallocate (ab)
      else
         numberofuniquepairs = 0
      end if
   end subroutine getuniqueindicesforpairs
!> this subroutine creates edge nodes for a given mesh.
!>
!> @param mesh the original mesh data structure.
!> @param newmesh the new mesh data structure to be populated with edge nodes.
!> @param mark an array indicating which edges need new nodes.
   subroutine createedgenodes(mesh, newmesh, mark)
      type(originalmeshtype)::mesh
      type(renewedmeshtype)::newmesh
      integer, dimension(:)::mark
!*** total number of nodes, including existing ones
      nar = mesh%nar
      newmesh%nno = mesh%nno
      write (*, *) "nno before=", mesh%nno
      do iar = 1, nar
         if (mark(iar) .ne. 0) newmesh%nno = newmesh%nno + 1
      end do
      allocate (newmesh%parentnodes(2, newmesh%nno))
      do ino = 1, newmesh%nno
         newmesh%parentnodes(1:2, ino) = ino
      end do
      newmesh%nno = mesh%nno
      do iar = 1, nar
         if (mark(iar) .ne. 0) then
            newmesh%nno = newmesh%nno + 1
            mark(iar) = newmesh%nno
            ino1 = mesh%arno(mesh%arni(iar))
            ino2 = mesh%arno(mesh%arni(iar) + 1)
            if (ino1 .eq. ino2) stop "same node in edges???"
            newmesh%parentnodes(1:2, newmesh%nno) = [ino1, ino2]
         end if
      end do
      write (*, *) "nno after=", newmesh%nno
   end subroutine createedgenodes
end module remeshsimplex
