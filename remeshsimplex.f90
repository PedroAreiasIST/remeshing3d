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
MODULE remeshsimplex
!> @brief type definition for the original mesh structure.
!> this type is used to represent the original mesh in the simulation.
!> it contains all the necessary information about the mesh elements,
!> nodes, and connectivity required for the remeshing process.
  TYPE :: originalmeshtype
!*** nr. nodes:
    INTEGER::nno
!*** nr. edges
    INTEGER::nar
!*** conectivities
    INTEGER::nel
    INTEGER, DIMENSION(:), ALLOCATABLE::elno, elni
!*** edge data
    INTEGER, DIMENSION(:), ALLOCATABLE::elir, elar, arni, arno, noir, noar
!*** list of marked edges defined from nodes
    INTEGER::nmarkednodepairs
    INTEGER, DIMENSION(:, :), ALLOCATABLE::markednodepairs
  END TYPE originalmeshtype
!> @brief type definition for a renewed mesh structure.
!>
!> this type is used to represent a mesh that has been renewed or updated.
!> it contains all the necessary information and data structures to describe
!> the new state of the mesh after a remeshing operation.
  TYPE :: renewedmeshtype
    INTEGER::nno, npoint, nbar2, ntria3, ntetra4
    INTEGER::nel
    INTEGER, DIMENSION(:), ALLOCATABLE::elno, elni
    INTEGER, DIMENSION(:, :), ALLOCATABLE::parentnodes
  END TYPE renewedmeshtype
CONTAINS
!> @brief splits the given mesh into a new mesh.
!>
!> this subroutine takes an existing mesh and splits it into a new mesh.
!> the details of how the mesh is split are determined by the implementation
!> within this subroutine.
!>
!> @param[in] mesh the original mesh to be split.
!> @param[out] newmesh the resulting new mesh after splitting.
  SUBROUTINE splitmesh(mesh, newmesh)
    USE OMP_LIB
    REAL(8)::starttime, endtime
    TYPE(originalmeshtype)::mesh
    TYPE(renewedmeshtype)::newmesh
    INTEGER, DIMENSION(:), ALLOCATABLE::mark
!*** determine edge relation indices for the original mesh
    starttime = omp_get_wtime()
    CALL edgerelations(mesh)
    CALL edgerelationsmarked(mesh, mark)
!*** create edge nodes
    CALL createedgenodes(mesh, newmesh, mark)
!*** nummber of elements by type
    newmesh%npoint = 0
    newmesh%nbar2 = 0
    newmesh%ntria3 = 0
    newmesh%ntetra4 = 0
!*** dimensions the point elements
    CALL splitpoint0(mesh, newmesh%npoint)
!*** dimensions the split of 2-node elements (bar2)
    CALL splitbar20(mark, mesh, newmesh%nbar2)
!*** dimensions the split of 3-node elements (tria3)
    CALL splittria30(mark, mesh, newmesh%ntria3)
!*** dimensions the split of 4-node elements (tetra4)
    CALL splittetra40(mark, mesh, newmesh%ntetra4)
!*** now dimensions the essentials of the new mesh
    CALL dimensionsnewmesh(newmesh)
!*** point
    CALL splitpoint1(mesh, newmesh)
!*** bar2
    CALL splitbar21(mark, mesh, newmesh)
!*** tria3
    CALL splittria31(mark, mesh, newmesh)
!*** tetra4
    CALL splittetra41(mark, mesh, newmesh)
    endtime = omp_get_wtime()
    PRINT *, "splitting took", endtime-starttime, "seconds"
  END SUBROUTINE splitmesh
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
  SUBROUTINE splitpoint0(mesh, npoint)
    TYPE(originalmeshtype)::mesh
    npoint = 0
    DO iel = 1, mesh%nel
      IF (mesh%elni(iel+1)-mesh%elni(iel) .EQ. 1) THEN
        npoint = npoint+1
      END IF
    END DO
    WRITE (*, *) "npoint before=", npoint
    WRITE (*, *) "npoint after=", npoint
  END SUBROUTINE splitpoint0
!> @brief splits a point in the mesh based on the provided mark.
!>
!> this subroutine takes a mark and the current mesh, and generates a new mesh
!> by splitting a point according to the mark. it is used in mesh refinement
!> processes.
!>
!> @param[in] mark an integer or logical value indicating the point to be split.
!> @param[in] mesh the original mesh data structure.
!> @param[out] newmesh the new mesh data structure after the point has been split.
  SUBROUTINE splitpoint1(mesh, newmesh)
    TYPE(originalmeshtype)::mesh
    TYPE(renewedmeshtype)::newmesh
    ipoint = 0
    last = 0
    DO iel = 1, mesh%nel
      IF (mesh%elni(iel+1)-mesh%elni(iel) .EQ. 1) THEN
        ipoint = ipoint+1
        ing = mesh%elno(mesh%elni(iel))! node 1 of edge
        newmesh%elno(newmesh%elni(last+ipoint):newmesh%elni(last+1+ipoint)-1) = [ing, iel]
      END IF
    END DO
  END SUBROUTINE splitpoint1
!> @brief subroutine to split a bar element in a mesh.
!>
!> this subroutine is responsible for splitting a bar element in the mesh
!> based on the provided mark and mesh data. the number of bar elements
!> after the split is updated in `nbar2`.
!>
!> @param mark an integer array indicating the elements to be split.
!> @param mesh a derived type or array containing the mesh data.
!> @param nbar2 an integer representing the number of bar elements after the split.
  SUBROUTINE splitbar20(mark, mesh, nbar2)
    TYPE(originalmeshtype)::mesh
    INTEGER, DIMENSION(*)::mark
    nbar2 = 0
    DO iel = 1, mesh%nel
      IF (mesh%elni(iel+1)-mesh%elni(iel) .EQ. 2) THEN
        nbar2 = nbar2+1
      END IF
    END DO
    WRITE (*, *) "nbar2 before=", nbar2
    nbar2 = 0
    DO iel = 1, mesh%nel
      IF (mesh%elni(iel+1)-mesh%elni(iel) .EQ. 2) THEN
        ia = mesh%elar(mesh%elir(iel))
        IF (mark(ia) .NE. 0) THEN
          nbar2 = nbar2+2
        ELSE
          nbar2 = nbar2+1
        END IF
      END IF
    END DO
    WRITE (*, *) "nbar2 after=", nbar2
  END SUBROUTINE splitbar20
!> @brief splits a bar element into two new elements.
!>
!> this subroutine takes a marked element from the mesh and splits it into two new elements,
!> updating the mesh structure accordingly.
!>
!> @param mark an integer array indicating which elements are marked for splitting.
!> @param mesh the original mesh structure containing the elements to be split.
!> @param newmesh the updated mesh structure after the elements have been split.
  SUBROUTINE splitbar21(mark, mesh, newmesh)
    TYPE(originalmeshtype)::mesh
    TYPE(renewedmeshtype)::newmesh
    INTEGER, DIMENSION(*)::mark
    INTEGER::ibar2, ing1, ing2, ing3, ia, iel, last
    ibar2 = 0
    last = newmesh%npoint
    DO iel = 1, mesh%nel
      IF (mesh%elni(iel+1)-mesh%elni(iel) .EQ. 2) THEN
        ia = mesh%elar(mesh%elir(iel))
        ing2 = mark(ia)
        SELECT CASE (ing2)
        CASE (0)
          ing1 = mesh%elno(mesh%elni(iel))! node 1 of edge
          ing3 = mesh%elno(mesh%elni(iel)+1)! node 2 of edge
          ibar2 = ibar2+1
          newmesh%elno(newmesh%elni(last+ibar2):newmesh%elni(last+1+ibar2)-1) = [ing1, ing3, iel]
        CASE default
          ing1 = mesh%elno(mesh%elni(iel))! node 1 of edge
          ing3 = mesh%elno(mesh%elni(iel)+1)! node 2 of edge
          ibar2 = ibar2+1
          newmesh%elno(newmesh%elni(last+ibar2):newmesh%elni(last+1+ibar2)-1) = [ing1, ing2, iel]
          ibar2 = ibar2+1
          newmesh%elno(newmesh%elni(last+ibar2):newmesh%elni(last+1+ibar2)-1) = [ing2, ing3, iel]
        END SELECT
      END IF
    END DO
  END SUBROUTINE splitbar21
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
  SUBROUTINE splittria30(mark, mesh, ntria3)
    TYPE(originalmeshtype)::mesh
    INTEGER, DIMENSION(*) ::mark
    ntria3 = 0
    DO iel = 1, mesh%nel
      IF (mesh%elni(iel+1)-mesh%elni(iel) .EQ. 3) THEN
        ntria3 = ntria3+1
      END IF
    END DO
    WRITE (*, *) "ntria3 before=", ntria3
    ntria3 = 0
    DO iel = 1, mesh%nel
      IF (mesh%elni(iel+1)-mesh%elni(iel) .EQ. 3) THEN
        kk = 1
        DO ik = mesh%elir(iel), mesh%elir(iel+1)-1
          IF (mark(mesh%elar(ik)) .NE. 0) kk = kk+1
        END DO
        ntria3 = ntria3+kk
      END IF
    END DO
    WRITE (*, *) "ntria3 after=", ntria3
  END SUBROUTINE splittria30
!> @brief splits a triangular element into three new elements.
!>
!> this subroutine takes a marked triangular element and splits it into
!> new triangular elements. the new elements are added to the mesh.
!>
!> @param mark an integer indicating the marked triangular element to be split.
!> @param mesh the original mesh containing the triangular elements.
!> @param newmesh the new mesh where the split triangular elements will be added.
  SUBROUTINE splittria31(mark, mesh, newmesh)
    IMPLICIT NONE
    TYPE(originalmeshtype), INTENT(in) :: mesh
    TYPE(renewedmeshtype), INTENT(inout) :: newmesh
    INTEGER, DIMENSION(:), INTENT(in) :: mark
    INTEGER :: elem, markededgecount, edge, localedge, node1, node2, node3
    INTEGER :: newnode1, newnode2, newnode3, newnode4, newnode5, newnode6, tricounter
    INTEGER :: edgenode1, edgenode2, prevnode1, prevnode2, prevnode3, mark1, mark2, newmeshstart
    tricounter = 0  ! initialize triangle counter
    newmeshstart = newmesh%npoint+newmesh%nbar2
! loop through each element in the mesh
    DO elem = 1, mesh%nel
! check if the element has 3 nodes (triangular)
      IF (mesh%elni(elem+1)-mesh%elni(elem) == 3) THEN
        markededgecount = 0
! count the number of marked edges
        DO edge = mesh%elir(elem), mesh%elir(elem+1)-1
          IF (mark(mesh%elar(edge)) /= 0) markededgecount = markededgecount+1
        END DO
        SELECT CASE (markededgecount)
        CASE (0)
! no marked edges, retain original triangle
          newnode1 = mesh%elno(mesh%elni(elem))
          newnode2 = mesh%elno(mesh%elni(elem)+1)
          newnode3 = mesh%elno(mesh%elni(elem)+2)
          tricounter = tricounter+1
          newmesh%elno(newmesh%elni(tricounter+newmeshstart):newmesh%elni(tricounter+newmeshstart+1)-1) = [newnode1, newnode2, newnode3, elem]
        CASE (1)
! one marked edge, split triangle
          DO edge = mesh%elir(elem), mesh%elir(elem+1)-1
            IF (mark(mesh%elar(edge)) /= 0) THEN
              localedge = localedgeinelement(mesh, elem, mesh%elar(edge))
              node1 = localedge
              node2 = iclock(3, node1+1)
              node3 = iclock(3, node1+2)
              newnode1 = mesh%elno(mesh%elni(elem)-1+node1)
              newnode2 = mark(mesh%elar(edge))
              newnode3 = mesh%elno(mesh%elni(elem)-1+node2)
              newnode4 = mesh%elno(mesh%elni(elem)-1+node3)
! create two new triangles
              tricounter = tricounter+1
              newmesh%elno(newmesh%elni(tricounter+newmeshstart):newmesh%elni(tricounter+newmeshstart+1)-1) = [newnode1, newnode2, newnode4, elem]
              tricounter = tricounter+1
              newmesh%elno(newmesh%elni(tricounter+newmeshstart):newmesh%elni(tricounter+newmeshstart+1)-1) = [newnode2, newnode3, newnode4, elem]
            END IF
          END DO
        CASE (2)
          DO prevnode1 = 1, 3
            prevnode2 = iclock(3, prevnode1+1)
            prevnode3 = iclock(3, prevnode1+2)
            edgenode1 = mesh%elar(mesh%elir(elem)-1+prevnode1)
            edgenode2 = mesh%elar(mesh%elir(elem)-1+prevnode2)
            mark1 = mark(edgenode1)
            mark2 = mark(edgenode2)
            IF (mark1 /= 0 .AND. mark2 /= 0) THEN
              newnode1 = mesh%elno(mesh%elni(elem)-1+prevnode1)!node1
              newnode2 = mark(edgenode1)
              newnode3 = mesh%elno(mesh%elni(elem)-1+prevnode2)
              newnode4 = mark(edgenode2)
              newnode5 = mesh%elno(mesh%elni(elem)-1+prevnode3)
              IF (newnode2 > newnode4) THEN
                tricounter = tricounter+1
                newmesh%elno(newmesh%elni(tricounter+newmeshstart):newmesh%elni(tricounter+newmeshstart+1)-1) = [newnode1, newnode2, newnode5, elem]
                tricounter = tricounter+1
                newmesh%elno(newmesh%elni(tricounter+newmeshstart):newmesh%elni(tricounter+newmeshstart+1)-1) = [newnode2, newnode4, newnode5, elem]
              ELSE
                tricounter = tricounter+1
                newmesh%elno(newmesh%elni(tricounter+newmeshstart):newmesh%elni(tricounter+newmeshstart+1)-1) = [newnode1, newnode2, newnode4, elem]
                tricounter = tricounter+1
                newmesh%elno(newmesh%elni(tricounter+newmeshstart):newmesh%elni(tricounter+newmeshstart+1)-1) = [newnode1, newnode4, newnode5, elem]
              END IF
              tricounter = tricounter+1
              newmesh%elno(newmesh%elni(tricounter+newmeshstart):newmesh%elni(tricounter+newmeshstart+1)-1) = [newnode2, newnode3, newnode4, elem]
            END IF
          END DO
        CASE (3)
! three marked edges, fully subsplit triangle
          newnode1 = mesh%elno(mesh%elni(elem))
          newnode2 = mark(mesh%elar(mesh%elir(elem)))
          newnode3 = mesh%elno(mesh%elni(elem)+1)
          newnode4 = mark(mesh%elar(mesh%elir(elem)+1))
          newnode5 = mesh%elno(mesh%elni(elem)+2)
          newnode6 = mark(mesh%elar(mesh%elir(elem)+2))
! create four new triangles
          tricounter = tricounter+1
          newmesh%elno(newmesh%elni(tricounter+newmeshstart):newmesh%elni(tricounter+newmeshstart+1)-1) = [newnode1, newnode2, newnode6, elem]
          tricounter = tricounter+1
          newmesh%elno(newmesh%elni(tricounter+newmeshstart):newmesh%elni(tricounter+newmeshstart+1)-1) = [newnode2, newnode3, newnode4, elem]
          tricounter = tricounter+1
          newmesh%elno(newmesh%elni(tricounter+newmeshstart):newmesh%elni(tricounter+newmeshstart+1)-1) = [newnode4, newnode5, newnode6, elem]
          tricounter = tricounter+1
          newmesh%elno(newmesh%elni(tricounter+newmeshstart):newmesh%elni(tricounter+newmeshstart+1)-1) = [newnode2, newnode4, newnode6, elem]
        END SELECT
      END IF
    END DO
  END SUBROUTINE splittria31
!> @brief splits tetrahedral elements in a mesh.
!>
!> this subroutine is responsible for splitting tetrahedral elements in a mesh
!> based on certain criteria. it modifies the mesh structure by dividing the
!> tetrahedral elements.
!>
!> @param mark an array indicating which tetrahedral elements need to be split.
!> @param mesh the mesh data structure containing the tetrahedral elements.
!> @param ntetra4 the number of tetrahedral elements in the mesh.
  SUBROUTINE splittetra40(mark, mesh, ntetra4)
    TYPE(originalmeshtype)::mesh
    INTEGER, DIMENSION(*)::mark
    INTEGER, DIMENSION(6)::midnodes
!*** determines the new number of elements
    ntetra4 = 0
    DO iel = 1, mesh%nel
      IF (mesh%elni(iel+1)-mesh%elni(iel) .EQ. 4) THEN
        ntetra4 = ntetra4+1
      END IF
    END DO
    WRITE (*, *) "ntetra4 before=", ntetra4
    ntetra4 = 0
    DO iel = 1, mesh%nel
      IF (mesh%elni(iel+1)-mesh%elni(iel) .EQ. 4) THEN
        midnodes = 0
        DO ik = mesh%elir(iel), mesh%elir(iel+1)-1
          iar = mesh%elar(ik)
          midnodes(ik-mesh%elir(iel)+1) = mark(iar)
        END DO
        CALL splittetwork0(midnodes, ntets)
        ntetra4 = ntetra4+ntets
      END IF
    END DO
    WRITE (*, *) "ntetra4 after=", ntetra4
  END SUBROUTINE splittetra40
!> splits a tetrahedron into smaller tetrahedra based on the provided
!> marking and mesh information.
!>
!> this subroutine takes an existing mesh and a marking array to
!> generate a new refined mesh by splitting the tetrahedra.
!>
!> @param mark an array indicating which tetrahedra to split.
!> @param mesh the original mesh data structure.
!> @param newmesh the new mesh data structure after splitting.
  SUBROUTINE splittetra41(mark, mesh, newmesh)
    TYPE(originalmeshtype)::mesh
    TYPE(renewedmeshtype)::newmesh
    INTEGER, DIMENSION(:) :: mark
    INTEGER, DIMENSION(4, 8)::tets
    INTEGER, DIMENSION(6)::midnodes
    itetra4 = 0
    last = newmesh%npoint+newmesh%nbar2+newmesh%ntria3
    DO iel = 1, mesh%nel
      IF (mesh%elni(iel+1)-mesh%elni(iel) .EQ. 4) THEN
        midnodes = 0
        DO ik = mesh%elir(iel), mesh%elir(iel+1)-1
          iar = mesh%elar(ik)
          midnodes(ik-mesh%elir(iel)+1) = mark(iar)
        END DO
        CALL splittetwork1(mesh%elno(mesh%elni(iel):mesh%elni(iel+1)-1), midnodes, ntets, tets)
        DO i = 1, ntets
          itetra4 = itetra4+1
          newmesh%elno(newmesh%elni(itetra4+last):newmesh%elni(itetra4+last+1)-1) = [tets(1, i), tets(2, i), tets(3, i), tets(4, i), iel]
        END DO
      END IF
    END DO
  END SUBROUTINE splittetra41
!> @brief splits tetrahedral elements in the mesh.
!>
!> this subroutine is responsible for splitting the work tetrahedral elements
!> of the mesh. it takes the current nodes, mid-nodes, and the number of tetrahedral
!> elements as input parameters and performs the necessary operations to split them.
!>
!> @param nodes an array containing the coordinates of the nodes in the mesh.
!> @param midnodes an array containing the coordinates of the mid-nodes in the mesh.
!> @param ntets the number of tetrahedral elements in the mesh.
  SUBROUTINE splittetwork0(midnodes, ntets)
    INTEGER, DIMENSION(6)::midnodes, marked
    INTEGER, DIMENSION(3, 4)::nodeedges
    INTEGER, DIMENSION(5)::imarks
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
    DO i = 1, 6
      IF (midnodes(i) .GT. 0) THEN
        nmarked = nmarked+1
        marked(nmarked) = i
      END IF
    END DO
!*** number of marks per node
!*** and local edges
    imarks = 1
    DO i = 1, 4
      imarks(i+1) = imarks(i)
      DO j = 1, 3
        ie = nodeedges(j, i)! all elementwise edges of each node
        IF (midnodes(ie) .GT. 0) THEN! if this edge is marked
          imarks(i+1) = imarks(i+1)+1! updates the pointer
        END IF
      END DO
    END DO
!*** now the division
!*** according to the cases
    SELECT CASE (nmarked)
    CASE (0)
      ntets = 1
    CASE (1)
      ntets = 2
    CASE (2)
      i1 = 0
      i2 = 0
      mm = 0
      DO i = 1, 4
        mm = max(mm, imarks(i+1)-imarks(i))
      END DO
!*** one node with two marked edges:
      IF (mm .EQ. 2) THEN
        ntets = 3
      ELSE
        ntets = 4
      END IF
    CASE (3)
!*** now with 3 marked edges
!*** find case
      n4 = 0
      n3 = 0
      DO i = 1, 4
        IF (imarks(i+1)-imarks(i) .EQ. 3) THEN
          n4 = i
        ELSE IF (imarks(i+1)-imarks(i) .EQ. 0) THEN
          n3 = i
        END IF
      END DO
      IF (n3 .NE. 0) THEN
        ntets = 4
      ELSE IF (n4 .NE. 0) THEN
        ntets = 4
      ELSE
        ntets = 5
      END IF
    CASE (4)
      ntets = 6
    CASE (5)
      ntets = 7
    CASE (6)
      ntets = 8
    END SELECT
  END SUBROUTINE splittetwork0
!> splits the network of tetrahedral elements.
!>
!> this subroutine takes the existing nodes and midnodes of the tetrahedral
!> elements and splits them into smaller tetrahedra.
!>
!> @param nodes an array containing the coordinates of the original nodes.
!> @param midnodes an array containing the coordinates of the midpoints of the edges.
!> @param ntets the number of tetrahedral elements.
!> @param tets an array containing the indices of the nodes that form each tetrahedral element.
  SUBROUTINE splittetwork1(nodes, midnodes, ntets, tets)
! inputs
    INTEGER, DIMENSION(4) :: nodes          ! node indices of the tetrahedron's vertices
    INTEGER, DIMENSION(6) :: midnodes       ! midpoint node indices for each edge (0 if not present)
! outputs
    INTEGER :: ntets                        ! number of tetrahedra after division
    INTEGER, DIMENSION(4, *) :: tets        ! node indices of the resulting tetrahedra
! internal variables
    INTEGER, DIMENSION(6) :: marked         ! indices of marked edges (edges with midnodes)
    INTEGER, DIMENSION(6) :: post, prev     ! next and previous nodes for each edge
    INTEGER, DIMENSION(2, 6) :: edgenodes    ! node indices for each edge
    INTEGER, DIMENSION(3, 4) :: nodeedges    ! edge indices connected to each node
    INTEGER, DIMENSION(5) :: imarks         ! start indices in jmarks for each node
    INTEGER, DIMENSION(12) :: jmarks        ! marked edges per node
    INTEGER :: nmarked                      ! number of marked edges
    INTEGER :: i, j, ie, mm, ie1, itemp
    INTEGER :: n1, n2, n3, n4, n5, n6, n7, n8, n9, n10
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
    DO i = 1, 6
      IF (midnodes(i) .GT. 0) THEN
        nmarked = nmarked+1
        marked(nmarked) = i
      END IF
    END DO
! determine the number of marked edges per node and record the edges
    imarks(1) = 1
    jmarks = 0
    DO i = 1, 4
      imarks(i+1) = imarks(i)
      DO j = 1, 3
        ie = nodeedges(j, i)            ! edge connected to node i
        IF (midnodes(ie) .GT. 0) THEN   ! if edge is marked
          jmarks(imarks(i+1)) = ie    ! record edge index
          imarks(i+1) = imarks(i+1)+1 ! update pointer
        END IF
      END DO
    END DO
! split the tetrahedron based on the number of marked edges
    SELECT CASE (nmarked)
    CASE (0)
! case #1: no marked edges (original tetrahedron)
      ntets = 1
      tets(1:4, 1) = nodes(1:4)
    CASE (1)
! case #2: one marked edge
! identify the nodes connected by the marked edge
      n1 = edgenodes(1, marked(1))
      n2 = edgenodes(2, marked(1))
! find the other two nodes (n3 and n4)
! of opposing edge
      CALL nodesfrom2nodes(edgenodes, post, prev, n3, n1, n2, n4)
! convert local node indices to global node indices
      n1 = nodes(n1)
      n2 = nodes(n2)
      n3 = nodes(n3)
      n4 = nodes(n4)
      n5 = midnodes(marked(1)) ! midpoint node of the marked edge
! build two new tetrahedra
      ntets = 2
      CALL splitetintoitself(n3, n4, n5, n2, tets(1:4, 1))
      CALL splitetintoitself(n4, n3, n5, n1, tets(1:4, 2))
    CASE (2)
! cases #3 and #4: two marked edges
! find the maximum number of marked edges connected to any node
      mm = 0
      DO i = 1, 4
        mm = max(mm, imarks(i+1)-imarks(i))
      END DO
      IF (mm == 2) THEN
! case #3: one node with two marked edges
! identify the node with 0 marked edges (n3) and the node with 2 marked edges (n4)
        n3 = 0
        n4 = 0
        DO i = 1, 4
          IF (imarks(i+1)-imarks(i) == 0) THEN
            n3 = i
          ELSE IF (imarks(i+1)-imarks(i) == 2) THEN
            n4 = i
          END IF
        END DO
        IF (n3 == 0 .OR. n4 == 0) THEN
          STOP "error: either n3 or n4 is zero in case #3"
        END IF
! get the remaining nodes (n1 and n2)
        CALL nodesfrom2nodes(edgenodes, post, prev, n1, n3, n4, n2)
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
        CALL splitpyramidintotets(n1, n2, n5, n6, n3, tets(1:4, 1), tets(1:4, 2))
        CALL splitetintoitself(n6, n5, n3, n4, tets(1:4, 3))
      ELSE
! case #4: each node has one marked edge
! identify two nodes connected by a marked edge
        ie1 = marked(1)
        n1 = edgenodes(1, ie1)
        n2 = edgenodes(2, ie1)
! find the other two nodes (n3 and n4)
        CALL nodesfrom2nodes(edgenodes, post, prev, n3, n1, n2, n4)
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
        CALL splitetintoitself(n2, n3, n5, n6, tets(1:4, 1))
        CALL splitetintoitself(n2, n6, n5, n4, tets(1:4, 2))
        CALL splitetintoitself(n1, n5, n6, n4, tets(1:4, 3))
        CALL splitetintoitself(n1, n5, n3, n6, tets(1:4, 4))
      END IF
    CASE (3)
! cases #5, #6, and #7: three marked edges
! identify nodes based on the number of marked edges connected
      n4 = 0
      n3 = 0
      DO i = 1, 4
        IF (imarks(i+1)-imarks(i) == 3) THEN
          n4 = i
        ELSE IF (imarks(i+1)-imarks(i) == 0) THEN
          n3 = i
        END IF
      END DO
      IF (n4 .NE. 0) THEN
! case #6: one node with three marked edges
        n1 = iclock(4, n4+1)   ! next node in cyclic order
        CALL nodesfrom2nodes(edgenodes, post, prev, n2, n1, n4, n3)
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
        CALL splitetintoitself(n4, n7, n6, n5, tets(1:4, 1))
        CALL splitprismintotets(n2, n7, n6, n1, n3, n5, tets(1:4, 2), tets(1:4, 3), tets(1:4, 4))
      ELSE IF (n3 .NE. 0) THEN
! case #5: one node with zero marked edges
        n4 = iclock(4, n3+1)
        CALL nodesfrom2nodes(edgenodes, post, prev, n1, n3, n4, n2)
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
        CALL splitetintoitself(n1, n6, n3, n5, tets(1:4, 1))
        CALL splitetintoitself(n2, n3, n6, n7, tets(1:4, 2))
        CALL splitetintoitself(n6, n3, n5, n7, tets(1:4, 3))
        CALL splitetintoitself(n5, n7, n3, n4, tets(1:4, 4))
      ELSE
! case #7: other case
        n1 = 0
        n2 = 0
        m1 = 0
        m2 = 0
        DO i = 1, 4
          IF (imarks(i+1)-imarks(i) == 1) THEN
            ie = jmarks(imarks(i))
            IF (n1 == 0) THEN
              n1 = i
              IF (edgenodes(1, ie) .EQ. n1) THEN
                n4 = edgenodes(2, ie)
                m1 = +1
              ELSE
                n4 = edgenodes(1, ie)
                m1 = -1
              END IF
            ELSE
              n2 = i
              IF (edgenodes(1, ie) .EQ. n2) THEN
                n3 = edgenodes(2, ie)
                m2 = +1
              ELSE
                n3 = edgenodes(1, ie)
                m2 = -1
              END IF
            END IF
          END IF
        END DO
        isymm = tetra4symmetry(n1, n2, n3, n4)
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
        CALL splitpyramidintotets(n6, n7, n4, n2, n5, tets(1:4, 1), tets(1:4, 2))
        CALL splitpyramidintotets(n3, n1, n5, n7, n6, tets(1:4, 3), tets(1:4, 4))
        CALL splitetintoitself(n1, n2, n6, n5, tets(1:4, 5))
        IF (isymm .EQ. -1) THEN
          DO itet = 1, 5
            n3 = tets(3, itet)
            tets(3, itet) = tets(4, itet)
            tets(4, itet) = n3
          END DO
        END IF
      END IF
    CASE (4)
! cases #8 and #9: four marked edges
! find the maximum number of marked edges connected to any node
      mm = 0
      DO i = 1, 4
        mm = max(mm, imarks(i+1)-imarks(i))
      END DO
      IF (mm == 3) THEN
! case #8: one node with three marked edges
        n3 = 0
        n4 = 0
        DO i = 1, 4
          IF (imarks(i+1)-imarks(i) == 1) THEN
            n3 = i
          ELSE IF (imarks(i+1)-imarks(i) == 3) THEN
            n4 = i
          END IF
        END DO
        CALL nodesfrom2nodes(edgenodes, post, prev, n1, n3, n4, n2)
! midpoint nodes
        n5 = splitnode(n1, n2, midnodes, edgenodes)
        n6 = splitnode(n2, n4, midnodes, edgenodes)
        n7 = splitnode(n3, n4, midnodes, edgenodes)
        n8 = splitnode(n1, n4, midnodes, edgenodes)
        IF (n5 .EQ. 0 .OR. n6 .EQ. 0 .OR. n7 .EQ. 0 .OR. n8 .EQ. 0) THEN
          WRITE (*, *) "n5,n6,n7,n8", n5, n6, n7, n8
          READ (*, *)
        END IF
! convert local node indices to global node indices
        n1 = nodes(n1)
        n2 = nodes(n2)
        n3 = nodes(n3)
        n4 = nodes(n4)
! build six new tetrahedra
        ntets = 6
        CALL splitetintoitself(n8, n6, n7, n4, tets(1:4, 1))
        CALL splitetintoitself(n7, n8, n5, n6, tets(1:4, 2))
        CALL splitpyramidintotets(n2, n3, n7, n6, n5, tets(1:4, 3), tets(1:4, 4))
        CALL splitpyramidintotets(n1, n8, n7, n3, n5, tets(1:4, 5), tets(1:4, 6))
      ELSE
! case #9: edges form a quadrilateral
        ie = 0
        DO i = 1, 6
          IF (midnodes(i) == 0) THEN
            ie = i
            EXIT
          END IF
        END DO
        IF (ie == 0) THEN
          STOP "error: no unmarked edge found in case #9"
        END IF
        n1 = edgenodes(1, ie)
        n2 = edgenodes(2, ie)
        CALL nodesfrom2nodes(edgenodes, post, prev, n3, n1, n2, n4)
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
        CALL splitprismintotets(n8, n7, n3, n4, n5, n6, tets(1:4, 1), tets(1:4, 2), tets(1:4, 3))
        CALL splitprismintotets(n1, n2, n8, n5, n6, n7, tets(1:4, 4), tets(1:4, 5), tets(1:4, 6))
      END IF
    CASE (5)
! case #10: five marked edges
      ie = 0
      DO i = 1, 6
        IF (midnodes(i) == 0) THEN
          ie = i
          EXIT
        END IF
      END DO
      IF (ie == 0) THEN
        STOP "error: no unmarked edge found in case #10"
      END IF
      n1 = edgenodes(1, ie)
      n2 = edgenodes(2, ie)
      CALL nodesfrom2nodes(edgenodes, post, prev, n3, n1, n2, n4)
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
      CALL splitetintoitself(n9, n8, n5, n4, tets(1:4, 1))
      CALL splitetintoitself(n3, n6, n7, n8, tets(1:4, 2))
      CALL splitpyramidintotets(n6, n7, n9, n5, n8, tets(1:4, 3), tets(1:4, 4))
      CALL splitprismintotets(n1, n2, n9, n5, n6, n7, tets(1:4, 5), tets(1:4, 6), tets(1:4, 7))
    CASE (6)
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
      CALL splitetintoitself(n1, n9, n8, n5, tets(1:4, 1))
      CALL splitetintoitself(n2, n10, n9, n6, tets(1:4, 2))
      CALL splitetintoitself(n3, n8, n10, n7, tets(1:4, 3))
      CALL splitetintoitself(n5, n6, n7, n4, tets(1:4, 4))
      CALL splitoctointotets(n8, n6, n7, n5, n9, n10, tets(1:4, 5), tets(1:4, 6), tets(1:4, 7), tets(1:4, 8))
    END SELECT
  END SUBROUTINE splittetwork1
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
  SUBROUTINE splitprismintotets(i1, i2, i3, i4, i5, i6, tet1, tet2, tet3)
    INTEGER, DIMENSION(4) :: tet1, tet2, tet3
    LOGICAL :: swap46, swap13, swap16, swap25
! determine swaps
    swap46 = max(i4, i6) > max(i3, i5)
    swap13 = max(i1, i3) > max(i2, i4)
    swap16 = max(i1, i6) > max(i2, i5)
    swap25 = .NOT. swap16
    IF (swap46) THEN
      IF (swap13) THEN
        tet1 = [i1, i2, i6, i3]
        tet2 = [i1, i6, i4, i3]
        tet3 = [i1, i6, i5, i4]
      ELSE IF (swap16) THEN
        tet1 = [i1, i2, i6, i4]
        tet2 = [i1, i6, i5, i4]
        tet3 = [i2, i6, i4, i3]
      ELSE
        tet1 = [i1, i2, i5, i4]
        tet2 = [i2, i3, i6, i4]
        tet3 = [i2, i6, i5, i4]
      END IF
    ELSE
      IF (swap13) THEN
        IF (swap16) THEN
          tet1 = [i1, i2, i6, i3]
          tet2 = [i1, i3, i5, i4]
          tet3 = [i1, i6, i5, i3]
        ELSE
          tet1 = [i1, i3, i5, i4]
          tet2 = [i2, i5, i1, i3]
          tet3 = [i2, i6, i5, i3]
        END IF
      ELSE
        tet1 = [i1, i2, i5, i4]
        tet2 = [i2, i3, i5, i4]
        tet3 = [i2, i6, i5, i3]
      END IF
    END IF
  END SUBROUTINE splitprismintotets
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
  SUBROUTINE splitoctointotets(i1, i2, i3, i4, i5, i6, tet1, tet2, tet3, tet4)
    INTEGER, DIMENSION(4) :: tet1, tet2, tet3, tet4
    LOGICAL :: swap3546
! determine if diagonal swap between (35) and (46) is needed
    swap3546 = max(i3, i5) > max(i4, i6)
    IF (swap3546) THEN
      tet1 = [i1, i3, i4, i5]
      tet2 = [i1, i3, i5, i6]
      tet3 = [i3, i4, i5, i2]
      tet4 = [i3, i5, i6, i2]
    ELSE
      tet1 = [i1, i4, i6, i3]
      tet2 = [i1, i6, i4, i5]
      tet3 = [i2, i4, i6, i5]
      tet4 = [i2, i6, i4, i3]
    END IF
  END SUBROUTINE splitoctointotets
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
  SUBROUTINE splitpyramidintotets(i1, i2, i3, i4, i5, tet1, tet2)
    INTEGER, DIMENSION(4) :: tet1, tet2
    LOGICAL :: swap13
! determine if i1-i3 swap is needed
    swap13 = max(i1, i3) > max(i2, i4)
    IF (swap13) THEN
      tet1 = [i1, i2, i5, i3]
      tet2 = [i1, i5, i4, i3]
    ELSE
      tet1 = [i1, i2, i5, i4]
      tet2 = [i2, i3, i5, i4]
    END IF
  END SUBROUTINE splitpyramidintotets
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
  SUBROUTINE splitetintoitself(i1, i2, i3, i4, tet)
    INTEGER, DIMENSION(4)::tet
    tet = [i1, i2, i3, i4]
  END SUBROUTINE splitetintoitself
!* @brief     determines the symmetry of a tetrahedron defined by four nodes.
!*
!* @param[in]  n1  integer representing the first node of the tetrahedron.
!* @param[in]  n2  integer representing the second node of the tetrahedron.
!* @param[in]  n3  integer representing the third node of the tetrahedron.
!* @param[in]  n4  integer representing the fourth node of the tetrahedron.
!*
!* @return     integer value indicating the symmetry of the tetrahedron.
  INTEGER FUNCTION tetra4symmetry(n1, n2, n3, n4)
    INTEGER::num
    INTEGER, DIMENSION(12)::list
    num = n1*1000+n2*100+n3*10+n4*1
    list = [1234, 2143, 3412, 4321, 2314, 3124, 2431, 4132, 3241, 4213, 1342, 1423]
    DO i = 1, 12
      IF (num .EQ. list(i)) THEN
        tetra4symmetry = 1
        RETURN
      END IF
    END DO
    tetra4symmetry = -1
  END FUNCTION tetra4symmetry
!*** number of edges
!*** chk0
  INTEGER FUNCTION nedges(ity)
    ne = 0
    SELECT CASE (ity)
    CASE (1)
      ne = 0
    CASE (2)
      ne = 1
    CASE (3)
      ne = 3
    CASE (4)
      ne = 6
    CASE default
      STOP "unavailable element topology"
    END SELECT
    nedges = ne
  END FUNCTION nedges
!*** list of edges according to topology
!*** chk0
  SUBROUTINE listedges(ity, ll)
    INTEGER, DIMENSION(2, 12)::ll
    ll = 0
    SELECT CASE (ity)
    CASE (1)
      RETURN
    CASE (2)
      ll(1, 1) = 1
      ll(2, 1) = 2
    CASE (3)
      ll(1, 1) = 1
      ll(2, 1) = 2
      ll(1, 2) = 2
      ll(2, 2) = 3
      ll(1, 3) = 3
      ll(2, 3) = 1
    CASE (4)
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
    CASE default
      STOP "wrong call to listedges"
    END SELECT
  END SUBROUTINE listedges
!*** other node given an edge
!*** chk0
  INTEGER FUNCTION othernodeofedge(mesh, ino, iar)
    TYPE(originalmeshtype)::mesh
    othernodeofedge = 0
    IF (mesh%arno(mesh%arni(iar)) .EQ. ino) THEN
      othernodeofedge = mesh%arno(mesh%arni(iar)+1)
    ELSE IF (mesh%arno(mesh%arni(iar)+1) .EQ. ino) THEN
      othernodeofedge = mesh%arno(mesh%arni(iar))
    END IF
  END FUNCTION othernodeofedge
!> @brief subroutine to define the dimensions of a new mesh.
!>
!> this subroutine is responsible for setting up the dimensions
!> of a new mesh structure. it initializes and configures the
!> necessary parameters for the new mesh.
!>
!> @param newmesh the new mesh structure to be initialized.
  SUBROUTINE dimensionsnewmesh(newmesh)
    INTEGER::last
    TYPE(renewedmeshtype)::newmesh
    newmesh%nel = newmesh%npoint+newmesh%nbar2+newmesh%ntria3+newmesh%ntetra4
    ALLOCATE (newmesh%elni(newmesh%nel+1))
    last = 0
    newmesh%elni(1) = 1
    DO ipoint = 1, newmesh%npoint
      newmesh%elni(ipoint+1+last) = newmesh%elni(ipoint+last)+2
    END DO
    last = newmesh%npoint
    DO ibar2 = 1, newmesh%nbar2
      newmesh%elni(ibar2+1+last) = newmesh%elni(ibar2+last)+3
    END DO
    last = last+newmesh%nbar2
    DO itria3 = 1, newmesh%ntria3
      newmesh%elni(itria3+1+last) = newmesh%elni(itria3+last)+4
    END DO
    last = last+newmesh%ntria3
    DO itetra4 = 1, newmesh%ntetra4
      newmesh%elni(itetra4+1+last) = newmesh%elni(itetra4+last)+5
    END DO
    ALLOCATE (newmesh%elno(newmesh%elni(newmesh%nel+1)-1))
  END SUBROUTINE dimensionsnewmesh
!* @brief     determines the local edge of a given element in the mesh.
!*
!* @param[in] mesh  the mesh data structure containing all elements.
!* @param[in] iel   the index of the element in the mesh.
!* @param[in] ing   the index of the node within the element.
!*
!* @return    the local edge index corresponding to the given node.
  INTEGER FUNCTION localedgeinelement(mesh, iel, ing)
    TYPE(originalmeshtype)::mesh
    inl = 0
    isucc = 0
    DO ik = mesh%elir(iel), mesh%elir(iel+1)-1
      inl = inl+1
      IF (mesh%elar(ik) .EQ. ing) THEN
        isucc = 1
        EXIT
      END IF
    END DO
    IF (isucc .EQ. 1) THEN
      localedgeinelement = inl
    ELSE
      localedgeinelement = 0
    END IF
  END FUNCTION localedgeinelement
!* @brief     determines if a node should be split between two given nodes.
!*
!* @param[in]  i1         the index of the first node.
!* @param[in]  i2         the index of the second node.
!* @param[in]  midnodes   an array containing the midpoints of edges.
!* @param[in]  edgenodes  an array containing the nodes on the edges.
!*
!* @return     an integer indicating whether the node should be split.
  INTEGER FUNCTION splitnode(i1, i2, midnodes, edgenodes)
    INTEGER, DIMENSION(2, 6)::edgenodes
    INTEGER, DIMENSION(6)::midnodes
    ie = edgefromnodeslocal(i1, i2, edgenodes)
    splitnode = midnodes(abs(ie))
  END FUNCTION splitnode
!* @brief     determines the local edge index from given node indices.
!*
!* @param[in]  i1         the index of the first node.
!* @param[in]  i2         the index of the second node.
!* @param[in]  edgenodes  an array containing the node indices for each edge.
!*
!* @return     the local edge index corresponding to the given node indices.
  INTEGER FUNCTION edgefromnodeslocal(i1, i2, edgenodes)
    INTEGER, DIMENSION(2, 6)::edgenodes
    iedge = 0
    DO i = 1, 6
      IF (i1 .EQ. edgenodes(1, i) .AND. i2 .EQ. edgenodes(2, i)) THEN
        iedge = i
      ELSE IF (i1 .EQ. edgenodes(2, i) .AND. i2 .EQ. edgenodes(1, i)) THEN
        iedge = -i
      END IF
    END DO
    edgefromnodeslocal = iedge
  END FUNCTION edgefromnodeslocal
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
  SUBROUTINE nodesfrom2nodes(edgenodes, post, prev, i1, i2, i3, i4)
    INTEGER, DIMENSION(2, 6)::edgenodes
    INTEGER, DIMENSION(6)::post, prev
    i1 = 0
    i4 = 0
    iedge = edgefromnodeslocal(i2, i3, edgenodes)
    IF (iedge .GT. 0) THEN
      i1 = prev(iedge)
      i4 = post(iedge)
    ELSE
      i4 = prev(-iedge)
      i1 = post(-iedge)
    END IF
  END SUBROUTINE nodesfrom2nodes
!*** an edge given nodes
!*** chk0
  INTEGER FUNCTION edgefromnodes(mesh, ino1, ino2)
    TYPE(originalmeshtype)::mesh
    jar = 0
    DO ik = mesh%noir(ino1), mesh%noir(ino1+1)-1
      iar = mesh%noar(ik)
      IF (othernodeofedge(mesh, ino1, iar) .EQ. ino2) THEN
        jar = iar
        EXIT
      END IF
    END DO
    edgefromnodes = jar
  END FUNCTION edgefromnodes
!*** modular arithmetic
  INTEGER FUNCTION iclock(n, i)
    itmp = mod(i, n)
    IF (itmp .EQ. 0) THEN
      iclock = n
    ELSE
      iclock = itmp
    END IF
    IF (iclock .LE. 0) iclock = iclock+ceiling(-1.0d00*iclock/n)*n
  END FUNCTION iclock
!> this subroutine establishes the relationships between edges in
!> the given mesh. it processes the mesh data and marks the edges
!> accordingly.
!>
!> @param mesh the mesh data structure containing the elements and
!>             nodes of the mesh.
!> @param mark an array or data structure used to mark or label the
!>             edges based on certain criteria.
  SUBROUTINE edgerelations(mesh)
    TYPE(originalmeshtype) :: mesh
    INTEGER :: edgenodes(2, 12)
    INTEGER, ALLOCATABLE :: itemp(:), originalpairfromunique(:), uniquepairfromoriginal(:)
    INTEGER, ALLOCATABLE :: edgepairs(:, :)
    INTEGER :: iel, ik, n1, n2, totaledges, uniqueedges, iar
    INTEGER :: numelements, numedgeselem, numnodeselem
    numelements = mesh%nel
! calculate edge indices for each element
    IF (.NOT. allocated(mesh%elir)) ALLOCATE (mesh%elir(numelements+1))
    mesh%elir(1) = 1
    DO iel = 1, numelements
      mesh%elir(iel+1) = mesh%elir(iel)+nedges(mesh%elni(iel+1)-mesh%elni(iel))
    END DO
    totaledges = mesh%elir(numelements+1)-1
    IF (allocated(mesh%elar)) DEALLOCATE (MESH%ELAR)
    ALLOCATE (mesh%elar(totaledges))
    IF (allocated(edgepairs)) DEALLOCATE (edgepairs)
    ALLOCATE (edgepairs(2, totaledges))
    IF (allocated(uniquepairfromoriginal)) DEALLOCATE (uniquepairfromoriginal)
    ALLOCATE (uniquepairfromoriginal(totaledges))
!*** trims it
    totaledges = 0
    DO iel = 1, numelements
      numnodeselem = mesh%elni(iel+1)-mesh%elni(iel)
      numedgeselem = nedges(numnodeselem)
      CALL listedges(numnodeselem, edgenodes)
      DO ik = 1, numedgeselem
        n1 = mesh%elno(mesh%elni(iel)+edgenodes(1, ik)-1)
        n2 = mesh%elno(mesh%elni(iel)+edgenodes(2, ik)-1)
        IF (n1 == n2) STOP "same nodes"
        totaledges = totaledges+1
        edgepairs(1:2, totaledges) = [n1, n2]
        mesh%elar(mesh%elir(iel)+ik-1) = totaledges
      END DO
    END DO
    IF (allocated(originalpairfromunique)) DEALLOCATE (originalpairfromunique)
    ALLOCATE (originalpairfromunique(totaledges))
    CALL getuniqueindicesforpairs(totaledges, edgepairs, originalpairfromunique, uniquepairfromoriginal, uniqueedges)
    IF (allocated(mesh%arno)) DEALLOCATE (mesh%arno)
    ALLOCATE (mesh%arno(2*uniqueedges))
    IF (allocated(mesh%arni)) DEALLOCATE (mesh%arni)
    ALLOCATE (mesh%arni(uniqueedges+1))
    mesh%arni(1) = 1
    DO iar = 1, uniqueedges
      mesh%arni(iar+1) = mesh%arni(iar)+2
    END DO
    DO iar = 1, uniqueedges
      mesh%arno(2*iar-1:2*iar) = edgepairs(1:2, originalpairfromunique(iar))
    END DO
    mesh%nar = uniqueedges
    DO iel = 1, numelements
      DO ik = mesh%elir(iel), mesh%elir(iel+1)-1
        mesh%elar(ik) = uniquepairfromoriginal(mesh%elar(ik))
      END DO
    END DO
! count edge connections per node and build row pointer array
    IF (allocated(mesh%noir)) DEALLOCATE (mesh%noir)
    ALLOCATE (mesh%noir(mesh%nno+1))
    mesh%noir = 0
    DO iar = 1, mesh%nar
      n1 = mesh%arno(2*iar-1)
      n2 = mesh%arno(2*iar)
      mesh%noir(n1) = mesh%noir(n1)+1
      mesh%noir(n2) = mesh%noir(n2)+1
    END DO
    lol = mesh%noir(1)
    mesh%noir(1) = 1
    DO in = 1, mesh%nno
      newv = mesh%noir(in)+lol
      in1 = in+1
      lol = mesh%noir(in1)
      mesh%noir(in1) = newv
    END DO
! build the edge index array
    IF (.NOT. allocated(mesh%noar)) ALLOCATE (mesh%noar(mesh%noir(mesh%nno+1)-1))
    ALLOCATE (itemp(mesh%nno))
    itemp = 0
    DO iar = 1, mesh%nar
      n1 = mesh%arno(2*iar-1)
      n2 = mesh%arno(2*iar)
      mesh%noar(mesh%noir(n1)+itemp(n1)) = iar
      itemp(n1) = itemp(n1)+1
      mesh%noar(mesh%noir(n2)+itemp(n2)) = iar
      itemp(n2) = itemp(n2)+1
    END DO
! deallocate temporary arrays
    DEALLOCATE (edgepairs, itemp)
  END SUBROUTINE edgerelations
!> this subroutine establishes the relationships between edges in
!> the given mesh. it processes the mesh data and marks the edges
!> accordingly.
!>
!> @param mesh the mesh data structure containing the elements and
!>             nodes of the mesh.
!> @param mark an array or data structure used to mark or label the
!>             edges based on certain criteria.
  SUBROUTINE edgerelationsmarked(mesh, mark)
    TYPE(originalmeshtype) :: mesh
    INTEGER ::  n1, n2, iar
    INTEGER, ALLOCATABLE ::  mark(:)
! set the marked edges
    ALLOCATE (mark(mesh%nar))
    mark = 0
    DO i = 1, mesh%nmarkednodepairs
      n1 = mesh%markednodepairs(1, i)
      n2 = mesh%markednodepairs(2, i)
      iar = edgefromnodes(mesh, n1, n2)
      mark(iar) = abs(iar)
    END DO
  END SUBROUTINE edgerelationsmarked
! Subroutine to perform iterative parallel merge sort on pairs
  SUBROUTINE sortpair2(arr, size)
    INTEGER, INTENT(inout) :: arr(3, *)
    INTEGER, INTENT(in) :: size
    INTEGER :: curr_size, left_start, mid, right_end
    INTEGER, ALLOCATABLE :: temp(:, :)
! Start with subarrays of size 1, then 2, 4, 8, and so on
    curr_size = 1
    DO WHILE (curr_size <= size-1)
!$omp parallel private(left_start, mid, right_end, temp)
!$omp do schedule(dynamic)
      DO left_start = 1, size, 2*curr_size
! Calculate mid and right_end for current subarray
        mid = min(left_start+curr_size-1, size)
        right_end = min(left_start+2*curr_size-1, size)
        IF (mid < right_end) THEN
! Allocate temp array based on the size of the current subarray
          ALLOCATE (temp(3, right_end-left_start+1))
! Call the merge subroutine
          CALL merge(arr, temp, left_start, mid, right_end)
! Deallocate temp after use
          DEALLOCATE (temp)
        END IF
      END DO
!$omp end do
!$omp end parallel
      curr_size = 2*curr_size
    END DO
  END SUBROUTINE sortpair2
! Subroutine to merge two sorted halves
  SUBROUTINE merge(arr, temp, left, mid, right)
    INTEGER, INTENT(inout) :: arr(3, *), temp(3, *)
    INTEGER, INTENT(in) :: left, mid, right
    INTEGER :: i, j, k
    i = left
    j = mid+1
    k = 1
! Merge the two halves into temp array
    DO WHILE (i <= mid .AND. j <= right)
      IF (arr(1, i) < arr(1, j) .OR. (arr(1, i) == arr(1, j) .AND. arr(2, i) <= arr(2, j))) THEN
        temp(:, k) = arr(:, i)
        i = i+1
      ELSE
        temp(:, k) = arr(:, j)
        j = j+1
      END IF
      k = k+1
    END DO
! Copy the remaining elements from the left half (if any)
    DO WHILE (i <= mid)
      temp(:, k) = arr(:, i)
      i = i+1
      k = k+1
    END DO
! Copy the remaining elements from the right half (if any)
    DO WHILE (j <= right)
      temp(:, k) = arr(:, j)
      j = j+1
      k = k+1
    END DO
! Copy the sorted temp array back into the original array
    DO i = 1, right-left+1
      arr(:, left+i-1) = temp(:, i)
    END DO
  END SUBROUTINE merge
!*** get unique pair indices from
!*** a list of pairs
!> @brief this subroutine identifies unique pairs from a given list of pairs.
!>
!> @param numberofpairs the total number of pairs provided.
!> @param pairs an array containing the pairs to be processed.
!> @param originalpairfromunique an array that maps each unique pair to its original pair index.
!> @param uniquepairfromoriginal an array that maps each original pair to its unique pair index.
!> @param numberofuniquepairs the total number of unique pairs identified.
  SUBROUTINE getuniqueindicesforpairs(numberofpairs, pairs, originalpairfromunique, uniquepairfromoriginal, numberofuniquepairs)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: numberofpairs
    INTEGER, INTENT(in) :: pairs(2, numberofpairs)
    INTEGER, INTENT(out), ALLOCATABLE :: originalpairfromunique(:)  ! originalpairfromunique(unique) = original index of the unique pair
    INTEGER, INTENT(out) :: uniquepairfromoriginal(numberofpairs)   ! uniquepairfromoriginal(original) = unique index of the original index
    INTEGER, INTENT(out) :: numberofuniquepairs
    INTEGER :: i
    INTEGER, ALLOCATABLE :: ab(:, :)  ! ab(1, :) = normalized first elements
! ab(2, :) = normalized second elements
! ab(3, :) = original indices
    IF (numberofpairs > 0) THEN
      ALLOCATE (ab(3, numberofpairs), originalpairfromunique(numberofpairs))   ! allocate arrays
! normalize pairs and store original indices
      DO i = 1, numberofpairs
        IF (pairs(1, i) <= pairs(2, i)) THEN  ! corrected variable name
          ab(1, i) = pairs(1, i)
          ab(2, i) = pairs(2, i)
        ELSE
          ab(1, i) = pairs(2, i)
          ab(2, i) = pairs(1, i)
        END IF
        ab(3, i) = i  ! store original index
      END DO
! sort the ab array (including original indices)
      CALL sortpair2(ab, numberofpairs)
! initialize mappings
      numberofuniquepairs = 1
      originalpairfromunique(numberofuniquepairs) = ab(3, 1)  ! map unique index to original index
      uniquepairfromoriginal(ab(3, 1)) = numberofuniquepairs  ! map original index to unique index
      DO i = 2, numberofpairs
        IF ((ab(1, i) /= ab(1, i-1) .OR. ab(2, i) /= ab(2, i-1))) THEN
! unique pair found
          numberofuniquepairs = numberofuniquepairs+1
          originalpairfromunique(numberofuniquepairs) = ab(3, i)  ! map unique index to original index
          uniquepairfromoriginal(ab(3, i)) = numberofuniquepairs  ! map original index to unique index
        ELSE
! duplicate found
          uniquepairfromoriginal(ab(3, i)) = uniquepairfromoriginal(ab(3, i-1))
        END IF
      END DO
      DEALLOCATE (ab)
    ELSE
      numberofuniquepairs = 0
    END IF
  END SUBROUTINE getuniqueindicesforpairs
!> this subroutine creates edge nodes for a given mesh.
!>
!> @param mesh the original mesh data structure.
!> @param newmesh the new mesh data structure to be populated with edge nodes.
!> @param mark an array indicating which edges need new nodes.
  SUBROUTINE createedgenodes(mesh, newmesh, mark)
    TYPE(originalmeshtype)::mesh
    TYPE(renewedmeshtype)::newmesh
    INTEGER, DIMENSION(:)::mark
!*** total number of nodes, including existing ones
    nar = mesh%nar
    newmesh%nno = mesh%nno
    WRITE (*, *) "nno before=", mesh%nno
    DO iar = 1, nar
      IF (mark(iar) .NE. 0) newmesh%nno = newmesh%nno+1
    END DO
    ALLOCATE (newmesh%parentnodes(2, newmesh%nno))
    DO ino = 1, newmesh%nno
      newmesh%parentnodes(1:2, ino) = ino
    END DO
    newmesh%nno = mesh%nno
    DO iar = 1, nar
      IF (mark(iar) .NE. 0) THEN
        newmesh%nno = newmesh%nno+1
        mark(iar) = newmesh%nno
        ino1 = mesh%arno(mesh%arni(iar))
        ino2 = mesh%arno(mesh%arni(iar)+1)
        IF (ino1 .EQ. ino2) STOP "same node in edges???"
        newmesh%parentnodes(1:2, newmesh%nno) = [ino1, ino2]
      END IF
    END DO
    WRITE (*, *) "nno after=", newmesh%nno
  END SUBROUTINE createedgenodes
END MODULE remeshsimplex
