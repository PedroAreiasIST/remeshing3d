subroutine rLEVELSET(XX,YY,ZZ,val,tang,istyle)
  implicit real(8) (a-h,o-z)
  real(8)::val
  real(8),dimension(3)::tang
  x=xx+2.0d00
  y=yy+2.0d00
  z=zz
  val=(x**2.0d00/(1.0d00**2.0d00))+(y**2.0d00/(0.5d00**2.0d00))+(z**2.0d00/(0.06d00**2.0d00))-1.0d00
END subroutine rLEVELSET
!
!*** determinesxi
!
SUBROUTINE DETERMINESXI(XNMAX,XNMIN,XTRIAL,ISTYLE)
  IMPLICIT REAL(8) (A-H,O-Z)
  REAL(8),DIMENSION(3)::XNMIN,XNMAX,X1,X2,XTRIAL,tang
  XI=0.0D00
  X1=XNMAX
  X2=XNMIN
  call RLEVELSET(X1(1),X1(2),X1(3),S1,tang,ISTYLE)
  call RLEVELSET(X2(1),X2(2),X2(3),S2,tang,ISTYLE)
  IF(S1*S2.LT.0.0D00) THEN
     XI=S1/(S1-S2)
  END IF
  DO ID=1,3           
     XTRIAL(ID)=X1(ID)*(1.0D00-XI)+X2(ID)*XI
  END DO
  DO I=1,10
     call RLEVELSET(XTRIAL(1),XTRIAL(2),XTRIAL(3),S,tang,istyle)
     IF(S*S1.LT.0.0D00)THEN
        S2=S
        X2=XTRIAL
        XI=S1/(S1-S2)
     ELSEIF(S*S2.LT.0.0D00)THEN
        S1=S
        X1=XTRIAL
        XI=S1/(S1-S2)
     END IF
     DO ID=1,3           
        XTRIAL(ID)=X1(ID)*(1.0D00-XI)+X2(ID)*XI
     END DO
     xi=max(1.0d-3,min(1.0d00-1.0d-3,xi))
  END DO
END SUBROUTINE DETERMINESXI

!--------------------------------------
!*** put output in ensight format
!*** think about this with care
!*** number: number of steps
!*** filename:
!--------------------------------------
SUBROUTINE ensmaterial( &     
     number,filename, &
     nnoe,x, &
     nele,el_ni,el_no,escalar,evector)
  implicit real(8)(a-h,o-z)
  INTEGER,PARAMETER::lstring=80
  real(8),dimension(*)::escalar
  real(8),dimension(3,*)::evector
  real(8),dimension(3)::ev
  CHARACTER(lstring)::texto
  INTEGER,PARAMETER::mtp=4
  INTEGER,PARAMETER::mstr=79
  REAL(8),DIMENSION(3,*)::x
  CHARACTER(mstr)::ext
  CHARACTER(*)::filename
  CHARACTER(3)::pred
  INTEGER,DIMENSION(*)::el_ni,el_no
  CHARACTER(100)::cescalar="material"
  CHARACTER(100)::cevector="fiber"
  INTEGER,DIMENSION(:),ALLOCATABLE::el_list
  CHARACTER(3)::advstring
  character(6),dimension(4)::names=["point","bar2","tria3","tetra4"]
!****************
!*** check bounds
!****************
  IF(allocated(el_list)) DEALLOCATE(el_list)
  IF(nnoe.GT.99999999.OR.nele.GT.99999999) THEN
     WRITE(*,*) "Too large problem, more than 99999999 elements or nodes"
     STOP
  ENDIF
  IF(number.GT.9999) THEN
     WRITE(*,*) "Too many output files, please limit them to 9999 steps"
  ENDIF
!**********************
!*** geometry file .geo
!**********************
  IF(number.LE.9) THEN
     pred="000"
  ELSEIF(number.LE.99) THEN
     pred="00"
  ELSEIF(NUMBER.LE.999) THEN
     pred="0"
  ELSE
     pred=""
  ENDIF
  iun=14!iunit()
  WRITE(texto,"(i4)") number
  texto=trim(adjustl(pred))//trim(adjustl(texto))
  ext=trim(texto)
  OPEN(iun,file=trim(adjustl(filename))//".geo"//trim(adjustl(ext)),status="unknown")
  texto="Linha 1"
  WRITE(iun,"(a)") texto
  texto="Linha 2"
  WRITE(iun,"(a)") texto
  WRITE(iun,"(a)") "node id given"
  WRITE(iun,"(a)") "element id given"
  WRITE(iun,"(a)") "coordinates"
  WRITE(iun,"(i8)") nnoe
  DO in=1,nnoe
     WRITE(iun,"(i8,3e12.5)") in,(x(i,in),i=1,3)
  ENDDO
  WRITE(iun,"(a)") "part 1"
  texto="Only one part"
  WRITE(iun,"(a)") texto
  allocate(el_list(nele))
  DO itp=1,mtp
     ikount=0
     DO i=1,nele
!        if(escalar(i).le.1.1d00)cycle
        it=el_ni(i+1)-el_ni(i)-1
        IF(it.EQ.itp) THEN
           ikount=ikount+1
           el_list(ikount)=i
        ENDIF
     ENDDO
     IF(ikount.NE.0) THEN
        texto=names(itp)
        WRITE(iun,"(A)") adjustl(texto)
        WRITE(iun,"(I8)") ikount
        iti=0
        DO i=1,ikount
           WRITE(iun,"(*(I8))") el_list(i),(el_no(ikk),ikk=el_ni(el_list(i)),el_ni(el_list(i)+1)-2)
        ENDDO
     ENDIF
  ENDDO
  CLOSE(iun)
!*************
!*** case file
!*************
  iun=15!iunit()
  OPEN(iun,file=trim(adjustl(filename))//".case",status="unknown")
  WRITE(iun,"(a)") "FORMAT"
  WRITE(iun,"(a)") "type:   ensight"
  WRITE(iun,"(a)") "GEOMETRY"
  WRITE(iun,"(a)") "model: 1"//" "//trim(adjustl(filename))//".geo"//"****"
  WRITE(iun,"(a)") "VARIABLE"
  WRITE(iun,"(a)") "scalar per element:"//" 1 "//trim(cescalar)//" " &
       //trim(adjustl(filename))//trim(cescalar)//".res"//"****"
  WRITE(iun,"(a)") "vector per element:"//" 1 "//trim(cevector)//" " &
       //trim(adjustl(filename))//trim(cevector)//".res"//"****"
  WRITE(iun,"(a)") "TIME"
  WRITE(iun,"(a,i8)") "time set:"//" ",1
  WRITE(iun,"(a,i8)") "number of steps:"//" ",number+1
  WRITE(iun,"(a,i8)") "filename start number:"//" ",0
  WRITE(iun,"(a,i8)") "filename increment:"//" ",1
  WRITE(iun,"(a)") "time values:"
  nloop=(number+1)/5
  ntemp=mod(number+1,5)
  ik=0
  DO i=1,nloop
     iz=ik
     WRITE(iun,"(5i8)")(j,j=iz,iz+4)
     ik=iz+5
  ENDDO
  iz=ik
  WRITE(iun,"(5i8)")(i,i=iz,iz+ntemp-1)
  CLOSE(iun)
  iun=15
  OPEN(iun,file=trim(adjustl(filename))//trim(cescalar) &
       //".res"//trim(adjustl(ext)),status="unknown")
  WRITE(iun,"(a)") trim(cescalar)
  WRITE(iun,"(a)") "part 1"
  DO itp=1,mtp
     ikount=0
     DO j=1,nele
!        if(escalar(j).le.1.1d00)cycle
        it=el_ni(j+1)-el_ni(j)-1
        IF(it.EQ.itp) ikount=ikount+1
     ENDDO
     IF(ikount.NE.0) THEN
        WRITE(iun,"(a)") " "//trim(names(itp))//" "
        ik=1
        advstring="no"
        DO j=1,nele
!           if(escalar(j).le.1.1d00)cycle           
           es=escalar(j)
           it=el_ni(j+1)-el_ni(j)-1                   
           IF(it.EQ.itp) THEN
              ik=ik+1
              WRITE(iun,"(e12.5)",advance=trim(advstring))es
              IF(ik.EQ.6) THEN
                 ik=0
                 advstring="yes"
              ELSE
                 advstring="no"
              ENDIF
           ENDIF
        ENDDO
        write(iun,*)
     ENDIF
  ENDDO
  iun=16
  OPEN(iun,file=trim(adjustl(filename))//trim(cevector) &
       //".res"//trim(adjustl(ext)),status="unknown")
  WRITE(iun,"(a)") trim(cevector)
  WRITE(iun,"(a)") "part 1"
  DO itp=1,mtp
     ikount=0
     DO j=1,nele
!        if(escalar(j).le.1.1d00)cycle
        it=el_ni(j+1)-el_ni(j)-1
        IF(it.EQ.itp) ikount=ikount+1
     ENDDO
     IF(ikount.NE.0) THEN                 
        WRITE(iun,"(a)") " "//trim(names(itp))//" "        
        ik=1
        advstring="no"
        DO j=1,nele
!           if(escalar(j).le.1.1d00)cycle           
           ev=evector(1:3,j)
           it=el_ni(j+1)-el_ni(j)-1          
           IF(it.EQ.itp) THEN
              ik=ik+1
              WRITE(iun,"(3e12.5)",advance=trim(advstring))(evector(k,j),k=1,3)
              IF(ik.EQ.2) THEN
                 ik=0
                 advstring="yes"
              ELSE
                 advstring="no"
              ENDIF
           ENDIF
        ENDDO
        write(iun,*)
     ENDIF
  ENDDO
  close(iun)
ENDSUBROUTINE ensmaterial

module gidmeshmodule
contains
  subroutine readmeshessentials(umal,mesh,NOCO)
    use remeshsimplex
    IMPLICIT REAL(8)(A-H,O-Z)
    integer::umal
    type(originalmeshtype)::mesh
    LOGICAL::ICHECK
    INTEGER,PARAMETER::NL=13
    CHARACTER(140)::LEITU
    REAL(8),DIMENSION(:,:),ALLOCATABLE::NOCO
    INTEGER,DIMENSION(:),ALLOCATABLE::ITEMP
    INTEGER,DIMENSION(:),ALLOCATABLE::ELTYP
    INTEGER,DIMENSION(2)::LISTP1,LISTP2
    INTEGER,DIMENSION(2,12)::LL,LL2
!--------------------------------
!*** READING OF NODE COORDINATES
!--------------------------------
    DO I=1,NL
       READ(UMAL,*)
    ENDDO
    READ(UMAL,*) mesh%NNO
    write(*,*)"mesh%nno",mesh%nno
    ALLOCATE(NOCO(3,MESH%NNO))
    DO IN=1,mesh%NNO
       READ(UMAL,*) INO,(NOCO(IKK,INO),IKK=1,3)
       IF(INO.NE.IN) STOP "WRONG NODE NUMBERING IN THE MESH FILE"
    ENDDO
!--------------------------------------
!*** READING OF ELEMENT CONNECTIVITIES
!--------------------------------------
    READ(UMAL,*)
    READ(UMAL,*) mesh%NEL
    write(*,*)"nel=",mesh%nel
    allocate(mesh%elni(mesh%nel+1))
    mesh%elni(1)=1
    DO IEL=1,mesh%NEL
       READ(UMAL,*) JEL,ITY
       IF(JEL.NE.IEL) STOP "WRONG ELEMENT ORDERING IN THE MESH FILE"
       select case(ity)
       case(1)
          mesh%elni(iel+1)=mesh%elni(iel)+1
       case(2)
          mesh%elni(iel+1)=mesh%elni(iel)+2
       case(3)
          mesh%elni(iel+1)=mesh%elni(iel)+3
       case(5)
          mesh%elni(iel+1)=mesh%elni(iel)+4
       case default
          stop "unnaceptable element"
       end select
    ENDDO
    WRITE(*,*)"ELNI DONE",MESH%ELNI(MESH%NEL+1)
    ALLOCATE(MESH%ELNO(MESH%ELNI(MESH%NEL+1)-1))
    DO I=mesh%NEL,1,-1
       BACKSPACE(UMAL)
    ENDDO
    DO IEL=1,mesh%NEL     
       READ(UMAL,*) JEL,ITY,(MESH%ELNO(K),K=MESH%ELNI(JEL),MESH%ELNI(JEL+1)-1)
    ENDDO
  end SUBROUTINE READMESHESSENTIALS
end module gidmeshmodule
program asneiras
  do istyle=1,1!0!5
     call ellipsecut(istyle)
  end do
end program asneiras
subroutine ELLIPSECUT(istyle)
  use basfun
  use gidmeshmodule
  use remeshsimplex
  IMPLICIT REAL(8)(A-H,O-Z)
  CHARACTER(100)::lnm
  integer,dimension(:),allocatable::upper
  integer::umal
  LOGICAL::fich1
  INTEGER(2)::result2
  integer,dimension(4)::mk
  real(8),dimension(4)::sg
  real(8),dimension(3)::x1,x2,x3,x4,xmed,tang
  INTEGER::verif
  character(100)::texto,nomef,nomefbase
  real(8),dimension(:),allocatable::escalar
  real(8),dimension(:,:),allocatable::evector
  real(8),DIMENSION(:,:),ALLOCATABLE::NOCO,NOCO2
  type(originalmeshtype)::mesh,updatedmesh
  type(originalmeshtype)::originalmesh
  type(renewedmesh)::newmesh
  CALL getarg(1,texto(1:len(texto)))
  IF(texto.EQ."") THEN
     result2=system("ls *.malha > listfiles")
     iu=iunit()
     OPEN(iu,file="listfiles",status="unknown")
     nf=0
     DO
        READ(iu,"(a)",iostat=io) lnm
        IF(io.NE.0) EXIT
        nf=nf+1
     ENDDO
     IF(nf.EQ.1) THEN
        texto=lnm
     ELSE
        WRITE(*,"(a)",advance="no") "please insert the name of input file:"
        DO
           READ(*,"(a)") texto
           IF(texto.NE."") EXIT
        ENDDO
     ENDIF
  ENDIF
  nomef(1:)=texto(1:)
  nomef=adjustl(nomef)
  texto=" "
  CALL minusc(nomef)
  ii=index(nomef,"."//trim("malha"))
  IF(ii.NE.0) nomef(ii:)=" "
  jj=len_trim(nomef)
  IF(nomef(jj:jj).EQ.".") THEN
     nomef(jj:)=" "
  ENDIF
  nomefbase=nomef
  INQUIRE(file=trim(trim(nomefbase)//"."//"malha"),exist=fich1)
  IF(.NOT.fich1) THEN
     CALL purgadig(nomefbase)
     CALL simple(nomefbase)
  ENDIF
  WRITE(*,"(a)") "we are now trying to open",trim(adjustl(nomef))//"."//trim(adjustl("malha"))
  INQUIRE(file=trim(trim(nomefbase)//"."//"malha"),exist=fich1)
  IF(.NOT.fich1) THEN
     WRITE(*,"(a)") "the file ""malha"" appears to be inexistent"
     WRITE(*,"(a)") " this will soon cause a fatal error"
  ENDIF
  INQUIRE(file=trim(trim(nomef)//"."//"malha"),exist=fich1)
  IF(.NOT.fich1) THEN
     WRITE(*,"(a)") "the file ""malha"" appears to be inexistent"
     STOP
  ENDIF
  umal=12
  CALL openfl(nomefbase,"malha",umal,"desconhecido",.FALSE.,ierr)
  call readmeshessentials(umal,mesh,noco)
  originalmesh%nel=mesh%nel
  originalmesh%nno=mesh%nno
  originalmesh%elni=mesh%elni
  originalmesh%elno=mesh%elno
  originalmesh%nmarkednodepairs=0
  call edgerelations(originalmesh)
  allocate(originalmesh%markednodepairs(2,originalmesh%nar))
  do iar=1,originalmesh%nar
     in1=originalmesh%arno(originalmesh%arni(iar))
     in2=originalmesh%arno(originalmesh%arni(iar)+1)
     if(in2.gt.in1)then
        itemp=in1 
        in1=in2 
        in2=itemp 
     end if
     CALL RLEVELSET(NOCO(1,IN1),NOCO(2,IN1),NOCO(3,IN1),S1,TANG,ISTYLE)
     CALL RLEVELSET(NOCO(1,IN2),NOCO(2,IN2),NOCO(3,IN2),S2,TANG,ISTYLE)
     IF(S1*S2.LT.0.0D00) THEN
        ORIGINALMESH.NMARKEDNODEPAIRS=ORIGINALMESH.NMARKEDNODEPAIRS+1
        ORIGINALMESH.MARKEDNODEPAIRS(1,ORIGINALMESH.NMARKEDNODEPAIRS)=IN1
        ORIGINALMESH.MARKEDNODEPAIRS(2,ORIGINALMESH.NMARKEDNODEPAIRS)=IN2
     END IF
  end do
  write(*,*)"Now splits"
  call splitmesh(originalmesh,newmesh)
  Close(umal)
  updatedmesh%nel=newmesh%nel
  updatedmesh%nno=newmesh%nno
  updatedmesh%elni=newmesh%elni
  updatedmesh%elno=newmesh%elno
  allocate(noco2(3,updatedmesh%nno))
  write(*,*)"one"
  do ino=1,updatedmesh%nno
     ino1=newmesh%parentnodes(1,ino)
     ino2=newmesh%parentnodes(2,ino)
     if(ino1.ne.ino2)then
        in1=max(ino1,ino2)
        in2=min(ino1,ino2)
        xi=0.0d00
        CALL DETERMINESXI(NOCO(1:3,IN1),NOCO(1:3,IN2),noco2(1:3,ino),ISTYLE)
     else
        do id=1,3
           noco2(id,ino)=noco(id,ino1)
        end do
     end if
  end do
  write(*,*)"two"  
  umal=19  
  CALL openfl(nomefbase,"malha",umal,"desconhecido",.FALSE.,ierr)  
  WRITE(umal,"(A)") "This is a REMESHED malha file"
  DO I=1,11
     WRITE(umal,*)
  ENDDO
  WRITE(umal,"(A)") "NOW THE NODES:"
  WRITE(umal,"(I8)") updatedmesh%nno
  DO INO=1,updatedmesh%nno
     WRITE(umal,"(I8,3E18.6)") INO,noco2(1:3,INO)
  ENDDO
  WRITE(umal,"(A)") "NOW THE ELEMENTS:"
  WRITE(umal,"(I8)") updatedmesh%nel
  allocate(escalar(updatedmesh%nel))
  allocate(evector(3,updatedmesh%nel))  
  write(*,*)"two---2"  
  DO IEL=1,updatedmesh%nel
     ngash=updatedmesh%elni(iel+1)-updatedmesh%elni(iel)-1
     select case(updatedmesh%elni(iel+1)-updatedmesh%elni(iel)-1)
     case(1)
        ityp=1
     case(2)
        ityp=2
     case(3)
        ist=updatedmesh%elni(iel)
        n1=updatedmesh%elno(ist)
        n2=updatedmesh%elno(ist+1)
        n3=updatedmesh%elno(ist+2)
        x1=noco2(1:3,n1)
        x2=noco2(1:3,n2)
        x3=noco2(1:3,n3)        
        ityp=3
        xmed=0.3333333333d00*(x1+x2+x3)
        call rlevelset(xmed(1),xmed(2),xmed(3),s,tang,istyle)
        IF(s.LT.0.0D00)THEN
           texto="inside-s"
           escalar(iel)=1.0d00*(1.0d00+istyle)
        else
           texto="outside-s"
           escalar(iel)=1.0d00           
        end if
     case(4)
        ityp=5
        ist=updatedmesh%elni(iel)
!----------------------------------- 
        n1=updatedmesh%elno(ist)
        n2=updatedmesh%elno(ist+1)
        n3=updatedmesh%elno(ist+2)
        n4=updatedmesh%elno(ist+3)
        x1=noco2(1:3,n1)
        x2=noco2(1:3,n2)
        x3=noco2(1:3,n3)
        x4=noco2(1:3,n4)
        volt=voltetra(x1,x2,x3,x4)
        call rLEVELSET(x1(1),x1(2),x1(3),s1,tang,istyle)
        call rLEVELSET(x2(1),x2(2),x2(3),s2,tang,istyle)
        call rLEVELSET(x3(1),x3(2),x3(3),s3,tang,istyle)
        call rLEVELSET(x4(1),x4(2),x4(3),s4,tang,istyle)
        sg=[s1,s2,s3,s4]
        mk=0
        do in=1,4
           ino=updatedmesh%elno(ist-1+in)
           ino1=newmesh%parentnodes(1,ino)
           ino2=newmesh%parentnodes(2,ino)
           if(ino1.ne.ino2)then
              mk(in)=1         
           end if
        end do
        smin=1.0d30
        smax=-1.0d30
        if(mk(1).ne.1)smin=min(smin,s1)
        if(mk(2).ne.1)smin=min(smin,s2)
        if(mk(3).ne.1)smin=min(smin,s3)
        if(mk(4).ne.1)smin=min(smin,s4)
        if(mk(1).ne.1)smax=max(smax,s1)
        if(mk(2).ne.1)smax=max(smax,s2)
        if(mk(3).ne.1)smax=max(smax,s3)
        if(mk(4).ne.1)smax=max(smax,s4)
        if((mk(1).ne.0).or.(mk(2).ne.0).or.(mk(3).ne.0).or.(mk(4).ne.0))then
           np=0
           ng=0
           do k=1,4
              if(mk(k).eq.0)then
                 if(sg(k).lt.-1.0d-20)ng=ng+1
                 if(sg(k).gt.1.0d-20)np=np+1
              end if
           end do
           if(np*ng.ne.0)then
              write(*,*)"mk=",mk
              write(*,*)"s=",sg
              write(*,*)"iel=",iel
           end if
        end if
        if(-smin.gt.smax)then
           inside=1
        else
           inside=0
        end if
        if(inside.eq.1)then
           texto="inside-v"
           escalar(iel)=2.0d00
        else
           texto="outside-v"           
           escalar(iel)=1.0d00
        end if
!        write(*,*)"4,iel=",iel          
     end select
     WRITE(umal,"(20I8)",advance="no") IEL,ityp,(updatedmesh%elno(k),K=updatedmesh%elni(iel),updatedmesh%elni(iel+1)-2)
     write(umal,"(A)")"  "//trim(adjustl(texto))
!     write(*,*)"iel=",iel  
  ENDDO
  write(*,*)"three"
  write(umal,*)
  write(umal,*)0
!*** evector
  evector=0.0d000
  do iel=1,updatedmesh%nel
     if(updatedmesh%elni(iel+1)-updatedmesh%elni(iel)-1.eq.4)then
        xmed=0.0d00
        do ik=1,4
           xmed(1:3)=xmed(1:3)+noco2(1:3,updatedmesh%elno(updatedmesh%elni(iel)-1+ik))
        end do
        xmed=0.25d00*xmed
!        dz=dzfunction(xmed(2))
        call rLEVELSET(xmed(1),xmed(2),xmed(3),val,evector(1:3,iel),istyle)
!        evector(1,iel)=0.0d00
!        evector(2,iel)=1.0d00
!        evector(3,iel)=dz
     end if
  end do
  write(texto,*)istyle
  CALL ENSMATERIAL(0,TRIM(ADJUSTL(NOMEFBASE))//trim(adjustl(texto)), &
     UPDATEDMESH%NNO,NOCO2, &
     UPDATEDMESH%NEL,UPDATEDMESH%ELNI,UPDATEDMESH%ELNO,ESCALAR,EVECTOR)
  WRITE(*,*)"FOUR"
END subroutine ELLIPSECUT
