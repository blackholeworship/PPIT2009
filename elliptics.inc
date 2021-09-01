c-----------------------------------------------------------------------
c The following fortran subtroutine solves the elliptic equations
c in the problem using LAPACK's banded matrix solver DGBSV.
c-----------------------------------------------------------------------
        integer i,Nr,l_iter,j,c_iter
        real*8 alpha_res_nrm,psi_res_nrm,beta_res_nrm,res_nrm
        real*8 dr2,drsqrd,termp,term1,term4,term5,bterm,termm1
        real*8 r_iter,correction_weight
        logical trace,vtrace,solve_by_relax
        integer relax_iters,max_iter
        parameter (trace=.false.,vtrace=.false.,solve_by_relax=.false.)
        parameter (relax_iters=10000,correction_weight=0.5d0)
        logical first_warn,first_AH,second_warn,do_relax
        data first_warn/.true./,first_AH/.true./,second_warn/.true./
        data c_iter/0/
        save first_warn,second_warn
        save c_iter

        !---------------------------------------------------------------
        ! variables needed for DGBSV (see dgbsv.f in the LAPACK 
        ! distribution for the definitions)
        !
        ! For simplicity declare all variables global. For compatability
        ! will all versions of f77 we must used a constat Nr_max to
        ! define the arrays 
        !---------------------------------------------------------------

        integer Nr_max
        parameter (Nr_max=100000)
        integer N,NRHS,LDAB,INFO,KL,KU,LDB
        parameter (KL=2,KU=2,NRHS=1,LDAB=2*KL+KU+1)
        real*8 AB(2*KL+KU+1,Nr_max),B(Nr_max)
        integer IPIV(Nr_max)

        if (fixed_background.ne.0) then
          return
        end if

        Nr=g1_Nr
        LDB=Nr
        N=Nr
        dr2=2*dr
        drsqrd=dr*dr

        if (Nr.gt.Nr_max) then
           write(*,*) 'elliptics.inc: Error ... Nr_max too small'
           write(*,*) '  Nr_max=',Nr_max,' Nr=',Nr
           write(*,*) '  increase and recompile'
           stop
        end if

        !--------------------------------------------------------------
        ! Solve the elliptics via Newton iteration to within
        ! specified tolerance.  The matrix A(i,j)
        ! holds the Jacobian, i.e. the dependence of the linearized
        ! equation for the unknown at location r_i=(i-1)*dr
        ! w.r.t the variable at location r_j=(j-1)*dr. B(i)
        ! is minus the residual for that equation before the
        ! call to DGBSV, and holds the correction afterwards.
        ! The mapping to LAPACK storage is AB(KL+KU+1+i-j,j) = A(i,j)
        !--------------------------------------------------------------

        res_nrm=elliptics_tol+1
        l_iter=0

        if (trace) write(*,*) 'Solving Elliptics ...'

        max_iter=elliptics_max_iter
        do_relax=.false.
        if (.not.first_warn.and.solve_by_relax) then
           do_relax=.true.
           max_iter=relax_iters
        end if

        do while(res_nrm.gt.elliptics_tol.and.l_iter.lt.max_iter)
           l_iter=l_iter+1
           psi_res_nrm=0
           alpha_res_nrm=0
           beta_res_nrm=0
           if (trace) write(*,*) ' l_iter ',l_iter

           !-----------------------------------------------------------
           ! First, the Hamiltonian constraint for psi_np1 if 
           ! constrained
           !-----------------------------------------------------------
           if (free_psi.ne.0) then
              do i=1,Nr
                 psi_np1(i)=psi_evo_np1(i)
              end do
           else
              do i=1,Nr
                 B(i)=0
                 do j=i-2,i+2
                    AB(KL+KU+1+i-min(max(j,1),Nr),min(max(j,1),Nr))=0
                 end do
              end do
           
              !--------------------------------------------------------
              ! psi: r=0 B.C. psi'=0
              !--------------------------------------------------------
              i=1
              B(i)=-(-3*psi_np1(1)+4*psi_np1(2)-psi_np1(3))
              psi_res_nrm=psi_res_nrm+abs(B(i))
              j=1
              AB(KL+KU+1+i-j,j)=-3
              j=2
              AB(KL+KU+1+i-j,j)=4
              j=3
              AB(KL+KU+1+i-j,j)=-1

              if (do_relax) psi_np1(i)=psi_np1(i)+B(i)/AB(KL+KU+1+i-i,i)

              !--------------------------------------------------------
              ! psi: interior points
              !--------------------------------------------------------
              do i=2,Nr-1
                 bterm=(beta_np1(i+1)-beta_np1(i-1))/dr2
     &          - beta_np1(i)/r(i)
                 term5=(bterm/alpha_np1(i))**2/12
                 term1=cpi*(xi_np1(i)**2+PI_np1(i)**2)
                 termp=(2/r(i))/dr2

                 B(i)=-((psi_np1(i+1)-2*psi_np1(i)+psi_np1(i-1))/drsqrd 
     &               + (psi_np1(i+1)-psi_np1(i-1))*termp
     &               + psi_np1(i)**5*term5
     &               + psi_np1(i)*term1 )
                 psi_res_nrm=psi_res_nrm + abs(B(i))

                 j=i-1
                 AB(KL+KU+1+i-j,j)=1/drsqrd-termp
                 j=i
                 AB(KL+KU+1+i-j,j)=-2/drsqrd+5*psi_np1(i)**4*term5+term1
                 j=i+1
                 AB(KL+KU+1+i-j,j)=1/drsqrd+termp

                 if (do_relax) 
     &              psi_np1(i)=psi_np1(i)+B(i)/AB(KL+KU+1+i-i,i)
              end do

              !--------------------------------------------------------
              ! psi: outer boundary: (psi'+(\psi-1)/r)=0
              !--------------------------------------------------------
              i=Nr
              B(i)=-( (3*psi_np1(i)-4*psi_np1(i-1)+psi_np1(i-2))/dr2
     &            +(psi_np1(i)-1)/r(i) )
              psi_res_nrm=psi_res_nrm + abs(B(i))
              j=Nr
              AB(KL+KU+1+i-j,j)=3/dr2+1/r(i)
              j=Nr-1
              AB(KL+KU+1+i-j,j)=-4/dr2
              j=Nr-2
              AB(KL+KU+1+i-j,j)=1/dr2

              if (do_relax) psi_np1(i)=psi_np1(i)+B(i)/AB(KL+KU+1+i-i,i)
   
              if (.not.do_relax) then
                 call DGBSV(N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,INFO)
                 if (INFO.ne.0) then
                    write(*,*) 'Error: elliptics.inc solve for psi'
                    write(*,*) '       DGBSV returned info=',INFO
                    stop
                 end if
   
                 do i=1,Nr
                    psi_np1(i)=psi_np1(i)+B(i)*correction_weight
                 end do
              end if

              psi_res_nrm=psi_res_nrm/Nr

              if (trace) write(*,*) '    psi_res_nrm=',psi_res_nrm
              r_iter=l_iter
              if (vtrace) call vsxynt('psi',r_iter,r,psi_np1,nr)
           end if

           !-----------------------------------------------------------
           ! Now the momentum constraint for beta_np1
           !-----------------------------------------------------------
           do i=1,Nr
              B(i)=0
              do j=i-2,i+2
                 AB(KL+KU+1+i-min(max(j,1),Nr),min(max(j,1),Nr))=0
              end do
           end do
           
           !-----------------------------------------------------------
           ! beta: r=0 B.C. beta=0
           !-----------------------------------------------------------
           i=1
           B(i)=-(beta_np1(1))
           beta_res_nrm=beta_res_nrm+abs(B(i))
           j=1
           AB(KL+KU+1+i-j,j)=1

           if (do_relax) beta_np1(i)=beta_np1(i)+B(i)/AB(KL+KU+1+i-i,i)

           !-----------------------------------------------------------
           ! beta: interior points
           !-----------------------------------------------------------
           do i=2,Nr-1
              bterm=(6*(psi_np1(i+1)-psi_np1(i-1))/dr2/psi_np1(i)
     &                -(alpha_np1(i+1)-alpha_np1(i-1))/dr2/alpha_np1(i))
              term1=(-bterm/r(i)-2/r(i)**2)
              termp=(2/r(i)+bterm)/dr2

              B(i)=-((beta_np1(i+1)-2*beta_np1(i)+beta_np1(i-1))/drsqrd
     &               +(beta_np1(i+1)-beta_np1(i-1))*termp
     &               + beta_np1(i)*term1
     &               + 12*cpi*alpha_np1(i)*xi_np1(i)
     &                   *PI_np1(i)/psi_np1(i)**2 )
              beta_res_nrm=beta_res_nrm + abs(B(i))

              j=i-1
              AB(KL+KU+1+i-j,j)=1/drsqrd-termp
              j=i
              AB(KL+KU+1+i-j,j)=-2/drsqrd+term1
              j=i+1
              AB(KL+KU+1+i-j,j)=1/drsqrd+termp

              if (do_relax) 
     &           beta_np1(i)=beta_np1(i)+B(i)/AB(KL+KU+1+i-i,i)
           end do

           !-----------------------------------------------------------
           ! outer boundary: (beta'+\beta/r)=0
           !-----------------------------------------------------------
           i=Nr
           B(i)=-( (3*beta_np1(i)-4*beta_np1(i-1)+beta_np1(i-2))/dr2
     &            +(beta_np1(i))/r(i) )
           beta_res_nrm=beta_res_nrm + abs(B(i))
           j=Nr
           AB(KL+KU+1+i-j,j)=3/dr2+1/r(i)
           j=Nr-1
           AB(KL+KU+1+i-j,j)=-4/dr2
           j=Nr-2
           AB(KL+KU+1+i-j,j)=1/dr2

           if (do_relax) beta_np1(i)=beta_np1(i)+B(i)/AB(KL+KU+1+i-i,i)

           if (.not.do_relax) then
              call DGBSV(N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,INFO)

                 if (INFO.ne.0) then
                 write(*,*) 'Error: elliptics.inc solve for beta'
                 write(*,*) '       DGBSV returned info=',INFO
                 stop
              end if

              do i=1,Nr
                 beta_np1(i)=beta_np1(i)+B(i)*correction_weight
              end do
           end if

           beta_res_nrm=beta_res_nrm/Nr

           if (trace) write(*,*) '    beta_res_nrm=',beta_res_nrm
           r_iter=l_iter
           if (vtrace) call vsxynt('beta',r_iter,r,beta_np1,nr)

           !-----------------------------------------------------------
           ! Now the slicing condition for alpha_np1
           !-----------------------------------------------------------
           do i=1,Nr
              B(i)=0
              do j=i-2,i+2
                 AB(KL+KU+1+i-min(max(j,1),Nr),min(max(j,1),Nr))=0
              end do
           end do
           
           !-----------------------------------------------------------
           ! alpha: r=0 B.C. alpha'=0
           !-----------------------------------------------------------
           i=1
           B(i)=-(-3*alpha_np1(1)+4*alpha_np1(2)-alpha_np1(3))
           alpha_res_nrm=alpha_res_nrm+abs(B(i))
           j=1
           AB(KL+KU+1+i-j,j)=-3
           j=2
           AB(KL+KU+1+i-j,j)=4
           j=3
           AB(KL+KU+1+i-j,j)=-1

           if (do_relax) 
     &        alpha_np1(i)=alpha_np1(i)+B(i)/AB(KL+KU+1+i-i,i)

           !-----------------------------------------------------------
           ! alpha: interior points
           !-----------------------------------------------------------
           do i=2,Nr-1
              term1=-8*cpi*PI_np1(i)**2
              termp=(2/r(i)+2*(psi_np1(i+1)-psi_np1(i-1))
     &              /psi_np1(i))/dr2
              termm1=2*psi_np1(i)**4/3*
     &           (2*beta_np1(i)*(beta_np1(i+1)-beta_np1(i-1))/dr2/r(i)
     &            -beta_np1(i)**2/r(i)**2-
     &            ((beta_np1(i+1)-beta_np1(i-1))/dr2)**2)

              B(i)=-( (alpha_np1(i+1)-2*alpha_np1(i)
     &                               +alpha_np1(i-1))/drsqrd
     &               +(alpha_np1(i+1)-alpha_np1(i-1))*termp
     &               + termm1/alpha_np1(i)
     &               + alpha_np1(i)*term1 )
              alpha_res_nrm=alpha_res_nrm + abs(B(i))

              j=i-1
              AB(KL+KU+1+i-j,j)=1/drsqrd-termp
              j=i
              AB(KL+KU+1+i-j,j)=-2/drsqrd-termm1/alpha_np1(i)**2+term1
              j=i+1
              AB(KL+KU+1+i-j,j)=1/drsqrd+termp

              if (do_relax) 
     &           alpha_np1(i)=alpha_np1(i)+B(i)/AB(KL+KU+1+i-i,i)
           end do

           !-----------------------------------------------------------
           ! alpha: outer boundary: (alpha'+(\alpha-1)/r)=0
           !-----------------------------------------------------------
           i=Nr
           B(i)=-( (3*alpha_np1(i)-4*alpha_np1(i-1)+alpha_np1(i-2))/dr2
     &            +(alpha_np1(i)-1)/r(i) )
           alpha_res_nrm=alpha_res_nrm + abs(B(i))
           j=Nr
           AB(KL+KU+1+i-j,j)=3/dr2+1/r(i)
           j=Nr-1
           AB(KL+KU+1+i-j,j)=-4/dr2
           j=Nr-2
           AB(KL+KU+1+i-j,j)=1/dr2

           if (do_relax) 
     &        alpha_np1(i)=alpha_np1(i)+B(i)/AB(KL+KU+1+i-i,i)

           if (.not.do_relax) then
              call DGBSV(N,KL,KU,NRHS,AB,LDAB,IPIV,B,LDB,INFO)

              if (INFO.ne.0) then
                 write(*,*) 'Error: elliptics.inc solve for alpha'
                 write(*,*) '       DGBSV returned info=',INFO
                 stop
              end if

              do i=1,Nr
                 alpha_np1(i)=alpha_np1(i)+B(i)*correction_weight
              end do
           end if

           alpha_res_nrm=alpha_res_nrm/Nr

           if (trace) write(*,*) '    alpha_res_nrm=',alpha_res_nrm
           r_iter=l_iter
           if (vtrace) call vsxynt('alpha',r_iter,r,alpha_np1,nr)


           res_nrm=max(psi_res_nrm,beta_res_nrm,alpha_res_nrm)
        end do

        if (res_nrm.gt.elliptics_tol.and.first_warn) then
           write(*,*) '***********************************************'
           write(*,*) 'elliptics.inc: Warning ... failed to solve'
           write(*,*) '  elliptic equations to within ',elliptics_tol
           write(*,*) '  res_nrm=',res_nrm
           write(*,*) '    alpha_res_nrm=',alpha_res_nrm
           write(*,*) '    psi_res_nrm=',psi_res_nrm
           write(*,*) '    beta_res_nrm=',beta_res_nrm
           write(*,*) '(ONLY FIRST WARNING PRINTED)'
           first_warn=.false.
           if (solve_by_relax) then
              write(*,*) ' switching to relaxation from now onwards'
           end if
           write(*,*) '***********************************************'
        else if (res_nrm.gt.elliptics_tol.and.second_warn
     &           .and.solve_by_relax) then
           write(*,*) '***********************************************'
           write(*,*) 'elliptics.inc: Warning ... failed to solve'
           write(*,*) '  elliptic equations to within ',elliptics_tol
           write(*,*) '  using relaxation'
           write(*,*) '  res_nrm=',res_nrm
           write(*,*) '    alpha_res_nrm=',alpha_res_nrm
           write(*,*) '    psi_res_nrm=',psi_res_nrm
           write(*,*) '    beta_res_nrm=',beta_res_nrm
           write(*,*) '(ONLY FIRST WARNING PRINTED)'
           write(*,*) '***********************************************'
           second_warn=.false.
        end if

        call update2(psi_np1,beta_np1,alpha_np1,mass,g1_Nr,r,dr)
        if (vtrace) call vsxynt('mass',0.0d0,r,mass,nr)
