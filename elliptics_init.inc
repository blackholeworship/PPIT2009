c----------------------------------------------------------------------
c simply call update routine (see updates.f), but with future (np1)
c and past (n) time levels reversed.
c----------------------------------------------------------------------
        integer i,c_free_psi

        do i=1,g1_Nr
           alpha_n(i)=1
           psi_n(i)=1
           beta_n(i)=0
        end do

        ! always solve constraints here
   
        c_free_psi=0

        call evolve_e(PI_n,PI_np1,xi_n,xi_np1,psi_n,psi_np1,beta_n
     &  ,beta_np1,alpha_n,alpha_np1,psi_evo_n,psi_evo_np1,mass,
     &  g1_Nr,r,t,dr,cpi,elliptics_max_iter,elliptics_tol,
     &  fixed_background,c_free_psi)

        do i=1,g1_Nr
           psi_evo_n(i)=psi_n(i)
        end do
