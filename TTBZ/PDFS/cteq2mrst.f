       SUBROUTINE CTEQ6(X,SCALE,MODE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
       implicit none
       double precision X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU
       double precision Q,xsave,qsave,Ctq6Pdf,D,U
       integer mode,iset

         Q=SCALE
         Iset=MODE
         xsave=X
         qsave=Q
         U =         X * Ctq6Pdf(1,X,Q)
         D =         X * Ctq6Pdf(2,X,Q)
         USEA =      X * Ctq6Pdf(-1,X,Q)
         DSEA =      X * Ctq6Pdf(-2,X,Q)
         STR =       X * Ctq6Pdf(3,X,Q)
         CHM =       X * Ctq6Pdf(4,X,Q)
         BOT =       X * Ctq6Pdf(5,X,Q)
         GLU  =      X * Ctq6Pdf(0,X,Q)
         UPV=U-USEA
         DNV=D-DSEA
         X=xsave
         Q=qsave
        return
       end



       SUBROUTINE CTEQ10(X,SCALE,MODE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU)
       implicit none
       double precision X,SCALE,UPV,DNV,USEA,DSEA,STR,CHM,BOT,GLU
       double precision Q,xsave,qsave,Ct10Pdf,D,U
       integer mode,iset

         Q=SCALE
         Iset=MODE
         xsave=X
         qsave=Q
         U =         X * Ct10Pdf(1,X,Q)
         D =         X * Ct10Pdf(2,X,Q)
         USEA =      X * Ct10Pdf(-1,X,Q)
         DSEA =      X * Ct10Pdf(-2,X,Q)
         STR =       X * Ct10Pdf(3,X,Q)
         CHM =       X * Ct10Pdf(4,X,Q)
         BOT =       X * Ct10Pdf(5,X,Q)
         GLU  =      X * Ct10Pdf(0,X,Q)
         UPV=U-USEA
         DNV=D-DSEA
         X=xsave
         Q=qsave
        return
       end

