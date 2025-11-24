C     EPSZPD
      subroutine epszpd (omr, omi, pars, u, v, flag)
C     PARAMETER LIST
C     omr     = REAL      PART OF OM
C     omi     = IMAGINARY PART OF OM
C     U      = REAL      PART OF eps(om)
C     V      = IMAGINARY PART OF eps(om)
C     FLAG   = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
C     OCCUR OR NOT; TYPE LOGICAL;
C     THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
C     MEANING :
C     FLAG=.FALSE. : NO ERROR CONDITION
C     FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
C     BECOMES INACTIVE
C     OMR, OMI     ARE THE INPUT-PARAMETERS
C     U, V, FLAG  ARE THE OUTPUT-PARAMETERS
C     
C     
      implicit none
      DOUBLE PRECISION omr,omi,ky,kpar,omdi,etai,u,v,
     *     epsabs,epsrel,alim,blim,abserr,
     *     work(40000),resr,resi,zbb,bbi,minomdlim,tau,
     *     pars(5),sqrttwo,omsi,limsingsm,limsinglg
      double complex zaa,i,w,om
      INTEGER n,m,np1,nu,j,l,iwork(10000),limit
      LOGICAL A, B, FLAG
      integer nlimit,mf,nf,last,neval,ier,npts2
      PARAMETER (MINOMDLIM = -1e-6,
     *     nlimit=10000,epsrel=1.0e-4,epsabs=1.0e-8,
     *     sqrttwo =1.414213562373095,npts2=3,
     *     limsingsm=1.0e-12,limsinglg=1.0e-5)
      double precision spoints(3)
      double precision fepspd_re,fepspd_im,resfepspd_im,resfepspd_re
      EXTERNAL DQAG, fepspd_re,fepspd_im,resfepspd_im,resfepspd_re
      common /epscom/ omdi,omsi,etai,tau,ky,kpar,zbb,bbi,zaa,w,om
      FLAG = .FALSE.
      i=dcmplx(0,1)
      om=dcmplx(omr,omi)
      omdi=pars(1)
      etai=pars(2)
      tau=pars(3)
      ky=pars(4)
      kpar=pars(5)
      if(omdi.LT.MINOMDLIM) then
         zaa=-dcmplx(omr,omi)/omdi
         zbb=-kpar*sqrttwo/omdi
         w=zbb**2/4-zaa
         bbi=ky**2;
         omsi=-ky;
         if(dabs(dimag(w)).LT.limsingsm) then
            zaa=dble(zaa)+i*epsrel
            w=zbb**2/4-zaa
         endif
         Alim=0.0
         CALL DQAGI(Fepspd_re,alim,1,epsabs,epsrel,resr,abserr,neval,
     *        ier,nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         CALL DQAGI(Fepspd_im,alim,1,epsabs,epsrel,resi,abserr,neval,
     *        ier,nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         u=resr+1.0+1.0/tau;
         v=resi;
         if(omi.LT.0.AND.dble(w).GT.0) then
            Alim=-1.0+limsingsm
            Blim=1.0-limsingsm
            CALL DQAG(resFepspd_re,alim,blim,epsabs,epsrel,6,resr,
     *           abserr,neval,ier, nlimit,40000,last,iwork,work)
            if(ier.ne.0) goto 100
            CALL DQAG(resFepspd_im,alim,blim,epsabs,epsrel,6,resi,
     *           abserr,neval,ier,nlimit,40000,last,iwork,work)
            if(ier.ne.0) goto 100
            u=u-resr
            v=v-resi
         else if(omi.EQ.0.AND.dble(w).GT.0) then
            Alim=-1.0+limsingsm
            Blim=1.0-limsingsm
            CALL DQAG(resFepspd_re,alim,blim,epsabs,epsrel,6,resr,
     *           abserr,neval,ier, nlimit,40000,last,iwork,work)
            if(ier.ne.0) goto 100
            CALL DQAG(resFepspd_im,alim,blim,epsabs,epsrel,6,resi,
     *           abserr,neval,ier,nlimit,40000,last,iwork,work)
            if(ier.ne.0) goto 100
            u=u-0.5*resr
            v=v-0.5*resi
         endif
      endif
      RETURN
*
 100  FLAG = .TRUE.
      RETURN

c      write (*,*) "om:",om
c      write (*,*) "w:",w
c      write (*,*) "zbb:",zbb
c      write (*,*) "bbi:",bbi
*
      END

      double precision function Fepspd_re(s)
      double complex Fepspd
      external Fepspd
      double precision s
      Fepspd_re=dble(Fepspd(s))
      return
      end function Fepspd_re

      double precision function Fepspd_im(s)
      double complex Fepspd
      external Fepspd
      double precision s
      Fepspd_im=dimag(Fepspd(s))
      return
      end function Fepspd_im

      double complex function Fepspd(s)
      double precision s,limsinglg,xbr,Jr0,zbb,limsingsm
      double precision bbi,omdi,omsi,etai,tau,ky,kpar
      double complex z1,z2,zaa,G0,G2,weidGm,w,om,zr
      common /epscom/ omdi,omsi,etai,tau,ky,kpar,zbb,bbi,zaa,w,om
      external weidGm
      parameter (limsinglg = 1.0e-5,limsingsm=1.0e-12)
      zr=cdsqrt(4.0*w-2.0*s)
      z1=0.5*(zbb+zr)
      z2=0.5*(zbb-zr)
      G0=weidGm(z1,z2,0)
      G2=weidGm(z1,z2,2)
      xbr=dsqrt(bbi*2.0*s)
      jr0=dbesj0(xbr)
      Fepspd=Jr0**2*dexp(-s)*(
     *     (om-omsi*(1.0+(s-1.5)*etai))*G0-omsi*etai*G2
     *     )/omdi
      return
      end function Fepspd

      double precision function resFepspd_re(s)
      double complex resFepspd
      external resFepspd
      double precision s
      resFepspd_re=dble(resFepspd(s))
      return
      end function resFepspd_re

      double precision function resFepspd_im(s)
      double complex resFepspd
      external resFepspd
      double precision s
      resFepspd_im=dimag(resFepspd(s))
      return
      end function resFepspd_im

      double complex function resFepspd(mu)
      double precision mu,limsinglg,xbr,xbi,Jr0,Ji0,zbb,fnu
      double precision bbi,sqrtpi,omdi,etai,tau,ky,kpar,omsi
      integer ierr,nz
      EXTERNAL CBJ0
      double complex zaa,i,w,xb,J0,om
      common /epscom/ omdi,omsi,etai,tau,ky,kpar,zbb,bbi,zaa,w,om
      parameter (limsinglg = 1.0e-5, sqrtpi = 1.77245385090552)
      xb=2.0*cdsqrt(bbi*(1.0d0-mu**2)*w)
      i=dcmplx(0,1)
      call cbj0(xb,J0)
      resFepspd=i*4.0*J0**2*sqrtpi/omdi*cdsqrt(w)*
     *     cdexp(-(mu*cdsqrt(w)+0.5*zbb)**2-2.0*(1.0-mu**2)*w)*
     *     (om-omsi*(1.0+etai*(2.0*w*(1-mu**2)
     *     +(mu*cdsqrt(w)+0.5*zbb)**2-1.5)))
      return
      end function resFepspd

C     epsweid
      subroutine epsweid(om,pars,numza,res)
Cf2py double complex dimension(numza) :: om
Cf2py double precision dimension(5) :: pars
Cf2py double complex dimension(numza) :: res
      double precision pars(5),xi,yi,u,v
      double complex om(numza),res(numza)
      integer k,numza
      logical flag
      external epszpd,prerr
      double precision largestxi,largestyi,largeval
      parameter (largestxi=1E5,largestyi=1E5,largeval=1d120)
      do 10 k=1,numza
         xi=dble(om(k))
         yi=dimag(om(k))
         if (dabs(xi).gt.largestxi.or.dabs(yi).gt.largestyi) goto 120
         call epszpd(xi,yi,pars,u,v,flag)
         if(flag) goto 110
         res(k)=dcmplx(u,v)
 10   continue
      return
 110  return
 120  write(*,*) "om value out of bounds (+-1e4)! returning:", largeval
      res(k)=largeval
      continue
      end

C     epsweid
      subroutine sigweid(za,zb,b,anm,numn,numm,numza,res)
Cf2py double complex dimension(numza) :: za
Cf2py double precision :: zb
Cf2py double precision :: b
Cf2py double complex dimension(numn,numm) :: anm
Cf2py double complex dimension(numza) :: res
      double precision xi,yi,u,v,zb,b
      double complex za(numza),res(numza),anm(numn,numm),i
      integer k,numza,numn,numm,j,l
      logical flag
      external inmzpd,prerr
      double precision largestxi,largestyi,largeval,limtn
      parameter (largestxi=1E5,largestyi=1E5,largeval=1d120,
     *     limtn=1e-20)
      i=dcmplx(0,1);
      do 30 l=1,numza
         xi=dble(za(l))
         yi=dimag(za(l))
         if (dabs(xi).gt.largestxi.or.dabs(yi).gt.largestyi) goto 140
         call sigmazpd(xi,yi,zb(l),b,anm,numn,numm,u,v,flag)
         if(flag) goto 130
         res(l)=u+i*v;
 30   continue
      return
 130  call prerr(za(l),zb(l),b)
      write (*,*) "za:",za
      return
 140  write(*,*) "za value out of bounds (+-1e4)! returning:", largeval
      res(l)=largeval
      continue
      end

      subroutine sigmazpd (xi, yi, zb, bi, anm, 
     *     numn, numm, u, v, flag)
      DOUBLE PRECISION xi,yi,bi,zb,z1,z2,u,v, 
     *     epsabs,epsrel,alim,blim,abserr,spoints(3),
     *     work(40000),resr,resi,zbb,bbi,limsingsm,limtn
      double complex zaa,za,i,w,anml(40),anm(numn,numm)
      INTEGER np1,nu,j,l,iwork(10000),nml(40,2),ier,neval,numn,numm
      LOGICAL A, B, FLAG
      integer nlimit,mf,nf,last,npts2,k,numnml
      PARAMETER (limsingsm=1.0e-12,npts2=3,limtn=1.0e-20,
     *     epsrel=1.0e-4,epsabs=1.0e-8)
      double precision fsigpd_re,fsigpd_im,resfsigpd_im,resfsigpd_re
      EXTERNAL fsigpd_re,fsigpd_im,resfsigpd_im,resfsigpd_re
      external dqagi,dqagp,dqag,prerr
      common /sigcom/ anml,zbb,bbi,zaa,w,nml,numnml
      FLAG = .FALSE.
      i=dcmplx(0,1)
      za=XI+i*YI
      zaa=za
      zbb=zb
      bbi=bi
      l=0
      do 60 j=1,numn
         do 60 k=1,numm
            if(cdabs(anm(j,k)).LT.limtn) goto 60
            l=l+1
            anml(l)=anm(j,k)
            nml(l,1)=j-1
            nml(l,2)=k-1
 60   continue
      numnml=l
      w=zbb**2/4-zaa
      nlimit=10000
      ier=0
      if(dabs(dimag(w)).LT.limsingsm) then
         za=dble(za)+i*epsrel
         w=zbb**2/4-zaa
      endif
      alim=0.0
      CALL DQAGI(Fsigpd_re,alim,1,epsabs,epsrel,resr,abserr,neval,
     *     ier,nlimit,40000,last,iwork,work)
      if(ier.ne.0) goto 100
      CALL DQAGI(Fsigpd_im,alim,1,epsabs,epsrel,resi,abserr,neval,
     *     ier,nlimit,40000,last,iwork,work)
      if(ier.ne.0) goto 100
      u=resr
      v=resi
      if(dimag(zaa).LT.0.AND.dble(w).GT.0) then
         Alim=-1.0
         Blim=1.0
         CALL DQAG(resFsigpd_re,alim,blim,epsabs,epsrel,6,resr,abserr,
     *        neval,ier, nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         CALL DQAG(resFsigpd_im,alim,blim,epsabs,epsrel,6,resi,abserr,
     *        neval,ier,nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         u=u-resr
         v=v-resi
      else if(dimag(zaa).EQ.0.AND.dble(w).GT.0) then
         Alim=-1.0
         Blim=1.0
         CALL DQAG(resFsigpd_re,alim,blim,epsabs,epsrel,6,resr,abserr,
     *        neval,ier, nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         CALL DQAG(resFsigpd_im,alim,blim,epsabs,epsrel,6,resi,abserr,
     *        neval,ier,nlimit,40000,last,iwork,work)
         if(ier.ne.0) goto 100
         u=u-0.5*resr
         v=v-0.5*resi
      endif
      RETURN
*
 100  FLAG = .TRUE.
      call prerr(zaa,zbb,bbi)
      RETURN
*
      END
      
      double precision function Fsigpd_re(s)
      double complex Fsigpd
      external Fsigpd
      double precision s
      Fsigpd_re=dble(Fsigpd(s))
      return
      end function Fsigpd_re

      double precision function Fsigpd_im(s)
      double complex Fsigpd
      external Fsigpd
      double precision s
      Fsigpd_im=dimag(Fsigpd(s))
      return
      end function Fsigpd_im

      double complex function Fsigpd(s)
      double precision s,limsingsm,xbr,Jr0,Ji0,zbb,bbi,xbi,limtn
      integer numnml,nml(40,2),ierr,nz,j,k
      double complex z1,z2,zaa,Gm,weidGm,w,anml(40),zr
      common /sigcom/ anml,zbb,bbi,zaa,w,nml,numnml
c      common /sigcom/ anm,numn,numm,zbb,bbi,zaa,w
      external weidGm
      parameter (limsingsm = 1.0e-12,limtn=1.0e-20)
      zr=cdsqrt(4.0*w-2.0*s)
      z1=0.5*(zbb+zr)
      z2=0.5*(zbb-zr)
      Fsigpd=0.0d0
c      if ((cdabs(z1).LT.limsingsm).and.(cdabs(z2).LT.limsingsm)) then
c         res=0.0d0
c      else
         xbr=dsqrt(bbi*2.0d0*s)
         JR0=DBESJ0(XBR)
         do 50 j=1,numnml
            Gm=weidGm(z1,z2,nml(j,2))
            Fsigpd=Fsigpd+anml(j)*dexp(-s)*Jr0**2*
     *           Gm*s**(0.5d0*(nml(j,1)-1))
 50      continue
c      endif
c      Fsigpd=res
      return
      end function Fsigpd

      double precision function resFsigpd_re(s)
      double complex resFsigpd
      external resFsigpd
      double precision s
      resFsigpd_re=dble(resFsigpd(s))
      return
      end function resFsigpd_re

      double precision function resFsigpd_im(s)
      double complex resFsigpd
      external resFsigpd
      double precision s
      resFsigpd_im=dimag(resFsigpd(s))
      return
      end function resFsigpd_im

      double complex function resFsigpd(mu)
      double precision mu,xbr,xbi,Jr0,Ji0,zbb,bbi,sqrtpi,limtn,fnu
      integer numnml,nml(40,2),ierr,nz,j,k,nf,mf
      double complex zaa,i,w,xb,J0,anml(40),res
      external cbj0
      common /sigcom/ anml,zbb,bbi,zaa,w,nml,numnml
c      common /sigcom/ anm,numn,numm,zbb,bbi,zaa,w
      parameter (sqrtpi = 1.77245385090552,limtn=1.0e-20)
      xb=2.0d0*cdsqrt(bbi*(1-mu**2)*w)
      i=dcmplx(0,1.0D0)
      call cbj0(xb,J0)
      res=0.0d0
      do 70 j=1,numnml
         nf=nml(j,1)
         mf=nml(j,2)
         res=res+anml(j)*i*2**(0.5d0*(nf+3))*J0**2*sqrtpi*w**(nf*0.5d0)*
     *        (1-mu**2)**((nf-1)*0.5d0)*(mu*cdsqrt(w)+0.5D0*zbb)**mf*
     *        cdexp(-(mu*cdsqrt(w)+0.5D0*zbb)**2-2.0D0*(1.0D0-mu*mu)*w)
 70   continue
      resFsigpd=res
      return
      end function resFsigpd
