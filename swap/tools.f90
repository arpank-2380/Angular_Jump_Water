 module tools
 implicit none
 contains

 subroutine dotp(a,b,p)
 real a(3),b(3),p
 integer i
 p=0.0d0
 do i=1,3
   p=p+a(i)*b(i)
 enddo
 return
 end subroutine
                                                                                                                                                        
 subroutine crossp(a,b,c)
 real  a(3),b(3),c(3)
 integer i
 do i=1,3
  c(i)=0.0d0
 enddo
 c(1)=c(1)+a(2)*b(3)-a(3)*b(2)
 c(2)=c(2)+a(3)*b(1)-a(1)*b(3)
 c(3)=c(3)+a(1)*b(2)-a(2)*b(1)
 return
 end subroutine
                                                                                                                                                        
 function vmod(a)
 real vmod
 real a(3), sum
 integer i
 sum=0.0d0
 do i=1,3
  sum=sum+a(i)*a(i)
 enddo
 vmod=sqrt(sum)
 return
 end

 function angle(a,b)
 real angle
 real a(3), b(3),prod, costheta
 call dotp(a,b,prod)
 costheta=prod/(vmod(a)*vmod(b))
 if ( costheta .gt. 1.0d0) then
    costheta = 1.0d0
 else if ( costheta .lt. -1.0d0) then
    costheta = -1.0d0
 endif
 angle=acos(costheta)
 return
 end
                                                                                                                                                        
 subroutine projection_on_plane(a,b,c,v,psi,phi,theta)
 real a(3),b(3),c(3),v(3),proj_vec(3),unorm(3),comp_norm(3),rab(3), rac(3),cp(3),bisec(3)
 integer i
 real psi, phi, theta, cpmod, sinphi, dp, alpha, beta
 do i=1,3
   rab(i)=b(i)-a(i)
   rac(i)=c(i)-a(i)
 enddo
 psi=(180./3.14159)*angle(rab,rac)
 call angle_bisec(rab,rac,bisec)
 call crossp(rab,rac,cp)
 cpmod=vmod(cp)
 do i=1,3
 unorm(i)=cp(i)/cpmod
 enddo
 call dotp(unorm,v,dp)
 sinphi=abs(dp/vmod(v))
 if ( sinphi .gt. 1.0d0 ) then
    sinphi = 1.0d0
 endif
 phi=(180./3.14159)*asin(sinphi)
 do i=1,3
 comp_norm(i)=dp*unorm(i)
 proj_vec(i)=v(i)-comp_norm(i)
 enddo
  alpha=(180./3.14159)*angle(rab,bisec)
  beta=(180./3.14159)*angle(rab,proj_vec)
  theta=(180./3.14159)*angle(bisec,proj_vec)
 if((beta.gt.alpha).and.(beta.gt.theta)) then
       theta=theta
  else
       theta=-theta
  endif
 return
 end subroutine
                                                                                                                                                        
                                                                                                                                                        
  subroutine angle_bisec(p,q,bisec)
  real p(3),q(3),bisec(3)
  integer i  
     do i=1,3
       bisec(i)=p(i)/vmod(p)+q(i)/vmod(q)
     enddo
  return
  end

 end module tools
