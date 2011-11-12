*** soft approximation for the matrix element qq->W+4g
*format fortran;
CF u,v,w,x,y,z,XC3,XC4,XC5,XC6,XC7;
CF tr,tt,t,t1,t2,t3,t4,tstring,Fstring,F,F1,trace,scr;
S ans,V,xn,q,g,qb,a0,CC4,CC5,CC6,xn2,one,perm2,perm3,perm4,b0,e;
S bit2,bit3,bit4,e,f,b,q,qb,Qb,Q,nf,[Nc*d_(ix2,ix3)*d_(ix1,ix4)]
S [Nc*t(C5,ix4,ix1)*d_(ix2,ix3)],[Nc*t(C5,ix2,ix3)*d_(ix4,ix1)]
I C1,C2,C3,C4,C5,C6,C7,C8,C9,Cx1,Cx2,Cx3,Cx4;
I C3a,C4a,C5a,C1a,C2a;
S cf1,cf2,cf3,cf4,cf1on2,cf2on2,cf3on2,cf4on2;
s A123,A132,A213,A231,A312,A321;
s Ac123,Ac132,Ac213,Ac231,Ac312,Ac321;
s [p1k],[p2k],[p3k],[p4k],[p5k],[pqk],[pbqk];
s s1k,s2k,s3k,sqk,sbqk;
v p1,p2,p3,p4,p5,pq,pbq;
i mu;
s xxx;
Autodeclare CF A,B,Bf;
Autodeclare I Cx=V,Cy=V;
Autodeclare I i=xn,ix=xn,jx=xn,iq=xn,jq=xn,it=xn,jt=xn,iy=xn,jy=xn,iz=xn,iw=xn;
Autodeclare I ia=xn,ib=xn;
Autodeclare S c;
I j1,j2,j3,j4;
CF AR,AL,ALO,cd,cd1,Pr;
cf Attqqg,Attgqq,Atgtqq,Attqgq;
.global

* glossary
********************************************************************
* Attqqg1g2 = A(1)
* Attqqg2g1 = A(2)

* Attg1g2qq = A(3)
* Attg2g1qq = A(4)

* Attg1qqg2 = A(5)
* Attg2qqg1 = A(6)

* Atg1g2tqq = A(7)
* Atg2g1tqq = A(8)

* Atg1tqg2q = A(9)
* Atg2tqg1q = A(10)

* Attqg1g2q = A(11)
* Attqg2g1q = A(12)
*********************************************************************

* all 1/xn terms have to be multiplied by (-1)

g expr =
(
 cd(it,jq)*cd(iq,jx1)*cd(ix1,jx2)*cd(ix2,jt)*A(1)
+cd(it,jq)*cd(iq,jx2)*cd(ix2,jx1)*cd(ix1,jt)*A(2)
+cd(iq,jt)*cd(it,jx1)*cd(ix1,jx2)*cd(ix2,jq)*A(3)
+cd(iq,jt)*cd(it,jx2)*cd(ix2,jx1)*cd(ix1,jq)*A(4)
+cd(it,jx1)*cd(ix1,jq)*cd(iq,jx2)*cd(ix2,jt)*A(5)
+cd(it,jx2)*cd(ix2,jq)*cd(iq,jx1)*cd(ix1,jt)*A(6)
+1/xn*cd(it,jx2)*cd(ix2,jx1)*cd(ix1,jt)*cd(iq,jq)*A(7)
+1/xn*cd(it,jx1)*cd(ix1,jx2)*cd(ix2,jt)*cd(iq,jq)*A(8)
+1/xn*cd(it,jx1)*cd(ix1,jt)*cd(iq,jx2)*cd(ix2,jq)*A(9)
+1/xn*cd(it,jx2)*cd(ix2,jt)*cd(iq,jx1)*cd(ix1,jq)*A(10)
+1/xn*cd(iq,jx2)*cd(ix2,jx1)*cd(ix1,jq)*cd(it,jt)*A(11)
+1/xn*cd(iq,jx1)*cd(ix1,jx2)*cd(ix2,jq)*cd(it,jt)*A(12)
)*
Pr(iy1,jy1,ix1,jx1)*Pr(iy2,jy2,ix2,jx2)*
(
 cd1(it,jq)*cd1(iq,jy1)*cd1(iy1,jy2)*cd1(iy2,jt)*Ac(1)
+cd1(it,jq)*cd1(iq,jy2)*cd1(iy2,jy1)*cd1(iy1,jt)*Ac(2)
+cd1(iq,jt)*cd1(it,jy1)*cd1(iy1,jy2)*cd1(iy2,jq)*Ac(3)
+cd1(iq,jt)*cd1(it,jy2)*cd1(iy2,jy1)*cd1(iy1,jq)*Ac(4)
+cd1(it,jy1)*cd1(iy1,jq)*cd1(iq,jy2)*cd1(iy2,jt)*Ac(5)
+cd1(it,jy2)*cd1(iy2,jq)*cd1(iq,jy1)*cd1(iy1,jt)*Ac(6)
+1/xn*cd1(it,jy2)*cd1(iy2,jy1)*cd1(iy1,jt)*cd1(iq,jq)*Ac(7)
+1/xn*cd1(it,jy1)*cd1(iy1,jy2)*cd1(iy2,jt)*cd1(iq,jq)*Ac(8)
+1/xn*cd1(it,jy1)*cd1(iy1,jt)*cd1(iq,jy2)*cd1(iy2,jq)*Ac(9)
+1/xn*cd1(it,jy2)*cd1(iy2,jt)*cd1(iq,jy1)*cd1(iy1,jq)*Ac(10)
+1/xn*cd1(iq,jy2)*cd1(iy2,jy1)*cd1(iy1,jq)*cd1(it,jt)*Ac(11)
+1/xn*cd1(iq,jy1)*cd1(iy1,jy2)*cd1(iy2,jq)*cd1(it,jt)*Ac(12)
);


id cd1(i1?,i2?) = cd(i2,i1);
.sort



id Pr(ix2?,jx2?,ix1?,jx1?) =
cd(ix2,ix1)*cd(jx1,jx2)-1/xn*cd(ix2,jx2)*cd(jx1,ix1);
.sort

repeat, id cd(i1?,i2?)*cd(i2?,i3?) = cd(i1,i3);
.sort


id cd(i1?,i1?) = d_(i1,i1);
.sort

id xn = 3;
id 1/xn = 1/3;
.sort


b A,Ac;
Format doublefortran;

print;
.end







Id,t(C5,ix4,ix1)*d_(ix2,ix3)*xn^-1=cf1on2;
Id,t(C5,ix2,ix3)*d_(ix1,ix4)*xn^-1=cf3on2;
Id,t(C5,ix2,ix1)*d_(ix3,ix4)*xn^-2=cf2on2;
Id,t(C5,ix4,ix3)*d_(ix1,ix2)*xn^-2=cf4on2;
Id,t(C5,ix4,ix1)*d_(ix2,ix3)*xn^1=cf1;
Id,t(C5,ix2,ix3)*d_(ix1,ix4)*xn^1=cf3;
Id,t(C5,ix2,ix1)*d_(ix3,ix4)=cf2;
Id,t(C5,ix4,ix3)*d_(ix1,ix2)=cf4;

Id,cf1on2=cf1/xn^2;
Id,cf2on2=cf2/xn^2;
Id,cf3on2=cf3/xn^2;
Id,cf4on2=cf4/xn^2;
***cf1,cf3 leading
***cf2,cf4 QED
.sort

id pq.pq=0;
id pbq.pbq=0;
id p1.p1 = 0;
id p2.p2 = 0;
id p3.p3 = 0;
id I^2 = -1;
.sort

*id p1?.p2?=scr(p1,p2);
*.sort

id 1/[p1k] = 1/s1k;
id 1/[p2k] = 1/s2k;
id 1/[p3k] = 1/s3k;
id 1/[pqk] = 1/sqk;
id 1/[pbqk] = 1/sbqk;


*if (count(xn,1) !=5) discard;
*.sort


*b p1,p2,pq,pbq,p3,[pqk],[pbqk],[p1k],[p2k],[p3k];
b p1,p2,pq,pbq,p3,s1k,s2k,s3k,sqk,sbqk;
*B A123,A132,A213,A231,A312,A321;

Print;
.end
