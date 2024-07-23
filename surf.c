#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#define XX   static 
#include "surf.h"

#define    OK_m(i,j) (   OK[((i) - 1)*K    + (j)])
#define   ULV_m(i,j) (  ULV[((i) - 1)*m    + (j)])
#define     A_m(i,j) (    A[((i) - 1)*NN   + (j)])
#define    A1_m(i,j) (   A1[((i) - 1)*NN   + (j)])
#define     B_m(i,j) (    B[((i) - 1)*m    + (j)])
#define    B1_m(i,j) (   B1[((i) - 1)*m    + (j)])
#define    TINY      1.e-7



void assign_z(int nsg,int b_seg[],int ii, int no,double *pz, int *ptest_2);
void b_n_seg(int *pnsg, int b_seg[]);
void loc_node(int ii, int nsg, int b_seg[], int *pres, int *pno);
void comp_bg(int ii, double x, double y, int *ptest);
void Wi(double x, double y, double b[], double e[], double W[]);
int set_eq(int k, double D,double Nq, double ok[], double f[]);
int search_bary( double xa, double ya, double l[], double e[]);
void areacord(double x1, double y1, double x2, double y2,
double x3, double y3, double x, double y, double *pl1,
double *pl2,double *pl3);
void LUbs(int n, int NN, double *A, int m, double *B, int *ch);
void max_D(double *pD);
void rinput_bg();
void seglen(double x1,double y1, double x2, double y2, double *pL);
void rinput();
void routput();
void open_write();
void tri_len(int i,double tL[]);
void cabl_len(int i, double *prL);
void out_cord(double x2,double y2,double x1,double y1,double del1,
double *px, double *py);
char* trimwhitespace(char *str);
void chkf(char *str);


int main()
{
int ii, i, j, k, check;
int elem_no;
int test_1, nsg, *b_seg,res,no, test_2;
char buf[80];
double D, *OK, *F, det,Nq;
double b[3], e[3];
double x,y,xk,yk,fk, Q[3], W[3];
double z;
double x_bg, y_bg;

   rinput_bg();

   rinput();

   if((OK = (double *)calloc(5*5, sizeof(double))) == NULL) 
   {printf("out of core for OK\n");exit(1);}

   if((F = (double *)calloc(5, sizeof(double))) == NULL) 
   {printf("out of core for F\n");exit(1);}

   if((b_seg = (int *)calloc(2*nn_bg, sizeof(int))) == NULL) 
   {printf("out of core for b_seg\n");exit(1);}

   b_n_seg(&nsg, b_seg); /* load the boundary segments into b_seg[] */

   max_D(&D);

   printf("Dmax = %f\n",D);

   Nq = 0.55*D;

   printf("Nq = %f\n", Nq);

   printf("Higher Nq value can improve the quality of mesh\n");
   printf("Enter Nq or RETURN to accept computed value\n");
   if (fgets(buf,80,stdin) != NULL) sscanf(buf,"%lf",&Nq);
   printf("Nq is set to: %f\n",Nq);

   /* enter the interpolation point */

   for(ii = 0; ii < nn; ii++)
   {


      x = cord[ii][0]; y = cord[ii][1];

      loc_node(ii,nsg,b_seg,&res, &no);/* if res = 1 then node is on b_seg */

      if(res == 1){
      assign_z(nsg,b_seg,ii,no,&z,&test_2);
      if(test_2 == 1){ cord[ii][2] = z; continue;}
      }

      comp_bg(ii,x,y,&test_1);/* function skips bg nodes */
      if(test_1 == 1) continue;
/*
     printf("Enter (x,y) \n");gets(buf);sscanf(buf,"%F %F",&x, &y);
*/

     elem_no=search_bary(x,y,b,e);/* locate (x,y) in bg mesh*/
   
     for(i = 0; i < 3; i++)
     {

      k = nod_bg[elem_no - 1][i] - 1;

      fk = cord_bg[k][2]; xk = cord_bg[k][0]; yk = cord_bg[k][1];

      check = set_eq(k,D,Nq,OK,F);   /*set up the minimizing problem equations*/ 

      if(check < 5){printf("less than 5 points taken\n"); exit(1);}

      det = LUdec(5,5,OK-1, 1, F-1,1.e-200);  /* solve equations */       
/*   
     for(j = 0; j < 5; j++) printf("%e\n",F[j]);printf("det = %e\n", det);    
*/
   
       Q[i] = fk + F[0]*(x-xk)+F[1]*(y-yk)+F[2]*pow((x-xk),2.)
          + F[3]*(x-xk)*(y-yk)+F[4]*pow((y-yk),2.);
       }

      Wi(x,y,b,e, W);

      z = 0.;
      for(i = 0; i < 3; i++)
      z += Q[i]*W[i];
      
      cord[ii][2] = z;

      /*printf(" node = %d x = %f y = %f z = %f\n",ii+1, x,y,z);*/
   }
   
   ndime = 3;
   routput();

   printf("***input.dat file has been generated***\n");

return 0;
}

/* function finds the neighbouring boundary nodes of the bg mesh
   close to the boundary node of the mesh to be mapped and checks
if the neighbouring z cords are constant then the same z value is
returned and test_2 = 1 */

void assign_z(int nsg,int b_seg[],int ii, int no,double *pz, int *ptest_2)
{
int i;
int n1,n2,n3, jk;
double x1,y1,z1,x2,y2,z2,x,y,L1,L2,z3,z_diff;
double tol = TINY;
   n1 = b_seg[2*(no-1)+0];  n2 = b_seg[2*(no-1)+1];
   x1 = cord_bg[n1-1][0];  y1 = cord_bg[n1-1][1]; z1 = cord_bg[n1-1][2];
   x2 = cord_bg[n2-1][0];  y2 = cord_bg[n2-1][1]; z2 = cord_bg[n2-1][2];
    x = cord[ii][0];        y = cord[ii][1];

   seglen(x1,y1,x,y,&L1);  seglen(x2,y2,x,y,&L2);

   if(L1 <= L2){
   jk = no-2; if(jk < 0) jk = nsg-1;
   n3 = b_seg[2*jk+0];  /* 1st node of previous segment */
   z3 = cord_bg[n3-1][2];
   if(fabs(z3-z1) < tol && fabs(z1-z2) < tol)
   {*ptest_2 = 1; *pz = z1; return;}
   }

   if(L1 > L2){
   jk = no; if(jk > (nsg-1)) jk = 0;
   n3 = b_seg[2*jk+1];  /* 2nd node of the next segment */
   z3 = cord_bg[n3-1][2];
   if(fabs(z3-z2) < tol && fabs(z2-z1) < tol)
   {*ptest_2 = 1; *pz = z2; return;}
   }

*ptest_2 = -1;
z2 = 1.e20;
return;
}

/* function to load and count the boundary segments */

void b_n_seg(int *pnsg, int b_seg[])
{
int i,j,k,l;
int nsg, knt;
int ns1, ns2;
int ne1, ne2, ne3;
double xs1,ys1,xs2,ys2,seg_L,xs,ys;
double xe1,ye1,xe2,ye2,xe3,ye3;
double a1,a2,a3,area,l1,l2,l3;

   nsg = 0;                /* the no. of boundary segments */

   for(i=0; i < ne_bg; i++)
   {
      for(j = 0; j < 3; j++)
      {
      k = j + 1; if(k > 2) k = 0;
      ns1 = nod_bg[i][j];
      ns2 = nod_bg[i][k];
      xs1 = cord_bg[ns1-1][0]; ys1 = cord_bg[ns1-1][1];
      xs2 = cord_bg[ns2-1][0]; ys2 = cord_bg[ns2-1][1];
      seglen(xs1,ys1,xs2,ys2,&seg_L);  /* length of side */
      out_cord(xs1,ys1,xs2,ys2,seg_L,&xs,&ys);

         /* this loop checks the side as common with any other element */
         for(l = 0; l < ne_bg; l++){          knt = 0;
         ne1 = nod_bg[l][0]; ne2 = nod_bg[l][1]; ne3 = nod_bg[l][2];
         xe1 = cord_bg[ne1-1][0]; ye1 = cord_bg[ne1-1][1];
         xe2 = cord_bg[ne2-1][0]; ye2 = cord_bg[ne2-1][1];
         xe3 = cord_bg[ne3-1][0]; ye3 = cord_bg[ne3-1][1];
         areacord(xe1,ye1,xe2,ye2,xe3,ye3, xs,ys,&l1,&l2,&l3);
         if(fabs(l1) <= 1.e-10) l1 = 0.;
         if(fabs(l2) <= 1.e-10) l2 = 0.;
         if(fabs(l3) <= 1.e-10) l3 = 0.;
         if(fabs(l1 - 1.) <= 1.e-10) l1 = 1.;
         if(fabs(l2 - 1.) <= 1.e-10) l2 = 1.;
         if(fabs(l3 - 1.) <= 1.e-10) l3 = 1.;
         if(l1 >= 0. && l2 >= 0. && l3 >= 0. && l1 <= 1. 
         && l2 <= 1. && l3 <= 1.)
         {knt++; break;}
         }
      if(knt == 0){         
      b_seg[2*nsg + 0] = ns1; 
      b_seg[2*nsg + 1] = ns2; nsg++;}
      }
  }
*pnsg = nsg;
return;
}

/* function locates the nodes from the mesh file which are on the
   boundary segments.  If the node is on the boundary then *pres = 1
   else *pres = -1 the function returns the no. of b_seg as *pno*/

 void loc_node(int ii, int nsg, int b_seg[], int *pres, int *pno)
{
int i;
int n1,n2;
int res;
double x1,y1,x2,y2,x,y,t,L;
double tol = TINY;
   x = cord[ii][0];   y = cord[ii][1];
   for (i = 0; i < nsg; i++)
   {
   n1 = b_seg[2*i+0];      n2 = b_seg[2*i+1];
   x1 = cord_bg[n1-1][0];  y1 = cord_bg[n1-1][1];
   x2 = cord_bg[n2-1][0];  y2 = cord_bg[n2-1][1];

   *pres = -1;

   if(fabs(x2-x1) < tol)
   {
   if((fabs(x-x1) < tol) && ((y >= y1 && y <= y2) || (y <= y1 && y >= y2)))
   {*pres = 1; *pno = i + 1; return;}
   *pres = -1; continue;
   }

   if(fabs(y2-y1) < tol)
   {
   if((fabs(y-y1) < tol) && ((x >= x1 && x <= x2) || (x <= x1 && x >= x2)))
   {*pres = 1; *pno = i +1; return;}
   *pres = -1; continue;
   }

   seglen(x1,y1,x2,y2,&L);
   t = ((x-x1)/(x2-x1))/L;

   if(fabs(t) < tol) {*pres = 1; *pno = i + 1; return;}
   if(fabs(t-1) < tol) {*pres = 1; *pno = i + 1; return;}
   if(t >= 0. && t <= 1.) {*pres = 1; *pno = i + 1; return;}
   }
*pres = -1;
*pno = -1;
return;
}
    
/* function compares the current node with all the bg nodes and
   if a match is found then test is returned +1 else -1.  If a
   match has been found then the z-cord of the bg node is 
   loaded into the z-cord of the node under reference */

void comp_bg(int ii, double x, double y, int *ptest)
{
int i,j, test;
double tol_x, tol_y;

   for(i = 0; i < nn_bg; i++)
   {
     tol_x = fabs(cord_bg[i][0] - x);
     tol_y = fabs(cord_bg[i][1] - y);
     
     if(tol_x < TINY && tol_y < TINY)
     {
     cord[ii][2] = cord_bg[i][2];
     *ptest = test = 1;            
     return;
     }
   }
*ptest = test = -1;
return;
}

/* function calculates the weight functions for the located triangular
   element and stores the Wi, Wj, Wk values in array W[] */

void Wi(double x, double y, double b[], double e[], double W[])
{
int i,j,k;
double tol =  1.e-10;

   for(i = 0; i < 3; i++)
   {
   j = i + 1; if(j > 2) j = 0;
   k = j + 1; if(k > 2) k = 0;
  
   W[i] = pow(b[i],2.)*(3. - 2.*b[i]) + 3.*(pow(b[i],2.)*b[j]*b[k])
          /(b[i]*b[j]+b[i]*b[k]+b[j]*b[k] + tol) *
          (b[j]*(pow(e[i],2.)+pow(e[k],2.)-pow(e[j],2.))/pow(e[k],2.) +
           b[k]*(pow(e[i],2.)+pow(e[j],2.)-pow(e[k],2.))/pow(e[j],2.));
   }
return;
}


/* function sets up the minimization problem equation */
int set_eq(int k, double D,double Nq, double ok[], double f[])
{
int i,j, res, knt; 
double X, Y, C[5], di, fi, fk, det;
double Rq, rho;  
double a1,a2,a3,a4,a5,a6;

       for(i = 0; i < 5; i++)        /* zero the ok and f matrices */
       {f[i] = 0.;for(j = 0; j < 5; j++) ok[5*i+j] = 0.;}

       fk = cord_bg[k][2];             /* form the coefficient matrix */
       knt = 0;        
       for(i = 0; i < nn_bg; i++)
       {
       if(i == k) continue;

       fi = cord_bg[i][2];
       X = cord_bg[i][0] - cord_bg[k][0];  Y = cord_bg[i][1] - cord_bg[k][1];

       seglen(cord_bg[i][0],cord_bg[i][1],cord_bg[k][0],cord_bg[k][1], &di);


       Rq = 0.5*D*sqrt(Nq/nn);  

       rho = (Rq-di)/(Rq*di); 
       if(rho < 0) {rho = 0.; knt++;}
       else rho = pow(rho,2.);
/*
       if((Rq-di) < 0){ rho = 0.; knt++;}
       else rho = 1./di;
*/
       if((nn-knt) < 5) return(nn-knt);



       C[0] = X*rho;  C[1] = Y*rho; C[2] = pow(X,2.)*rho;
       C[3] = X*Y*rho; C[4] = pow(Y,2.)*rho;
          for(j = 0; j < 5; j++)
          {
          a1 = ok[5*j+0] += C[j]*X; 
          a2 = ok[5*j+1] += C[j]*Y; 
          a3 = ok[5*j+2] += C[j]*pow(X,2.); 
          a4 = ok[5*j+3] += C[j]*X*Y; 
          a5 = ok[5*j+4] += C[j]*pow(Y,2.); 
          a6 = f[j] += C[j]*(fi - fk);        
          }          
       }
/*       for(i=0; i < 5; i++) 
       printf("%8.4e %8.4e %8.4e %8.4e %8.4e %8.4e\n",
       ok[5*i+0],ok[5*i+1],ok[5*i+2],ok[5*i+3],ok[5*i+4],f[i]);
*/

return(nn-knt);
}


/* function searches the element/elements of bg mesh, which contain the
  reference coordinates under reference and returns the area coords
  and the length of sides of the element*/

int search_bary( double xa, double ya, double l[], double e[])
{
int i,m,n1[3];
double x[3], y[3];
double ei, ej, ek;   
   for(i=0; i < ne_bg; i++)
   {
      for(m=0; m < 3; m++)
      {
       n1[m] = nod_bg[i][m];
       x[m] = cord_bg[n1[m]-1][0];
       y[m] = cord_bg[n1[m]-1][1];
      }
   
   areacord(x[0],y[0], x[1], y[1], x[2], y[2], xa,ya,
   &l[0], &l[1], &l[2]);
   if(fabs(l[0]) <= pow(10.,-4.))
   l[0] = 0.;
   if(fabs(l[1]) <= pow(10.,-4.))
   l[1] = 0.;
   if(fabs(l[2]) <= pow(10.,-4.))
   l[2] = 0.;

   if(fabs(l[0] - 1.) <= pow(10.,-4.))
   l[0] = 1.;
   if(fabs(l[1] - 1.) <= pow(10.,-4.))
   l[1] = 1.;
   if(fabs(l[2] - 1.) <= pow(10.,-4.))
   l[2] = 1.;

   if(l[0] >= 0. && l[1] >= 0. && l[2] >= 0. && l[0] <= 1. && l[1] <= 1.
   && l[2] <= 1.)
   {

   seglen(x[1],y[1],x[2],y[2],&ei);  e[0] = ei;/* load sides */
   seglen(x[2],y[2],x[0],y[0],&ej);  e[1] = ej;
   seglen(x[0],y[0],x[1],y[1],&ek);  e[2] = ek;
   return(i+1);

   }
   }

return(0);
}


void j_errorexit(char *s)

{

printf(s);
/*perror(s);*/
exit(1);

} /* end of funct. j_errorexit  */


double LU(int n,int NN,double *A,int *ch,double eps)

{
register int i,j,k;
int lmax;
double p,t,smax,*q;
double det;

q = (double *)calloc(n+5,sizeof(double));
for(i = 1;i <= n;i++)
    {
    p = 0.0;
    for(j = 1;j <= n;j++)
        {
        t = fabs(A_m(i,j));
        if(p < t)
            p = t;
        }
    if(p == 0)
j_errorexit("row of the coeff. matr. is zero in LU() increase Nq");
    q[i] = p;
    }
for(i = 1;i <= n;i++)
    ch[i] = i;
det = 1.0;
for(k = 1;k <= n;k++)
    {
    smax = 0.0;
    lmax = k;
    for(i = k;i <= n;i++)
        {
        p = A_m(i,k);
        for(j = 1;j < k;j++)
            p -= A_m(i,j)*A_m(j,k);
        if(fabs(q[ch[i]]) < TINY) q[ch[i]] += TINY;
        t = fabs(p)/q[ch[i]];
        if(smax < t)
            {
            smax = t;
            lmax = i;
            }
        A_m(i,k) = p;
        }         /* for -i- */
    if(smax < eps)
    j_errorexit("The coeff. matr. is almost singular increase Nq ");
    if(lmax != k)
        {
        j = ch[k];
        det *= -1.0;
        ch[k] = ch [lmax];
        ch[lmax] = j;
        for(j = 1;j <= n;j++)
            {
            p = A_m(k,j);
            A_m(k,j) = A_m(lmax,j);
            A_m(lmax,j) = p;
            }     /* for -j- */
        }         /* if */
    for(i = k+1;i <= n;i++)
        {
        if(fabs(A_m(k,k)) < TINY) A_m(k,k) += TINY;
        A_m(i,k) /= A_m(k,k);
        p = A_m(k,i);
        for(j = 1;j < k;j++)
            p -= A_m(k,j) * A_m(j,i);
        A_m(k,i) = p;
        }         /* for -i- */
    }             /* for -k- */
for(i = 1;i <= n;i++)
    det *= A_m(i,i);
free((char *)q);
return(det);

}

double LUdec(int n,int NN,double *A,int m,double *B,
                             double eps)

{
register int i,j,k;
double *X,p,det;
double *A1, *B1, *r, sdp;
double tol = 1.e-4, tol_s;
int *ch;

   if((A1 = (double *)calloc(n*n + 5, sizeof(double))) == NULL) 
   {printf("out of core for A1\n");exit(1);}

   if((B1 = (double *)calloc(n*m + 5, sizeof(double))) == NULL) 
   {printf("out of core for B1\n");exit(1);}

   if((r = (double *)calloc(n + 5, sizeof(double))) == NULL) 
   {printf("out of core for r\n");exit(1);}

   for(i = 1; i <= n*n; i++) A1[i] = A[i];/* save A in A1 */
   for(i = 1; i <= n*m; i++)   B1[i] = B[i]; /* save B in B1 */

ch = (int *)calloc(n+5,sizeof(int));

det = LU(n,NN,A,ch,eps);

LUbs(n, NN, A, m, B, ch); /* first approx. solution */

   /* solution improvement starts here */
   tol_s = 0.;
   for(k = 1; k <= m; k++)
   {
    do{
    for(i = 1; i <= n; i++)  /* calculate the residuals */
    {
    sdp = -B1_m(i,k);
    for(j = 1; j <= n; j++) sdp += A1_m(i,j)*B_m(j,k);
    r[i] = sdp;
    tol_s += pow(sdp,2.);
    }   
    LUbs(n, NN, A, m, r, ch);
    for(i = 1; i <= n*m; i++) B_m(i,k) -= r[i];
    tol_s = sqrt(tol_s/((double) n));/* R.M.S. residuals*/
    } while( tol_s > tol);
   }

free((char *)ch);
return(det);

}  /* end of funct LUdec  */

/* function to perform back subsitution */
void LUbs(int n, int NN, double *A, int m, double *B, int *ch)
{
int i, j, k;
double p, *X; 

X = (double *)calloc(n+5,sizeof(double));

for(k = 1; k <= m ;k++)
    {
    for(i = 1;i <= n;i++)
        X[i] = B_m(ch[i],k);
    for(i = 1;i <= n;i++)
        {
        p = X[i];
        for(j = 1;j <= i-1;j++)
            p -= A_m(i,j)*X[j];
        X[i] = p;
        }
    for(i = n;i >= 1;i--)
        {
        p = X[i];
        for(j = i+1;j <= n;j++)
            p -= A_m(i,j)*X[j];
        if(fabs(A_m(i,i)) < TINY) A_m(i,i) += TINY;
        X[i] = p/A_m(i,i);
        }
    for(i = 1;i <= n;i++)
        B_m(i,k) = X[i];
    }                 /* for -k- */
free((char *)X);
return;
}

/* function calculates the maximum distance between any two nodes */
void max_D(double *pD)
{
int i, j;
double xi,yi,xj, yj, D, Dmax;
   Dmax = 0.;
   for(i = 0; i < nn_bg; i++)
   {
   xi = cord_bg[i][0]; yi = cord_bg[i][1];
       for(j = 0; j < nn; j++)
       {
       if(i == j) continue;
       xj = cord_bg[j][0]; yj = cord_bg[j][1];
       seglen(xi,yi,xj,yj, &D);
       Dmax = (Dmax > D) ? Dmax : D;
       }
   }
*pD = Dmax;
return;
}



void rinput_bg()
{
FILE *fp1;
char buf[256];
int i,j,k,n,nq,nodq[1][5],ngs[1],nsnod[1][8];
float val,x1,y1,z1;

    printf("Enter the (3D) coarse mesh file name (without extension)\n");
    chkf(fgets(buf,256,stdin));
    strcpy(buf,trimwhitespace(buf));
    strcat(buf,".dat");
    fp1  = fopen(buf,"r");

   chkf(fgets(buf,256,fp1));
   fputs(buf,stdout);
   chkf(fgets(buf,256,fp1));
   sscanf(buf,"%F%F",&dt,&damp);
   printf("%f %f \n",dt,damp);
   chkf(fgets(buf,256,fp1));
   sscanf(buf,"%d %d %d %d %d %d %d %d %d %d %d",&nn_bg,&ne_bg,&nb,&nm,&nl,
   &ndime, &nstr, &nbe, &kmax,&nmax, &iln);

   chkf(fgets(buf,256,fp1));
   sscanf(buf,"%d",&nq);



   for(i = 0; i < nn_bg; i++)
   {
   chkf(fgets(buf,256,fp1));
   sscanf(buf,"%d %f %f %f",&n, &x1, &y1,&z1);
   val = cord_bg[i][0] = x1;
   val = cord_bg[i][1] = y1;
   val = cord_bg[i][2] = z1;
   }



   if(ne_bg > 0)
   {
      for(i = 0; i < ne_bg; i++)
      {
      chkf(fgets(buf,256,fp1));
      sscanf(buf,"%d %d %d %d %d",&n,&nod_bg[i][0],&nod_bg[i][1],
      &nod_bg[i][2],&nod_bg[i][3]);
      }
   }



   if(nl > 0)
   {
      for(i = 0; i < nl;i++)
      {
      chkf(fgets(buf,256,fp1));
      sscanf(buf,"%d %d %d %d",&n,&nodl[i][0],&nodl[i][1],
      &nodl[i][2]);
      }
   }



   if(nq > 0)
   {
      for(i = 0; i < nq; i++)
      {
      chkf(fgets(buf,256,fp1));
      sscanf(buf,"%d %d %d %d %d %d",&n,&nodq[i][0],&nodq[i][1],
      &nodq[i][2],&nodq[i][3],&nodq[i][4]);
      }
   }



   for(i = 0; i < nm; i++)
   {
   chkf(fgets(buf,256,fp1));
   sscanf(buf,"%d %d %F %E %E %E %E",&n,&mat[i],&amat[i][0],&amat[i][1],
   &amat[i][2],&amat[i][3],&amat[i][4]);
   }



   for(i = 0; i < nb; i++)
   {
   chkf(fgets(buf,256,fp1));
   sscanf(buf,"%d %d %d %d",&nco[i][0],&nco[i][1],&nco[i][2],&nco[i][3]);
   }



   if(iln > 0)
   {
   iln = 0;
more:
   iln++;
   i = iln - 1;

   chkf(fgets(buf,256,fp1));
   sscanf(buf,"%d %F %F %F",&nc[i],&aload[i][0],&aload[i][1],&aload[i][2]);

   if(nc[i] < nn)
   goto more;
   }

   if(nstr > 0)
   {
   for(i = 0; i < nstr; i++)
   {
   chkf(fgets(buf,256,fp1));
   sscanf(buf,"%d",&ngs[i]);
   chkf(fgets(buf,256,fp1));
   sscanf(buf,"%d %d %d %d %d %d %d %d",&nsnod[i][0],&nsnod[i][1],
   &nsnod[i][2],&nsnod[i][3],&nsnod[i][4],&nsnod[i][5],&nsnod[i][6],
   &nsnod[i][7]);
   }
   }
return;
}



/*--------------------- lib ---------------------------*/
      
/* Function areacord returns area cord of a trian elem */
void areacord(double x1, double y1, double x2, double y2,
double x3, double y3, double x, double y, double *pl1,
double *pl2,double *pl3)
{
double area,a1,a2,a3;

      area=(x2*y3-x3*y2)-(x1*y3-x3*y1)+(x1*y2-x2*y1);

      a1=(x2*y3-x3*y2)+x*(y2-y3)+y*(x3-x2);
      a2=(x3*y1-x1*y3)+x*(y3-y1)+y*(x1-x3);
      a3=(x1*y2-x2*y1)+x*(y1-y2)+y*(x2-x1);

      *pl1=a1/area;
      *pl2=a2/area;
      *pl3=a3/area;
return;
}

      
/* Function to calculate the length of a segment */
void seglen(double x1,double y1, double x2, double y2, double *pL)
{
double L;
L=sqrt(pow(x2-x1,2.)+pow(y2-y1,2.));
*pL = L;
return;
}


void rinput()
{
FILE *fp1;
char buf[256];
int i,j,k,n,nodq[1][5],ngs[1],nsnod[1][8];
float val,x1,y1,z1;

    printf("Enter the (2D) refined mesh file name (without extension)\n");
    chkf(fgets(buf,256,stdin));
    strcpy(buf,trimwhitespace(buf));
    strcat(buf,".dat");
    fp1  = fopen(buf,"r");
    
    open_write();

   fgets(buf,256,fp1);
   fputs(buf,stdout);
   fputs(buf, fp0);
   fgets(buf,256,fp1);
   fputs(buf,fp0);
   sscanf(buf,"%F%F",&dt,&damp);
   printf("%f %f \n",dt,damp);
   fgets(buf,256,fp1);
   sscanf(buf,"%d %d %d %d %d %d %d %d %d %d %d",&nn,&ne,&nb,&nm,&nl,
   &ndime, &nstr, &nbe, &kmax,&nmax, &iln);

   fgets(buf,256,fp1);
   sscanf(buf,"%d",&nq);



   for(i = 0; i < nn; i++)
   {
   fgets(buf,256,fp1);
   sscanf(buf,"%d %f %f %f",&n, &x1, &y1,&z1);
   val = cord[i][0] = x1;
   val = cord[i][1] = y1;
   val = cord[i][2] = z1;
   }



   if(ne > 0)
   {
      for(i = 0; i < ne; i++)
      {
      fgets(buf,256,fp1);
      sscanf(buf,"%d %d %d %d %d",&n,&nod[i][0],&nod[i][1],
      &nod[i][2],&nod[i][3]);
      }
   }



   if(nl > 0)
   {
      for(i = 0; i < nl;i++)
      {
      fgets(buf,256,fp1);
      sscanf(buf,"%d %d %d %d",&n,&nodl[i][0],&nodl[i][1],
      &nodl[i][2]);
      }
   }



   if(nq > 0)
   {
      for(i = 0; i < nq; i++)
      {
      fgets(buf,256,fp1);
      sscanf(buf,"%d %d %d %d %d %d",&n,&nodq[i][0],&nodq[i][1],
      &nodq[i][2],&nodq[i][3],&nodq[i][4]);
      }
   }



   for(i = 0; i < nm; i++)
   {
   fgets(buf,256,fp1);
   sscanf(buf,"%d %d %F %E %E %E %E",&n,&mat[i],&amat[i][0],&amat[i][1],
   &amat[i][2],&amat[i][3],&amat[i][4]);
   }



   for(i = 0; i < nb; i++)
   {
   fgets(buf,256,fp1);
   sscanf(buf,"%d %d %d %d",&nco[i][0],&nco[i][1],&nco[i][2],&nco[i][3]);
   }



   if(iln > 0)
   {
   iln = 0;
more:
   iln++;
   i = iln - 1;

   fgets(buf,256,fp1);
   sscanf(buf,"%d %F %F %F",&nc[i],&aload[i][0],&aload[i][1],&aload[i][2]);

   if(nc[i] < nn)
   goto more;
   }

   if(nstr > 0)
   {
   for(i = 0; i < nstr; i++)
   {
   fgets(buf,256,fp1);
   sscanf(buf,"%d",&ngs[i]);
   fgets(buf,256,fp1);
   sscanf(buf,"%d %d %d %d %d %d %d %d",&nsnod[i][0],&nsnod[i][1],
   &nsnod[i][2],&nsnod[i][3],&nsnod[i][4],&nsnod[i][5],&nsnod[i][6],
   &nsnod[i][7]);
   }
   }
fclose(fp1);
return;
}



/* function to write the data to the primary data file */

void routput()
{
int i,j,k,n,nodq[1][5],ngs[1],nsnod[1][8];
double tL[3],rL;


   fprintf(fp0,"%6d%6d%5d%5d%5d%5d%5d%5d%5d%5d%5d\n",nn,ne,nb,nm,nl,
   ndime, nstr, nbe, kmax,nmax, iln);

   fprintf(fp0,"%5d\n",nq);



   for(i = 0; i < nn; i++)
   {
   fprintf(fp0,"%5d%15.5f%15.5f%15.5f\n",i+1, cord[i][0], cord[i][1], 
   cord[i][2]);
   }



   if(ne > 0)
   {
      for(i = 0; i < ne; i++)
      {
      tri_len(i, tL);
      fprintf(fp0,"%5d%5d%5d%5d%5d%15.8f%15.8f%15.8f\n",i+1,nod[i][0],nod[i][1],
      nod[i][2],nod[i][3],tL[0],tL[1],tL[2]);
      }
   }



   if(nl > 0)
   {
      for(i = 0; i < nl;i++)
      {
      cabl_len(i, &rL);
      fprintf(fp0,"%5d%5d%5d%5d%15.8f\n",i+1,nodl[i][0],nodl[i][1],
      nodl[i][2],rL);
      }
   }



   if(nq > 0)
   {
      for(i = 0; i < nq; i++)
      {
      fprintf(fp0,"%5d%5d%5d%5d%5d%5d",i+1,nodq[i][0],nodq[i][1],
      nodq[i][2],nodq[i][3],nodq[i][4]);
      }
   }



   for(i = 0; i < nm; i++)
   {
   fprintf(fp0,"%5d%5d%12.5e%12.5e%12.5e%12.5e%12.5e\n",i+1,mat[i],amat[i][0],
   amat[i][1],amat[i][2],amat[i][3],amat[i][4]);
   }



   for(i = 0; i < nb; i++)
   {
   fprintf(fp0,"%5d%5d%5d%5d\n",nco[i][0],nco[i][1],nco[i][2],nco[i][3]);
   }



   if(iln > 0)
   {
   iln = 0;
more:
   iln++;
   i = iln - 1;

   fprintf(fp0,"%5d%15.5f%15.5f%15.5f\n",nc[i],aload[i][0],aload[i][1],
   aload[i][2]);

   if(nc[i] < nn)
   goto more;
   }

   if(nstr > 0)
   {
   for(i = 0; i < nstr; i++)
   {
   fprintf(fp0,"%5d\n",ngs[i]);
   fprintf(fp0,"%5d%5d%5d%5d%5d%5d%5d%5d\n",nsnod[i][0],nsnod[i][1],
   nsnod[i][2],nsnod[i][3],nsnod[i][4],nsnod[i][5],nsnod[i][6],
   nsnod[i][7]);
   }
   }
fclose(fp0);
return;
}

/* Function to read input.dat file format*/
void open_write()
{
fp0=fopen("input.dat","w");
	if(fp0==(FILE *)NULL)
	{
	printf("cannot open input.dat\n");
	exit(1);
	}
return;
}

/* function to calculate the lengths of the three sides of each element
for storage in tlen[] */
void tri_len(int i,double tL[])
{
int j,k;
double x1,y1,x2,y2,L;

   k = 0;
      for(j = 0; j < 3; j++)
      {
      k = j + 1;
      if(k > 2)
      k = 0;
      x1 = cord[nod[i][j] - 1][0];
      y1 = cord[nod[i][j] - 1][1];
      x2 = cord[nod[i][k] - 1][0];
      y2 = cord[nod[i][k] - 1][1];
      seglen(x1,y1,x2,y2,&L);
      tL[j] = L;
      k++;
      }
return;
}

/* function to calculate the length of the cable segments */
void cabl_len(int i, double *prL)
{
double x1,y1,x2,y2,L; 

   x1 = cord[nodl[i][0] - 1][0];
   y1 = cord[nodl[i][0] - 1][1];
   x2 = cord[nodl[i][1] - 1][0];
   y2 = cord[nodl[i][1] - 1][1];
   seglen(x1,y1,x2,y2,&L);
   *prL = L;
return;
}

/* Function to calculate the a perpendicular cord at right 
of segment and at a distance .001*L from the segment */
void out_cord(double x2,double y2,double x1,double y1,double del1,
double *px, double *py)
{
double x3,y3,L,s;
seglen(x2,y2,x1,y1,&L);
s=.001;
x3=(x1+x2)/2.;
y3=(y1+y2)/2.;
*px=x3-s*(y2-y1);
*py=y3+s*(x2-x1);
return;
}

char* trimwhitespace(char *str)
{
  char *end;

  /* Trim leading space */
  while(isspace((unsigned char)*str)) str++;

  if(*str == 0)  /* All spaces? */
    return str;

  /* Trim trailing space */
  end = str + strlen(str) - 1;
  while(end > str && isspace((unsigned char)*end)) end--;

  /* Write new null terminator character */
  end[1] = '\0';

  return str;
}
void chkf(char *str){
   if( str == NULL){
     printf("line not read\n");
     exit(1);
   }
}
