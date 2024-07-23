
/*                         RINPUT.H                               */

/*               Header file for RINPUT.C                         */



#define    MN    4000
#define    ME    5000
#define    MB     500
#define    ML      2
#define    MM      5
#define    MDIME   3
#define    MSTR    2
#define    MBE    200
#define    MLN    500
#define    MGS     50


XX FILE *fp0;

XX int      nn, ne, nb, nm, nl, ndime, nstr, nbe, kmax, nmax, iln;
XX int      nq;
XX int      nod[ME][4], nodl[ML][3], mat[MM], nco[MB][MDIME+1], nc[MLN];
XX double   dt, damp;
XX double   cord[MN][MDIME], amat[MM][5], aload[MLN][MDIME];

XX int nn_bg,ne_bg, nod_bg[ME][4];
XX double cord_bg[MN][MDIME];


void j_errorexit(char *s);

double LU(int n,int NN,double *A,int *ch,double eps);

double LUdec(int n,int NN,double *A,int m,double *B,double eps);


