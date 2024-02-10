//CurvaEstiCAN.h
//Hongjun Li,
//

typedef struct{
	double x,y,z;
}Point_t, Vector_t;

typedef struct{
	int		  vNum;
    Point_t   *V;//所有的结点的坐标
//	int       *vLbl;//每个点的标识号
    //
	Vector_t   *Vn;
	Vector_t   *Dir1;
	Vector_t   *Dir2;
	double     *K1;
	double     *K2;
	double     K1Box[2], K2Box[2];//K1和 K2 的最小最大值
}PCDAnaOBJ_t;

//--------------input, initiate,output, free memory------
int InitTObjFile( PCDAnaOBJ_t *pOBJ);
int FreeTObjFile( PCDAnaOBJ_t *pOBJ);

// xyz file
int ReadFilePtsxyz( char *strFile,  PCDAnaOBJ_t *pPts );//读取xyz file
int SavePtsCurvaDirsK1K2( char *strFile, PCDAnaOBJ_t *pPts );//

//-----------------main algorithm -----------------------
int CreateGlobKdTree( void );//int neibNum
int FreeGlobKdTree(void) ; 
int NormalVectorCalculatingWithQFit( PCDAnaOBJ_t *pPtCdd, int numNeighbor);
int curvatureCalculatingWithCAN( PCDAnaOBJ_t *pPtCdd, int numNeighbor);