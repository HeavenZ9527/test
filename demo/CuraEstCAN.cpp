//PCDAnalysis.cpp
//
/*
*
*
*/
//#include"StdAfx.h"
#include<stdio.h>
#include<math.h>
#include<fstream>
#include<time.h>


#include "CurvaEstiCAN.h"

#include "ANN/ANN.h"
#include "SolveEquation/svdSolveEquation.h"



char buf[100];//屏幕输出


#define PI 3.14159265358979

			  //1 plane 2 柱面, 3 球面,4 马鞍面  其他
			  //-----------------------------------------

PCDAnaOBJ_t   *pScanD = NULL;
ANNkd_tree    *GkdTree;
ANNpointArray dataPts;

//--------------------------------------------------------------
int main()
{
	double start, stop, durationTime;
	start = clock();
	int numNeighbor = 2400;
	char strFile[500] = "TorusSurf.xyz";
	char svStrFile[500] = "TorusSurfDirK1K2.xyz";
	pScanD = new PCDAnaOBJ_t;
	InitTObjFile(pScanD);
	int opend = ReadFilePtsxyz(strFile, pScanD);
	if (opend == 0) {
		printf("!! Error: cannot open:%s\n", strFile);
		return 0;
	}
	CreateGlobKdTree();//int neibNum
	pScanD->Vn = new Vector_t[pScanD->vNum];
	NormalVectorCalculatingWithQFit(pScanD, numNeighbor);
	pScanD->Dir1 = new Vector_t[pScanD->vNum];
	pScanD->Dir2 = new Vector_t[pScanD->vNum];
	pScanD->K1 = new double[pScanD->vNum];
	pScanD->K2 = new double[pScanD->vNum];
	curvatureCalculatingWithCAN(pScanD, numNeighbor);

	SavePtsCurvaDirsK1K2(svStrFile, pScanD);

	FreeTObjFile(pScanD);
	stop = clock();

	durationTime = ((double)(stop - start)) / CLK_TCK;
	printf("程序耗时：%f ms. \n", durationTime);
	system("pause");
	return 1;
}


int InitTObjFile(PCDAnaOBJ_t *pOBJ)
{
	pOBJ->vNum = 0;
	pOBJ->V = NULL;//所有的结点的坐标
				   //	pOBJ->vLbl = NULL;//每个点的标识号
				   //	pOBJ->Box[0] = pOBJ->Box[1] = 0;
				   //	pOBJ->Box[2] = pOBJ->Box[3] = 0;
				   //	pOBJ->Box[4] = pOBJ->Box[5] = 0;
				   //
	pOBJ->Vn = NULL;
	pOBJ->Dir1 = NULL;
	pOBJ->Dir2 = NULL;
	pOBJ->K1 = NULL;
	pOBJ->K2 = NULL;
	//
	return 1;
}

int FreeTObjFile(PCDAnaOBJ_t *pOBJ)
{
	if (pOBJ->V)
		delete[]pOBJ->V;//所有的结点的坐标
	if (pOBJ->Vn)
		delete[]pOBJ->Vn;
	if (pOBJ->Dir1)
		delete[]pOBJ->Dir1;
	if (pOBJ->Dir2)
		delete[]pOBJ->Dir2;
	if (pOBJ->K1)
		delete[]pOBJ->K1;
	if (pOBJ->K2)
		delete[]pOBJ->K2;
	return 1;

}


//# respond
int ReadFilePtsxyz(char *strFile, PCDAnaOBJ_t *pPts)
{//读取xyz file
	FILE *filep;
	char ch[300];
	int i;
	float x, y, z;

	printf("<========%s\n", strFile);
	pPts->vNum = 0;

	filep = fopen(strFile, "r");
	if (filep == NULL) {
		printf("Error: can not open\n");
		return 0;
	}
	//-----------------------------------------------------------
	while (fgets(ch, 300, filep) != NULL) {
		if (ch[0] == '#' || strlen(ch) < 2) {
			continue;
		}
		else
			pPts->vNum++;
	}

	if (pPts->vNum == 0) {
		printf("Error: with v as begin for each point.\n");
		fclose(filep);
		return 0;
	}
	//------------------------------------------------------------	
	pPts->V = new Point_t[pPts->vNum];

	fseek(filep, 0, 0);
	i = 0;
	while (fgets(ch, 300, filep) != NULL) {
		if (ch[0] == '#' || strlen(ch) < 2)
			continue;
		else {
			sscanf(ch, "%g %g %g", &x, &y, &z);
			pPts->V[i].x = x;
			pPts->V[i].y = y;
			pPts->V[i].z = z;
			i++;
		}
	}//while
	fclose(filep);
	printf("  %d points \n", pPts->vNum);
	return pPts->vNum;
}

int SavePtsCurvaDirsK1K2(char *strFile, PCDAnaOBJ_t *pPts)
{//
	int     i;
	FILE	*filep;

	printf("========>%s ", strFile);
	filep = fopen(strFile, "wt");
	if (filep == NULL)
	{
		printf("can not be opened for save!\n");
		exit(0);
	}
	//-------------------
	/*
	for(i=0;i< pPts->vNum;i++){
	fprintf(filep, "%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f  %.8f  %.8f  %.8f \n",
	pPts->V[i].x,    pPts->V[i].y,    pPts->V[i].z,
	pPts->Dir1[i].x, pPts->Dir1[i].y, pPts->Dir1[i].z,
	pPts->Dir2[i].x, pPts->Dir2[i].y, pPts->Dir2[i].z,
	pPts->K1[i], pPts->K2[i], abs(pPts->K1[i]*pPts->K2[i]));
	}
	fclose(filep);
	printf(" %d Vertices \n", pPts->vNum );
	return 1;
	*/
	for (i = 0; i< pPts->vNum; i++) {
		fprintf(filep, "%.8f %.8f %.8f %.8f \n",
			//pPts->V[i].x, pPts->V[i].y, pPts->V[i].z, abs(pPts->K1[i] * pPts->K2[i]));
			pPts->V[i].x, pPts->V[i].y, pPts->V[i].z, pPts->K1[i] * pPts->K2[i] );
			//pPts->V[i].x, pPts->V[i].y, pPts->V[i].z, (pPts->K1[i] + pPts->K2[i])/2 );

	}
	fclose(filep);
	printf(" %d Vertices \n", pPts->vNum);
	return 1;
}


//-----------------basic algorithms -----------------------
int CreateGlobKdTree(void) //int neibNum
{
	//ANNpointArray dataPts;
	dataPts = annAllocPts(pScanD->vNum, 3);
	for (int i = 0; i< pScanD->vNum; i++) {
		dataPts[i][0] = pScanD->V[i].x;
		dataPts[i][1] = pScanD->V[i].y;
		dataPts[i][2] = pScanD->V[i].z;
	}
	GkdTree = NULL;
	GkdTree = new ANNkd_tree(dataPts, pScanD->vNum, 3);
	if (GkdTree)
		return 1;
	else
		return 0;
}

int FreeGlobKdTree(void) //int neibNum
{
	//ANNpointArray dataPts;
	delete GkdTree;
	GkdTree = NULL;
	return 1;
}



int NormalVectorCalculatingWithQFit(PCDAnaOBJ_t *pPtCdd, int numNeighbor)
//void WeightSVDCalNrms(Point_t *Pts, int numberPts,ANNkd_tree *kdTree,Point_t *nrmPts, int numNeighbor )
{
	//---------------参考czlPMA2006 ----------  

	int i, j;
	int numberPts = pPtCdd->vNum;
	Point_t *nrmPts = pPtCdd->Vn;
	Point_t *Pts = pPtCdd->V;

	Point_t *nrm0 = new Point_t[numberPts];  //暂存初始法向量
	Point_t *PtNeighbor = new Point_t[numNeighbor];//用于存近邻的数据
												   //	PtNeighbor=new Point_t[numNeighbor];
	Point_t *nrmNeighbor = new Point_t[numNeighbor];//暂存近邻点坐标
	int   *qArraySN = new int[numNeighbor];//暂存近邻们的的点在数组的序号

	int rowOfM = numNeighbor; // M 的行数
	float **M = new float*[numNeighbor + 1];  //svdcmp(M,  rowOfM, 3, ww, vv);
	for (i = 0; i<numNeighbor + 1; i++)
		M[i] = new float[4];
	float *ww = new float[4];
	float **vv = new float*[numNeighbor + 1];
	for (i = 0; i<numNeighbor + 1; i++)
		vv[i] = new float[4];

	int k;
	float minW;

	Point_t p;
	float *nrmTp = new float[3];

	float vectorLen = 0, delta_norm = (numNeighbor - 1) / 2 * (4.0e-1), sumWeight = 0;
	float *weight = new float[numNeighbor];

	//--------  step 1 ------------------------------

	ANNpoint      queryPt;
	ANNidxArray   nnIdx;
	ANNdistArray  dists;
	double epsKdtree = 0;
	queryPt = annAllocPt(3);

	nnIdx = new ANNidx[numNeighbor + 1];
	dists = new ANNdist[numNeighbor + 1];


	for (i = 0; i<numberPts; i++) { //for(i=0;i<numberPts;i++)
		nrm0[i].x = 0;  nrm0[i].y = 0;  nrm0[i].z = 0;
		//找出numNeighbor个最近邻,并且计算出PQi,再转化为单位向量
		p = Pts[i];
		nrmTp[0] = 0; nrmTp[1] = 0; nrmTp[2] = 0; sumWeight = 0;
		//		findNearNeighbors( i,Pts, numberPts, numNeighbor, qArraySN, disM );//返回近邻点的序号
		queryPt[0] = Pts[i].x;
		queryPt[1] = Pts[i].y;
		queryPt[2] = Pts[i].z;
		GkdTree->annkSearch(queryPt, numNeighbor + 1, nnIdx, dists, epsKdtree);
		for (j = 0; j<numNeighbor; j++)
			qArraySN[j] = nnIdx[j + 1];//若不含自身，需加找一个近邻点

		for (j = 0; j<numNeighbor; j++) {
			PtNeighbor[j].x = Pts[qArraySN[j]].x - p.x;
			PtNeighbor[j].y = Pts[qArraySN[j]].y - p.y;
			PtNeighbor[j].z = Pts[qArraySN[j]].z - p.z;
			vectorLen = sqrtf(PtNeighbor[j].x*PtNeighbor[j].x + PtNeighbor[j].y*PtNeighbor[j].y + PtNeighbor[j].z*PtNeighbor[j].z);
			weight[j] = vectorLen / delta_norm;
			if (weight[j] >= 1) {
				weight[j] = 0;
			}
			else {
				weight[j] = (1 - weight[j])*(1 - weight[j])*(1 - weight[j])*(1 - weight[j])*(1 + 4 * weight[j]);
				sumWeight = sumWeight + weight[j];
				nrmTp[0] = nrmTp[0] + Pts[qArraySN[j]].x*weight[j];
				nrmTp[1] = nrmTp[1] + Pts[qArraySN[j]].y*weight[j];
				nrmTp[2] = nrmTp[2] + Pts[qArraySN[j]].z*weight[j];
			}
		}
		nrmTp[0] = nrmTp[0] / sumWeight;//mean
		nrmTp[1] = nrmTp[1] / sumWeight;
		nrmTp[2] = nrmTp[2] / sumWeight;

		rowOfM = 0;
		for (j = 0; j<numNeighbor; j++) {
			if (weight[j]>0) {
				rowOfM++;
				M[rowOfM][1] = (Pts[qArraySN[j]].x - nrmTp[0]) * weight[j];
				M[rowOfM][2] = (Pts[qArraySN[j]].y - nrmTp[1]) * weight[j];
				M[rowOfM][3] = (Pts[qArraySN[j]].z - nrmTp[2]) * weight[j];
			}
		}
		for (j = 1; j <= 3; j++) {
			vv[j][1] = 1.0; vv[j][2] = 1.0; vv[j][3] = 1.0;
			ww[j] = 0.0;
		}
		svdcmp(M, rowOfM, 3, ww, vv);
		//选出最小奇异值SVD_W[0],和次小奇异值minW及其对应的特征向量

		if (ww[1] <= ww[2]) {
			ww[0] = ww[1]; minW = ww[2]; k = 1;
		}
		else {
			ww[0] = ww[2]; minW = ww[1]; k = 2;
		}
		if (ww[3] <= ww[0]) {
			minW = ww[0]; ww[0] = ww[3]; k = 3;
		}
		else if (ww[3] <= minW)
			minW = ww[3];

		//		cout<<"最小奇异值："<<SVD_W[0]<<endl;;

		if (minW>1e-16) {
			nrm0[i].x = -vv[1][k];
			nrm0[i].y = -vv[2][k];
			nrm0[i].z = -vv[3][k];
		}
		//		nrmPts[i]=nrm0[i];		
	}

	//--------  step 2 normalvoting ------------------------
	float *weight2 = new float[numNeighbor];
	float *yy = new float[numNeighbor];
	float *alpha = new float[numNeighbor];
	float *opi = new float[numNeighbor];
	Point_t *nrv = new Point_t[numberPts];  //
	int *idx = new int[numNeighbor];
	int idx2;
	float eps = 1.0e-12;

	for (i = 0; i<numberPts; i++) { //for(i=0;i<numberPts;i++)
		nrmPts[i].x = nrm0[i].x;
		nrmPts[i].y = nrm0[i].y;
		nrmPts[i].z = nrm0[i].z;

		//找出numNeighbor个最近邻,并且计算出PQi,再转化为单位向量
		p = Pts[i];
		nrmTp[0] = 0;
		nrmTp[1] = 0;
		nrmTp[2] = 0;
		sumWeight = 0;
		idx2 = 0;

		//		findNearNeighbors( i,Pts, numberPts, numNeighbor, qArraySN, disM );//返回近邻点的序号
		queryPt[0] = Pts[i].x;
		queryPt[1] = Pts[i].y;
		queryPt[2] = Pts[i].z;
		GkdTree->annkSearch(queryPt, numNeighbor + 1, nnIdx, dists, epsKdtree);
		for (j = 0; j<numNeighbor; j++)
			qArraySN[j] = nnIdx[j + 1];//若不含自身，需加找一个近邻点

		for (j = 0; j<numNeighbor; j++) {
			PtNeighbor[j].x = Pts[qArraySN[j]].x - p.x;//dif[][]
			PtNeighbor[j].y = Pts[qArraySN[j]].y - p.y;
			PtNeighbor[j].z = Pts[qArraySN[j]].z - p.z;

			nrmNeighbor[j].x = nrm0[qArraySN[j]].x;//nr[][]
			nrmNeighbor[j].y = nrm0[qArraySN[j]].y;
			nrmNeighbor[j].z = nrm0[qArraySN[j]].z;

			yy[j] = PtNeighbor[j].x*nrmNeighbor[j].x + PtNeighbor[j].y*nrmNeighbor[j].y + PtNeighbor[j].z*nrmNeighbor[j].z;

			vectorLen = PtNeighbor[j].x*PtNeighbor[j].x + PtNeighbor[j].y*PtNeighbor[j].y + PtNeighbor[j].z*PtNeighbor[j].z;
			alpha[j] = asin(fabs(yy[j]) / sqrtf(vectorLen + eps));

			weight[j] = sqrtf(vectorLen) / delta_norm;
			if (weight[j] >= 1)
				weight[j] = 1;
			weight[j] = (1 - weight[j])*(1 - weight[j])*(1 - weight[j])*(1 - weight[j])*(1 + 4 * weight[j]);

			if (alpha[j] > 0 && alpha[j] <= PI / 4) {
				idx[j] = 1;
				nrv[idx2].x = vectorLen / yy[j] / 2 * nrmNeighbor[j].x - PtNeighbor[j].x;
				nrv[idx2].y = vectorLen / yy[j] / 2 * nrmNeighbor[j].y - PtNeighbor[j].y;
				nrv[idx2].z = vectorLen / yy[j] / 2 * nrmNeighbor[j].z - PtNeighbor[j].z;

				vectorLen = sqrtf(nrv[idx2].x*nrv[idx2].x + nrv[idx2].y*nrv[idx2].y + nrv[idx2].z*nrv[idx2].z);
				nrv[idx2].x = nrv[idx2].x / vectorLen;
				nrv[idx2].y = nrv[idx2].y / vectorLen;
				nrv[idx2].z = nrv[idx2].z / vectorLen;

				weight2[idx2] = weight[j];

				idx2++;
			}
			else
				idx[j] = 0;

		}
		for (j = 0; j<numNeighbor; j++) {
			if (alpha[j] < eps) {
				nrv[idx2].x = nrmNeighbor[j].x;
				nrv[idx2].y = nrmNeighbor[j].y;
				nrv[idx2].z = nrmNeighbor[j].z;
				weight2[idx2] = weight[j];
				idx2++;
				if (idx2 == numNeighbor)
					break;
			}
		}

		rowOfM = 0;
		for (j = 0; j<numNeighbor; j++) {
			if (weight2[j]>0) {
				rowOfM++;
				M[rowOfM][1] = nrv[j].x * weight2[j];
				M[rowOfM][2] = nrv[j].y * weight2[j];
				M[rowOfM][3] = nrv[j].z * weight2[j];
			}
		}
		for (j = 1; j <= 3; j++) {
			vv[j][1] = 1.0; vv[j][2] = 1.0; vv[j][3] = 1.0;
			ww[j] = 0.0;
		}
		svdcmp(M, rowOfM, 3, ww, vv);
		//选出最大奇异值wW[0]及其对应的特征向量

		if (ww[1] >= ww[2]) {
			ww[0] = ww[1]; k = 1;
		}
		else {
			ww[0] = ww[2]; k = 2;
		}
		if (ww[3] > ww[0]) {
			ww[0] = ww[3]; k = 3;
		}

		/*
		if( vv[2][k] >0 ){
		nrmPts[i].x= vv[1][k];
		nrmPts[i].y= vv[2][k];
		nrmPts[i].z= vv[3][k];
		}
		else{
		nrmPts[i].x= -vv[1][k];
		nrmPts[i].y= -vv[2][k];
		nrmPts[i].z= -vv[3][k];
		}*/

	}

	delete[]PtNeighbor;
	delete[]nrmNeighbor;//暂存近邻点坐标
	delete[]qArraySN;//暂存近邻们的的点在数组的序号
	delete[]nrmTp;
	//	delete []disM;//距离矩阵
	for (i = 0; i<numNeighbor; i++)
		delete M[i];
	delete[]M;
	delete[]ww;
	for (i = 0; i<numNeighbor + 1; i++)
		delete vv[i];
	delete[]vv;
	delete[]weight;

	delete[]weight2;//此地方溢出
	delete[]yy;
	delete[]alpha;
	delete[]opi;
	delete[]nrv;  //
	delete[]idx;

	delete[]nnIdx;
	delete[] dists;
	return 1;
}


/*-------------------04 curvature ---------------------*/
/*
int curvatureCalculatingWithCubic( PCDAnaOBJ_t *pPtCdd,int numNeighbor   )
//void Adj_normal_cubicMethod(Point_t *Pts, int numOfPt,Point_t *nrmPts, int numNeighbor,
//							double *K1, double *K2,
//							Point_t *PriDirection1,Point_t *PriDirection2 ) //Tog
{//二次曲面拟合,只作列满秩情形
//绝对值小的在K1--P1
Point_t *Pts=pPtCdd->V;
int numOfPt=pPtCdd->vNum;
Point_t *nrmPts=pPtCdd->Vn;
float *K1 = pPtCdd->K1;
float *K2 = pPtCdd->K2;
Point_t *PriDirection1 = pPtCdd->Dir1;
Point_t *PriDirection2 = pPtCdd->Dir2;
float  distTp, adjPtsDist=0;//最近邻的距离,3个最近邻的平均距离,其倒数用于法曲率最大绝对值
//-----------------------------------
int i,j,k;

Point_t *PtNeighbor=new Point_t[numNeighbor+1];//用于存近邻的数据，Point_t[0] 存p，其余[1..num]存近邻
Point_t *nrmNeighbor=new Point_t[numNeighbor];
int     *qArraySN=new int[ numNeighbor ];//暂存近邻们的的点在数组的序号

float **M=new float*[ 3*numNeighbor ];
for(i=0;i<3*numNeighbor;i++)
M[i]=new float[7];
float *b=new float[ 3*numNeighbor ];
float *x=new float[7];//x=[A B C D E F G]
float A,B,C;//and Weingarten matrix [A, B; B, C]
float **aa=new float*[ 3*numNeighbor+1 ];
for(i=0;i<3*numNeighbor+1;i++)
aa[i]=new float[8];
float *ww=new float[8];
float **vv=new float*[ 3*numNeighbor+1 ];
for(i=0;i<3*numNeighbor+1;i++)
vv[i]=new float[8];
float *xT=new float[8];

float *ox=new float[3];
float *oy=new float[3];
float *jubuhua=new float[3];
float theta,phi;

float teXianglian1[2],teXianglian2[2];
float vectorNorm,lamada1,lamada2;
Point_t p,pNrm;

ANNpoint      queryPt;
ANNidxArray   nnIdx;
ANNdistArray  dists;

double epsKdtree = 0 ;
queryPt = annAllocPt(3);

nnIdx = new ANNidx[ numNeighbor+1 ];
dists = new ANNdist[ numNeighbor+1 ];

adjPtsDist = 100000000000000000;
for(i=0;i<numOfPt;i++){
//找出numNeighbor个最近邻
p=Pts[i];
pNrm=nrmPts[i];
queryPt[0] = Pts[i].x;
queryPt[1] = Pts[i].y;
queryPt[2] = Pts[i].z;
GkdTree->annkSearch( queryPt, numNeighbor+1,nnIdx, dists,epsKdtree);
for(j=0;j<numNeighbor;j++)
qArraySN[j] = nnIdx[j+1];//若不含自身，需加找一个近邻点
for(j=0;j<numNeighbor;j++){
PtNeighbor[j].x=Pts[ qArraySN[j] ].x;
PtNeighbor[j].y=Pts[ qArraySN[j] ].y;
PtNeighbor[j].z=Pts[ qArraySN[j] ].z;
nrmNeighbor[j].x=nrmPts[ qArraySN[j] ].x;
nrmNeighbor[j].y=nrmPts[ qArraySN[j] ].y;
nrmNeighbor[j].z=nrmPts[ qArraySN[j] ].z;
}
distTp = (dists[1]+dists[2])/2;
if( adjPtsDist > distTp )
adjPtsDist = distTp ;

//点坐标和法方向都转化为局部坐标系的坐标
theta=acos(nrmPts[i].z);
phi=atan2(nrmPts[i].y,nrmPts[i].x);
ox[0]=-sin(phi);  ox[1]=cos(phi); ox[2]=0;
oy[0]= cos(phi)*cos(theta);
oy[1]= sin(phi)*cos(theta);
oy[2]=-sin(theta);
for(j=0;j<numNeighbor;j++){
jubuhua[0]=PtNeighbor[j].x-p.x;//近邻点
jubuhua[1]=PtNeighbor[j].y-p.y;
jubuhua[2]=PtNeighbor[j].z-p.z;
PtNeighbor[j].x=jubuhua[0]*ox[0]+jubuhua[1]*ox[1]+jubuhua[2]*ox[2];
PtNeighbor[j].y=jubuhua[0]*oy[0]+jubuhua[1]*oy[1]+jubuhua[2]*oy[2];
PtNeighbor[j].z=jubuhua[0]*nrmPts[i].x+jubuhua[1]*nrmPts[i].y+jubuhua[2]*nrmPts[i].z;

jubuhua[0]=nrmNeighbor[j].x;//法方向
jubuhua[1]=nrmNeighbor[j].y;
jubuhua[2]=nrmNeighbor[j].z;
nrmNeighbor[j].x=jubuhua[0]*ox[0]+jubuhua[1]*ox[1]+jubuhua[2]*ox[2];
nrmNeighbor[j].y=jubuhua[0]*oy[0]+jubuhua[1]*oy[1]+jubuhua[2]*oy[2];
nrmNeighbor[j].z=jubuhua[0]*nrmPts[i].x+jubuhua[1]*nrmPts[i].y+jubuhua[2]*nrmPts[i].z;

}
//构造M*x=b
for(j=0;j<numNeighbor;j++){
M[j][0] = PtNeighbor[j].x*PtNeighbor[j].x/2;
M[j][1] = PtNeighbor[j].x*PtNeighbor[j].y;
M[j][2] = PtNeighbor[j].y*PtNeighbor[j].y/2;
M[j][3] = PtNeighbor[j].x*PtNeighbor[j].x*PtNeighbor[j].x;
M[j][4] = PtNeighbor[j].x*PtNeighbor[j].x*PtNeighbor[j].y;
M[j][5] = PtNeighbor[j].x*PtNeighbor[j].y*PtNeighbor[j].y;
M[j][6] = PtNeighbor[j].y*PtNeighbor[j].y*PtNeighbor[j].y;
b[j] = PtNeighbor[j].z;
}
for(j=0;j< numNeighbor;j++){ //i.e.j=numNeighbor;j< 2*numNeighbor;j++
k=j+numNeighbor;
M[k][0] = PtNeighbor[j].x;
M[k][1] = PtNeighbor[j].y;
M[k][2] = 0;
M[k][3] = 3*PtNeighbor[j].x*PtNeighbor[j].x;
M[k][4] = 2*PtNeighbor[j].x*PtNeighbor[j].y;
M[k][5] =   PtNeighbor[j].y*PtNeighbor[j].y;
M[k][6] = 0;
if( fabs(nrmNeighbor[j].z) < 1.0e-12 )
b[k]=0;
else
b[k] = -nrmNeighbor[j].x/nrmNeighbor[j].z;
}
for(j=0;j< numNeighbor;j++){//i.e.j=2*numNeighbor;j< 3*numNeighbor;j++
k=j+2*numNeighbor;
M[k][0] = 0;
M[k][1] = PtNeighbor[j].x;
M[k][2] = PtNeighbor[j].y;
M[k][3] = 0;
M[k][4] =   PtNeighbor[j].x*PtNeighbor[j].x;
M[k][5] = 2*PtNeighbor[j].x*PtNeighbor[j].y;
M[k][6] = 3*PtNeighbor[j].y*PtNeighbor[j].y;
if( fabs(nrmNeighbor[j].z) < 1.0e-12 )
b[k]=0;
else
b[k] = -nrmNeighbor[j].y/nrmNeighbor[j].z;
}

//------SVD 求解M*x=b-----------------------------

for(k=1;k<=3*numNeighbor;k++){
for(j=1;j<=7;j++)	{
aa[k][j]=M[k-1][j-1];
//cout<<aa[k][j]<<" ";
}
//			cout<<b[k-1]<<endl;
}
for(k=1;k<=7;k++)		{
for(j=1;j<=7;j++)
vv[k][j]=1.0;
ww[k]=0.0;
xT[k-1]=0.0;
x[k-1]=0.0;
}

svdcmp(aa,  3*numNeighbor, 7, ww, vv);

//xx=(v*diag(1/wi))*(u'*b)
for(j=1;j<=7;j++)		{
for(k=1;k<=7;k++)			{
if(fabs(ww[j])<1e-12)
vv[k][j]=vv[k][j]*0;
else
vv[k][j]=vv[k][j]/ww[j];
}
}


for(j=1;j<=7;j++)		{
for(k=1;k<=3*numNeighbor;k++)			{
xT[j-1]=xT[j-1]+aa[k][j]*b[k-1];
}

}

//		cout<<"\n\nthe solution is:"<<endl;
for(k=1;k<=7;k++)		{
for(j=1;j<=7;j++)			{
x[k-1]=x[k-1]+vv[k][j]*xT[j-1];
}
//			cout<<x[k-1]<<endl;
}
//---------------------------------------

A=x[0];
B=x[1];
C=x[2];
//写出Weingarten matrix [A, B; B, C]的特征值和特征向量
lamada1=(A+C+sqrtf((A-C)*(A-C)+4*B*B))/2;
lamada2=A+C-lamada1;
if( fabs( lamada1 ) <= fabs( lamada2 ) ){//加个绝对值试一试20151215
K1[i] = lamada1; //写出各个点的主曲率
K2[i] = lamada2;
}
else{
K1[i]=lamada2;
K2[i]=lamada1;
}
if(0){//20160403 limit the max and min of Ki
distTp = fabs( K1[i] );
if( distTp > adjPtsDist   ){
K1[i] = K1[i] / distTp * adjPtsDist;
}
distTp = fabs( K2[i] );
if( distTp > adjPtsDist   ){
K2[i] = K2[i] / distTp * adjPtsDist;
}
}
if( (K1[i]-A)== 0 && B==0 ){
vectorNorm=sqrtf( (K1[i]-C)*(K1[i]-C)+B*B );
teXianglian1[0] = (K1[i]-C)/vectorNorm;               //特征向量1
teXianglian1[1] = B/vectorNorm;
}
else{
vectorNorm=sqrtf( (K1[i]-A)*(K1[i]-A)+B*B );
teXianglian1[0] = B/vectorNorm;       //特征向量1
teXianglian1[1] = (K1[i]-A)/vectorNorm;
}
if( (K2[i]-A)== 0 && B==0 ){
vectorNorm=sqrtf( (K2[i]-C)*(K2[i]-C)+B*B );
teXianglian2[0] = (K2[i]-C)/vectorNorm;              //特征向量2
teXianglian2[1] = B/vectorNorm;
}
else{
vectorNorm=sqrtf( (K2[i]-A)*(K2[i]-A)+B*B );
teXianglian2[0] = B/vectorNorm;      //特征向量2
teXianglian2[1] = (K2[i]-A)/vectorNorm;
}

PriDirection1[i].x=teXianglian1[0]*ox[0]+teXianglian1[1]*oy[0]; //大的主方向
PriDirection1[i].y=teXianglian1[0]*ox[1]+teXianglian1[1]*oy[1];
PriDirection1[i].z=teXianglian1[0]*ox[2]+teXianglian1[1]*oy[2];

PriDirection2[i].x=teXianglian2[0]*ox[0]+teXianglian2[1]*oy[0]; //小的主方向
PriDirection2[i].y=teXianglian2[0]*ox[1]+teXianglian2[1]*oy[1];
PriDirection2[i].z=teXianglian2[0]*ox[2]+teXianglian2[1]*oy[2];

}


/*
cout<<"origi matrix "<<endl;
for(i=0;i<3;i++){
cout<<M[i][0]<<" "<<M[i][1]<<" "<<M[i][2]<<endl;
}
cout<<"detM="<<detM<<endl;
cout<<"inv matrix "<<endl;
for(i=0;i<3;i++){
cout<<invM[i][0]<<" "<<invM[i][1]<<" "<<invM[i][2]<<endl;
}
cout<<"A  B  C= "<<A<<" "<<B<<" "<<C<<endl;
cout<<K1[i]<<" "<<K2[i]<<endl;
cout<<" teXianglian1: "<<teXianglian1[0]<<" "<<teXianglian1[1]<<endl;
cout<<" teXianglian2: "<<teXianglian2[0]<<" "<<teXianglian2[1]<<endl;
*/

/*
for(i=0;i<3;i++)
delete M[i];
delete []M;

delete []b;
delete []ox;
delete []oy;
delete []jubuhua;

delete []qArraySN;

for(i=0;i<3*numNeighbor+1;i++)
delete aa[i];
delete []aa;
delete []ww;
for(i=0;i<3*numNeighbor+1;i++)
delete vv[i];
delete []vv;
delete []xT;

return 1;
}
*/

int curvatureCalculatingWithCAN(PCDAnaOBJ_t *pPtCdd, int numNeighbor)
//void EdgeAngleMethod(Point_t *Pts, int numberPts,Point_t *nrmPts, int numNeighbor,
//					 ANNkd_tree *kdTree,
//					 double *K1, double *K2,
//					 Point_t *PriDirection1,Point_t *PriDirection2 )
{//二次曲面拟合,只作列满秩情形
	Point_t *Pts = pPtCdd->V;
	int     numberPts = pPtCdd->vNum;
	Point_t *nrmPts = pPtCdd->Vn;
	double *K1 = pPtCdd->K1;
	double *K2 = pPtCdd->K2;
	Point_t *PriDirection1 = pPtCdd->Dir1;
	Point_t *PriDirection2 = pPtCdd->Dir2;
	int i, j;
	//	Point_t *PtNeighbor=new Point_t[numNeighbor+1];//用于存近邻的数据，Point_t[0] 存p，其余[1..num]存近邻
	Point_t *PtNeighbor = new Point_t[numNeighbor];
	Point_t *nrmNeighbor = new Point_t[numNeighbor];
	int   *qArraySN = new int[numNeighbor];//暂存近邻们的的点在数组的序号

	float **M = new float*[3];
	for (i = 0; i<3; i++)
		M[i] = new float[3];
	float **invM = new float*[3];
	for (i = 0; i<3; i++)
		invM[i] = new float[3];
	float *b = new float[3];
	float detM = 0;

	float A, B, C;//Weingarten matrix [A, B; B, C]
	float *ox = new float[3];
	float *oy = new float[3];
	float *jubuhua = new float[3];
	float theta, phi, xy, nxy, u, v, sdif, cosb;

	float teXianglian1[2], teXianglian2[2];
	float vectorNorm, lamada1, lamada2;
	Point_t p, qNrm;

	//计算点之间的距离矩阵，便于找到最近邻
	//	float *disM=new float[ numberPts*(numberPts-1)/2 ];//距离矩阵
	//	ptArray2distanceMatrix( Pts, numberPts, disM );//计算任意两点之间的距离

	ANNpoint      queryPt;
	ANNidxArray   nnIdx;
	ANNdistArray  dists;

	double epsKdtree = 0;
	queryPt = annAllocPt(3);

	nnIdx = new ANNidx[numNeighbor + 1];
	dists = new ANNdist[numNeighbor + 1];


	for (i = 0; i<numberPts; i++) { //for(i=0;i<numberPts;i++)
									//找出numNeighbor个最近邻
		p = Pts[i]; qNrm = nrmPts[i];
		//		findNearNeighbors( i,Pts, numberPts, numNeighbor, qArraySN, disM );//返回近邻点的序号
		queryPt[0] = Pts[i].x;
		queryPt[1] = Pts[i].y;
		queryPt[2] = Pts[i].z;
		//GkdTree->annkSearch(queryPt, numNeighbor + 1, nnIdx, dists, epsKdtree);
		GkdTree->annkSearch(queryPt, numNeighbor + 1, nnIdx, dists, epsKdtree);
		for (j = 0; j<numNeighbor; j++)
			qArraySN[j] = nnIdx[j + 1];//若不含自身，需加找一个近邻点

		for (j = 0; j<numNeighbor; j++) {
			PtNeighbor[j].x = Pts[qArraySN[j]].x;
			PtNeighbor[j].y = Pts[qArraySN[j]].y;
			PtNeighbor[j].z = Pts[qArraySN[j]].z;
			//			cout<<qArray[i]<<endl;
			nrmNeighbor[j].x = nrmPts[qArraySN[j]].x;
			nrmNeighbor[j].y = nrmPts[qArraySN[j]].y;
			nrmNeighbor[j].z = nrmPts[qArraySN[j]].z;
		}

		//转化为局部坐标系的坐标
		theta = acos(nrmPts[i].z);
		phi = atan2(nrmPts[i].y, nrmPts[i].x);
		ox[0] = -sin(phi);  ox[1] = cos(phi); ox[2] = 0;
		oy[0] = cos(phi)*cos(theta);
		oy[1] = sin(phi)*cos(theta);
		oy[2] = -sin(theta);
		for (j = 0; j<numNeighbor; j++) {
			jubuhua[0] = PtNeighbor[j].x - p.x;
			jubuhua[1] = PtNeighbor[j].y - p.y;
			jubuhua[2] = PtNeighbor[j].z - p.z;
			PtNeighbor[j].x = jubuhua[0] * ox[0] + jubuhua[1] * ox[1] + jubuhua[2] * ox[2];
			PtNeighbor[j].y = jubuhua[0] * oy[0] + jubuhua[1] * oy[1] + jubuhua[2] * oy[2];
			PtNeighbor[j].z = jubuhua[0] * nrmPts[i].x + jubuhua[1] * nrmPts[i].y + jubuhua[2] * nrmPts[i].z;

			jubuhua[0] = nrmNeighbor[j].x;//法方向
			jubuhua[1] = nrmNeighbor[j].y;
			jubuhua[2] = nrmNeighbor[j].z;
			nrmNeighbor[j].x = jubuhua[0] * ox[0] + jubuhua[1] * ox[1] + jubuhua[2] * ox[2];
			nrmNeighbor[j].y = jubuhua[0] * oy[0] + jubuhua[1] * oy[1] + jubuhua[2] * oy[2];
			nrmNeighbor[j].z = jubuhua[0] * nrmPts[i].x + jubuhua[1] * nrmPts[i].y + jubuhua[2] * nrmPts[i].z;
			vectorNorm = sqrtf(nrmNeighbor[j].x*nrmNeighbor[j].x + nrmNeighbor[j].y*nrmNeighbor[j].y +
				nrmNeighbor[j].z*nrmNeighbor[j].z);
			nrmNeighbor[j].x = nrmNeighbor[j].x / vectorNorm;
			nrmNeighbor[j].y = nrmNeighbor[j].y / vectorNorm;
			nrmNeighbor[j].z = nrmNeighbor[j].z / vectorNorm;
		}
		//写出M=(A'A),  b=A'z 并求出x=[A,B,C]=inv(M) *b
		for (j = 0; j<3; j++) {
			M[j][0] = 0; M[j][1] = 0; M[j][2] = 0; b[j] = 0;
		}
		for (j = 0; j<numNeighbor; j++) {
			xy = sqrt(PtNeighbor[j].x*PtNeighbor[j].x + PtNeighbor[j].y*PtNeighbor[j].y);
			if (fabs(xy) < 1.0e-21)
				xy = 1.e-21;
			u = PtNeighbor[j].x / xy;
			v = PtNeighbor[j].y / xy;
			M[0][0] = M[0][0] + u*u*u*u;
			M[0][1] = M[0][1] + 2 * u*u*u*v;
			M[0][2] = M[0][2] + u*u*v*v;
			M[1][1] = M[1][1] + 4 * u*u*v*v;
			M[1][2] = M[1][2] + 2 * u*v*v*v;
			M[2][2] = M[2][2] + v*v*v*v;

			nxy = nrmNeighbor[j].x*u + nrmNeighbor[j].y*v;
			if (fabs(nxy)<1.0e-21)
				nxy = 1.0e-21;
			if ((PtNeighbor[j].z / xy - nrmNeighbor[j].z / nxy) < -(1.0e-21))
				sdif = -1;
			else if ((PtNeighbor[j].z / xy - nrmNeighbor[j].z / nxy) > (1.0e-21))
				sdif = 1;
			else
				sdif = 0;
			cosb = -sdif*nxy / fabs(nxy)*nrmNeighbor[j].z / sqrt(nrmNeighbor[j].z*nrmNeighbor[j].z + nxy*nxy);
			b[0] = b[0] + u*u*sdif*sqrt(1 - cosb*cosb) / xy;
			b[1] = b[1] + 2 * u*v*sdif*sqrt(1 - cosb*cosb) / xy;
			b[2] = b[2] + v*v*sdif*sqrt(1 - cosb*cosb) / xy;
		}
		M[1][0] = M[0][1]; M[2][0] = M[0][2]; M[2][1] = M[1][2];

		detM = M[0][0] * M[1][1] * M[2][2] + M[1][0] * M[2][1] * M[0][2] + M[2][0] * M[0][1] * M[1][2]
			- M[2][0] * M[1][1] * M[0][2] - M[1][0] * M[0][1] * M[2][2] - M[0][0] * M[2][1] * M[1][2];
		invM[0][0] = (M[1][1] * M[2][2] - M[1][2] * M[2][1]) / detM;
		invM[0][1] = -(M[1][0] * M[2][2] - M[1][2] * M[2][0]) / detM;
		invM[0][2] = (M[1][0] * M[2][1] - M[1][1] * M[2][0]) / detM;
		invM[1][1] = (M[0][0] * M[2][2] - M[0][2] * M[2][0]) / detM;
		invM[1][2] = -(M[0][0] * M[2][1] - M[0][1] * M[2][0]) / detM;
		invM[2][2] = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) / detM;
		invM[1][0] = invM[0][1]; invM[2][0] = invM[0][2]; invM[2][1] = invM[1][2];

		A = invM[0][0] * b[0] + invM[0][1] * b[1] + invM[0][2] * b[2];
		B = invM[1][0] * b[0] + invM[1][1] * b[1] + invM[1][2] * b[2];
		C = invM[2][0] * b[0] + invM[2][1] * b[1] + invM[2][2] * b[2];

		//写出Weingarten matrix [A, B; B, C]的特征值和特征向量
		lamada1 = (A + C + sqrtf((A - C)*(A - C) + 4 * B*B)) / 2;
		lamada2 = A + C - lamada1;
		if (lamada1 >= lamada2) {
			K1[i] = lamada1; //写出各个点的主曲率
			K2[i] = lamada2;
		}
		else {
			K1[i] = lamada2;
			K2[i] = lamada1;
		}
		if (fabs(K1[i] - A) < 1.0e-16 && fabs(B) < 1.0e-16) {
			vectorNorm = sqrtf((K1[i] - C)*(K1[i] - C) + B*B);
			if (vectorNorm < 1.0e-16) {
				teXianglian1[0] = 1;
				teXianglian1[1] = 0;
			}
			else {
				teXianglian1[0] = (K1[i] - C) / vectorNorm;               //特征向量1
				teXianglian1[1] = B / vectorNorm;
			}
		}
		else {
			vectorNorm = sqrtf((K1[i] - A)*(K1[i] - A) + B*B);
			if (vectorNorm < 1.0e-16) {
				teXianglian1[0] = 1;
				teXianglian1[1] = 0;
			}
			else {
				teXianglian1[0] = B / vectorNorm;       //特征向量1
				teXianglian1[1] = (K1[i] - A) / vectorNorm;
			}
		}
		if (fabs(K2[i] - A) <1.0e-16 && fabs(B) <1.0e-16) {
			vectorNorm = sqrtf((K2[i] - C)*(K2[i] - C) + B*B);
			if (vectorNorm < 1.0e-16) {
				teXianglian2[0] = 0;
				teXianglian2[1] = 1;
			}
			else {
				teXianglian2[0] = (K2[i] - C) / vectorNorm;              //特征向量2
				teXianglian2[1] = B / vectorNorm;
			}
		}
		else {
			vectorNorm = sqrtf((K2[i] - A)*(K2[i] - A) + B*B);
			if (vectorNorm < 1.0e-16) {
				teXianglian2[0] = 0;
				teXianglian2[1] = 1;
			}
			else {
				teXianglian2[0] = B / vectorNorm;      //特征向量2
				teXianglian2[1] = (K2[i] - A) / vectorNorm;
			}
		}
		/*		cout<<M[0][0]<<" "<<M[0][1]<<" "<<M[0][2]<<endl;
		cout<<M[1][0]<<" "<<M[1][1]<<" "<<M[1][2]<<endl;
		cout<<M[2][0]<<" "<<M[2][1]<<" "<<M[2][2]<<endl;

		cout<<invM[0][0]<<" "<<invM[0][1]<<" "<<invM[0][2]<<endl;
		cout<<invM[1][0]<<" "<<invM[1][1]<<" "<<invM[1][2]<<endl;
		cout<<invM[2][0]<<" "<<invM[2][1]<<" "<<invM[2][2]<<endl;

		cout<<K1[i]<<" "<<K2[i]<<endl;
		cout<<teXianglian1[0]<<" "<<teXianglian1[1]<<endl;
		cout<<teXianglian2[0]<<" "<<teXianglian2[1]<<endl;
		*/
		PriDirection1[i].x = teXianglian1[0] * ox[0] + teXianglian1[1] * oy[0]; //大的主方向
		PriDirection1[i].y = teXianglian1[0] * ox[1] + teXianglian1[1] * oy[1];
		PriDirection1[i].z = teXianglian1[0] * ox[2] + teXianglian1[1] * oy[2];

		PriDirection2[i].x = teXianglian2[0] * ox[0] + teXianglian2[1] * oy[0];//小的主方向
		PriDirection2[i].y = teXianglian2[0] * ox[1] + teXianglian2[1] * oy[1];
		PriDirection2[i].z = teXianglian2[0] * ox[2] + teXianglian2[1] * oy[2];

	}


	/*
	cout<<"origi matrix "<<endl;
	for(i=0;i<3;i++){
	cout<<M[i][0]<<" "<<M[i][1]<<" "<<M[i][2]<<endl;
	}
	cout<<"detM="<<detM<<endl;
	cout<<"inv matrix "<<endl;
	for(i=0;i<3;i++){
	cout<<invM[i][0]<<" "<<invM[i][1]<<" "<<invM[i][2]<<endl;
	}
	cout<<"A  B  C= "<<A<<" "<<B<<" "<<C<<endl;
	cout<<K1[i]<<" "<<K2[i]<<endl;
	cout<<" teXianglian1: "<<teXianglian1[0]<<" "<<teXianglian1[1]<<endl;
	cout<<" teXianglian2: "<<teXianglian2[0]<<" "<<teXianglian2[1]<<endl;
	*/


	for (i = 0; i<3; i++)
		delete M[i];
	delete[]M;

	for (i = 0; i<3; i++)
		delete invM[i];
	delete[]invM;

	delete[]b;
	delete[]ox;
	delete[]oy;
	delete[]jubuhua;
	//	delete []disM;//距离矩阵
	delete[]qArraySN;

	delete[]nnIdx;
	delete[] dists;

	return 1;
}
