///////////////////////////////////////////////////////////////////
//
// Written by Jianming Zhang (zhangjianm@gmail.com)
// College of Mechanical and Automotive Engineering, Hunan university
// Copyright (c) July 16, 2008
//
///////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "SRC\DLIM\DLIMIntp.h"
#include "SRC\Globals\globals.hpp"
#include "SRC\Globals\algebra.h"
#include "SRC\Globals\formula.h"
#include "SRC\Cells\IntpPtBase.h"
#include "SRC\Structure\Structure.h"
#include "SRC\Lapack\Inc\blas.h"
#include "SRC\Lapack\Inc\mkl_service.h"
#include "SRC\Shapes\LineElmtsShapeFun.h"
#include "Analysis.hpp"
#include "SRC\Analysis\Potential\Potential.hpp"
#include "SRC\Analysis\Elasticity\Elasticity.hpp"
////////////////////////////////////////////////////////////////////
#ifdef _DEBUG 
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif
////////////////////////////////////////////////////////////////////
//[LWC 20161217]
////////////////////////////////////////////////////////////////////
CAnalysis::CAnalysis()
{
	pCmptStructure=NULL;
	pSimulation=NULL;

	tlCmptFce=0;

	constructorInit();
}
CAnalysis::CAnalysis(CCmptStructure *pCmptStruct, CSimulation *pSimu)
{
	pCmptStructure=pCmptStruct;
	pSimulation=pSimu;

	tlCmptFce=pCmptStruct->tlCmptFce;
	vpCmptFce=new CCmptFace *[tlCmptFce];
	for(int i=0; i<tlCmptFce; i++) vpCmptFce[i]=pCmptStruct->vpCmptFce[i];
	vpDmProblem=new CProblem *[tlCmptFce];
	constructorInit();
}
void CAnalysis::constructorInit()
{
	analysisNo_=0;
	analysisKind=potentialAnal;

	tlCrvSrcPt=tlCrvCelPt=0;
	tlSrfSrcPt=tlSrfCelPt=0;

	vDomainNo_=NULL;
	vDomainIndex=NULL;

	aryPrescribedPt=NULL;
	pMatrixAll=NULL;
	vtFAll=NULL;

	logStream=tempStream=NULL;

	bInitialized=FALSE; bAssembled=FALSE; bSolved=FALSE; bResult=FALSE; bSubmitted=FALSE;

	baseFileName="";
	csUnkFile="";
	csResFile="";

	lSizeOfDiscFiles=0.0;
	lSizeOfClustreTree=0.0;//simply summed for all domains
	lSizeOfFullMatrix=0.0;//simply summed for all domains
	lSizeOfRkMatrix=0.0;//simply summed for all domains
	lSizeOfSuperMatrix=0.0;//simply summed for all domains
	lSizeOfSystemFullMatrix=0.0;//simply summed for all domains

	lSizeOfAssembledSuperAll=0.0;
	lSizeOfAssembledSuperAllFull=0.0;
	lSizeOfAssembledSuperAllRank=0.0;

	lSizeOfSolvedSuperAll=0.0;
	lSizeOfSolvedSuperAllFull=0.0;
	lSizeOfSolvedSuperAllRank=0.0;

	lSizeOfDiscFiles=0.0;
}
CAnalysis::~CAnalysis()
{
	finish();
	MKL_FreeBuffers();
}
void CAnalysis::finish()
{
	if(vpDmProblem)
	{
		for(int i=0; i!=tlCmptFce; i++)
		{
			if(vpDmProblem[i]) delete vpDmProblem[i]; vpDmProblem[i]=NULL;
		}
		delete [] vpDmProblem; vpDmProblem=NULL;
	}
	if(vpCmptFce) delete [] vpCmptFce; vpCmptFce=NULL;
	if(vDomainNo_) delete [] vDomainNo_; vDomainNo_=NULL;
	if(vDomainIndex) delete [] vDomainIndex; vDomainIndex=NULL;
	if(tlCmptFce>1)
	{
		if(pMatrixAll) delete [] pMatrixAll; pMatrixAll=NULL;
		if(vtFAll) delete [] vtFAll; vtFAll=NULL;
		if(pMatrixAll_i) delete [] pMatrixAll_i; pMatrixAll_i=NULL;
		if(vtFAll_i) delete [] vtFAll_i; vtFAll_i=NULL;
	}
	if(aryPrescribedPt) delete [] aryPrescribedPt; aryPrescribedPt=NULL;
	
	bInitialized=FALSE; bAssembled=FALSE; bSolved=FALSE; bResult=FALSE; bSubmitted=FALSE;
}
void CAnalysis::finishRecovery()
{
	for(int i=0; i<tlCmptFce; i++)
	{
		vpDmProblem[i]->deleteRecovery();
	}
}
////////////////////////////////////////////////////////////////////
int CAnalysis::initialize()
{
	int i;

	tlCrvSrcPt=pCmptStructure->tlCrvSrcPt;
	tlCrvCelPt=pCmptStructure->tlCrvCelPt;
	tlSrfSrcPt=pCmptStructure->tlSrfSrcPt;
	tlSrfCelPt=pCmptStructure->tlSrfCelPt;
	ASSERT( tlCrvSrcPt != 0);

	vDomainNo_=new int[tlCmptFce];
	vDomainIndex=new int[tlCmptFce+1];

	vDomainIndex[0]=-1;
	for(i=0; i<tlCmptFce; i++)
	{
		vDomainNo_[i]=vpCmptFce[i]->pGeoFce->fceNo_;
		vDomainIndex[vpCmptFce[i]->pGeoFce->fceNo_]=i;
	}

	for(i=0; i<tlCmptFce; i++)
	{
        if(vpDmProblem[i]->initProblem())
		{
			::ErrDisplay("initialize wrong");
			return 1;
		}
	}

	setSingularShpMtx_Vtx0(ctrlDLIM.offset);
	setSingularShpMtx_Vtx1(ctrlDLIM.offset);

	bInitialized=TRUE;

	getSizeInformation();
	outputInfoFile();

	return 0;
}
////////////////////////////////////////////////////////////////////
int CAnalysis::getMatrices()
{
	CTime timeMxStrat;
	Routine_BEGIN("Calculating matrices", &timeMxStrat);
	WriteOneNumberToLog("    The total number of domains", tlCmptFce);

	for(int i=0; i<tlCmptFce; i++)
	{
		Subroutine_BEGIN("Domain No.", i+1);
        if(vpDmProblem[i]->getMatrices())
		{
			::ErrDisplay("Get matrices wrong in CAnalysis");
			return 1;
		}
        if(vpDmProblem[i]->assembler())
		{
			::ErrDisplay("Assembler() fail");
			return 1;
		}
	}

	Routine_END("Calculating matrices", timeMxStrat);

	return 0;
}
int CAnalysis::assembler()//create pSuperAll and get right hand vector vtF
{
	if(assembler_A()) return 1;
	if(assembler_X()) return 1;
	////////////////
	bAssembled=TRUE;
	getSizeInformation();

	return 0;
}
int CAnalysis::assembler_A()
{
	//create pSuperAll
	if(tlCmptFce==1)
	{
		pMatrixAll=vpDmProblem[0]->pDomainMX->pMatrixA;
		pMatrixAll_i=vpDmProblem[0]->pDomainMX->pMatrixA_i;
	}
	else if(tlCmptFce>1)
	{
		//To be completed 
		ASSERT(FALSE);
	}
	else
	{
		ASSERT(FALSE);
		return 1;
	}
	
	return 0;
}
int CAnalysis::assembler_X()//get right hand vector vtF
{
	if(tlCmptFce==1)
	{
		vtFAll=vpDmProblem[0]->pDomainMX->vtF;
		vtFAll_i=vpDmProblem[0]->pDomainMX->vtF_i;
	}
	else if(tlCmptFce>1)
	{
		//To be completed 
		ASSERT(FALSE);
	}
	else
	{
		ASSERT(FALSE);
		return 1;
	}

	return 0;
}
////////////////////////////////////////////////////////////////////
int CAnalysis::solver()
{
	CTime solverBegin;
	Routine_BEGIN("Solver", &solverBegin);
	
	if(InterpolationTest==0)
	{
		switch(analysisKind)
		{
		case potentialAnal:
			{
				if(blas_eqSolver(pMatrixAll, vtFAll, tlCrvSrcPt, 1))
				{
					::ErrDisplay("Solver fail");
					return 1;
				}
			}break;

		case planeStrainAnal:
		case planeStressAnal:
			{
				if(blas_eqSolver(pMatrixAll, vtFAll, tlCrvSrcPt*2, 1))//
				{
					::ErrDisplay("Solver fail");
					return 1;
				}
			}break;

		case acousticsAnal:
			{
				//if(blas_eqSolver_BySVD_c(pMatrixAll_i, vtFAll_i, tlCrvSrcPt, 1))
				if(blas_eqSolver_c(pMatrixAll_i, vtFAll_i, tlCrvSrcPt, 1))
				{
					::ErrDisplay("Solver fail");
					return 1;
				}
			}break;

		default: 
			{
				ASSERT(FALSE);
				::ErrDisplay("Such kind of analysis has not been defined");
				::ErrDisplay("Solver fail");
				return 1;
			}
		}
	}

    Routine_END("  Get right hand vector", solverBegin);

	if(tlCmptFce>1)
	{
		//To be completed 
		ASSERT(FALSE);
	}

	Routine_BEGIN("  Get boundary variables for subdomains", &solverBegin);
	
	if(pCmptStructure->createSearchTree())//add 9.14 �������������ӿ���߽�������һ���ֵ
	{
		::ErrDisplay("create searching tree fail");
		return 1;
	}
    if(vpDmProblem[0]->solver())//������߽��ϵ�ƽ��ֵ
	{
		::ErrDisplay("Solver fail");
		return 1;
	}

	Routine_END("  Get boundary variables for subdomains", solverBegin);
	bSolved=TRUE;
	getSizeInformation();
	//outputUnknownsFile();

	Routine_END("Solver", solverBegin);

	return 0;
}
////////////////////////////////////////////////////////////////////
int CAnalysis::getResults()
{
	for(int k=0; k<tlCmptFce; k++)
	{
        if(vpDmProblem[k]->getResults()) return 1;
	}
	bResult=TRUE;

	return 0;
}
////////////////////////////////////////////////////////////////////
void CAnalysis::outputSource(FILE *stream)
{
	int i;

	if(ctrlOutCore.bOutOfCore==TRUE) fprintf(stream,"The way of memory management is OUT_OF_CORE\n");
	else fprintf(stream,"The way of memory management is IN_CORE\n");
	fprintf(stream,"   The FullPoolNum      =%4d\n",ctrlOutCore.FullPoolNum);
	fprintf(stream,"   The RankPoolNum      =%4d\n",ctrlOutCore.RankPoolNum);
	fprintf(stream,"   The MaxFileSize      =%15d\n",ctrlOutCore.MaxFileSize);
	fprintf(stream,"   The INT_SuperPoolSize=%15d\n",ctrlOutCore.INT_SuperPoolSize);
	fprintf(stream,"   The LUD_SuperPoolSize=%15d\n",ctrlOutCore.LUD_SuperPoolSize);

	fprintf(stream,"\n\nThe Hierarchical matrix control parameters:\n");
	fprintf(stream,"   The maximum rank for Low-rank matrix=%2d.\n",ctrlHMTRX.RKM_k);
	fprintf(stream,"   The max number of nodes in a leaf=%4d\n",ctrlHMTRX.MAX_nodes);
	fprintf(stream,"   The criteria to judge neighbourhood(eta)=%10.4g\n",ctrlHMTRX.RKM_eta);
	fprintf(stream,"   The tolerance for ACA to converge(epson)=%10.4g\n",ctrlHMTRX.ACA_epson);
	fprintf(stream,"   The maximum tolerance for ACA to converge(lowest_epson)=%10.4g\n",ctrlHMTRX.ACA_lowest_epson);

	fprintf(stream,"\nThe GMRES control parameters:\n");
	fprintf(stream,"   The given number of iteration of inner loop=%3d.\n",ctrlGMRES.IT_inner);
	fprintf(stream,"   The given number of iteration of outer loop=%3d.\n",ctrlGMRES.IT_outer);
	fprintf(stream,"   The desired tolerance of convergence=%10.4g\n",ctrlGMRES.tolerance);

	fprintf(stream,"\nThe Accuracy control parameters:\n");
	fprintf(stream,"   FULL_SVD_EPS1=%12.4g\n",svdEPS.Full_1);
	fprintf(stream,"   FULL_SVD_EPS2=%12.4g\n",svdEPS.Full_2);
	fprintf(stream,"   FULL_SVD_EPS3=%12.4g\n\n",svdEPS.Full_3);
	fprintf(stream,"   RK_SVD_EPS1=%12.4g\n",svdEPS.RK_1);
	fprintf(stream,"   RK_SVD_EPS2=%12.4g\n",svdEPS.RK_2);
	fprintf(stream,"   RK_SVD_EPS3=%12.4g\n\n",svdEPS.RK_3);

	fprintf(stream,"The total number of SVDs=%12d\n\n",svdEPS.rk_nrsvd);

/////
	fprintf(stream,"The total domain number=%d.\n\n",tlCmptFce);
	for(i=0; i<tlCmptFce; i++)
	{
		vpDmProblem[i]->outputSource(stream);
	}
/////
	int tlPrescribedPt=0;
	fprintf(stream,"\n");
	fprintf(stream,"The total PRESCRIBED point number=%d.\n\n",tlPrescribedPt);
	if(tlPrescribedPt>0)
	{
		fprintf(stream,"No. Domain No. Face No.            t1             t2      Potential\n");
		for(i=0; i<tlPrescribedPt; i++)
		{
			fprintf(stream,"%2d%10d%10d%15.5g%15.5g%15.5g\n",i+1,aryPrescribedPt[i].zoneIndex,
				aryPrescribedPt[i].faceIndex,aryPrescribedPt[i].pt.t1,aryPrescribedPt[i].pt.t2,aryPrescribedPt[i].u);
		}
		fprintf(stream,"\n");
	}
/////
	fprintf(stream,"\n\n\n");
}
////////////////////////////////////////////////////////////////////
int CAnalysis::import(FILE *stream)
{
	int i, n, tlBemDmn, iFaceNo;
	char str[200];
	CString cs;
	int anaNo_, INTanalysis;

	fscanf(stream,"%s%d", str, &n);
	fscanf(stream,"%s%d", str, &anaNo_);
	fscanf(stream,"%s%d", str, &INTanalysis);
	analysisNo_=anaNo_;
	analysisKind=(enum AnalysisType) INTanalysis;

	fscanf(stream,"%s%d", str, &tlBemDmn);
	ASSERT(tlBemDmn == tlCmptFce);
	cs.Format( _T("%s"), &str );
	cs.MakeUpper();
	if(cs.Compare(_T("TLBEMFACE:")))
	{
		ASSERT(FALSE); return 1;
	}

	/////import domain problems
	int tlOutDM=tlCmptFce-pCmptStructure->tlBondEdg;
	for(i=0; i<tlCmptFce; i++)
	{
		/////import material parameter
		fscanf(stream,"%s%d", str, &iFaceNo);
		ASSERT(iFaceNo == vpCmptFce[i]->pGeoFce->fceNo_);
		cs.Format( _T("%s"), &str );
		cs.MakeUpper();
		if(cs.Compare(_T("FACENO:")))
		{
			ASSERT(FALSE); return 1;
		}

		if(vpCmptFce[i]->material.import(stream))
		{
			vpCmptFce[i]->errorReport("Import material properties fail");
			return 1;
		}

		if(analysisKind == acousticsAnal)
		{
			fscanf(stream,"%s%lf", str, &tlDofForUnitLen);
			ASSERT(tlDofForUnitLen > 0);

			int elmtOrder;
			double waveNumber=vpCmptFce[i]->material.glb_WaveNum;
			double waveLen=TWO_PI/waveNumber;
			for(int j=0; j<vpCmptFce[i]->tlCmptEdg; j++)
			{
				CCurve *pCurve=vpCmptFce[i]->vpCmptEdg[j]->pGeoEdg->pCurve;
				switch(pCurve->NodeSpace.crvIntgCelType)
				{
				case TradIntgCelL1: elmtOrder=1; break;
				case TradIntgCelL2: elmtOrder=2; break;
				case TradIntgCelL3: elmtOrder=3; break;

				case DLIMIntgCelL1: elmtOrder=1; break;
				case DLIMIntgCelL2: elmtOrder=2; break;
				case DLIMIntgCelL3: elmtOrder=3; break;

				default: ASSERT(FALSE);
				}
				double crvLen=pCurve->getCurveLength(0.0, 1.0);
				int nElmt=(int)(tlDofForUnitLen*crvLen/waveLen/elmtOrder+0.5);
				pCurve->NodeSpace.nCrvIntgCel=nElmt;
			}
		}

		//LEAK_CHECK(CProblem, 1)
		//input according to the order of domains and then inclusions
		//if(i < tlOutDM) //BC for out domains
		//{
		//	if(i+1 != k)
		//	{
		//		::ErrDisplay("Wrong boundary condition");
		//		sprintf(str, "\ndomain NO.=%4d", i+1);
		//		::ErrDisplay(str);
		//		finish();
		//		return 1;
		//	}
		//}
		//else //BC for inclusions
		//{
		//	if(i-tlOutDM+1 != k)
		//	{
		//		::ErrDisplay("Wrong boundary condition");
		//		sprintf(str, "\ndomain NO.=%4d", i+1);
		//		::ErrDisplay(str);
		//		finish();
		//		return 1;
		//	}
		//}

		CProblem *pProblem=new CProblem(vpCmptFce[i], this);
		if(pProblem->importInitData(stream))
		{
			::ErrDisplay("Wrong boundary condition");
			sprintf(str, "\ndomain NO.=%4d", i+1);
			::ErrDisplay(str);
			finish();
			return 1;
		}
		vpDmProblem[i]=pProblem;
	}

	fscanf(stream,"%s%d", str, &tlBemDmn);
	ASSERT(tlBemDmn == tlCmptFce);
	cs.Format( _T("%s"), &str );
	cs.MakeUpper();
	if(cs.Compare(_T("TLRESFACE:")))
	{
		ASSERT(FALSE); return 1;
	}
	for(i=0; i<tlCmptFce; i++)
	{
		fscanf(stream,"%s%d", str, &iFaceNo);
		ASSERT(iFaceNo == vpCmptFce[i]->pGeoFce->fceNo_);
		cs.Format( _T("%s"), &str );
		cs.MakeUpper();
		if(cs.Compare(_T("FACENO:")))
		{
			ASSERT(FALSE); return 1;
		}

		if(vpDmProblem[i]->importResData(stream))
		{
			::ErrDisplay("Wrong boundary condition");
			sprintf(str, "\ndomain NO.=%4d", i+1);
			::ErrDisplay(str);
			finish();
			return 1;
		}
	}

	/////  import prescribed boundary conditions at points 
	bInitialized=FALSE; bAssembled=FALSE; bSolved=FALSE; bResult=FALSE; bSubmitted=FALSE;

	baseFileName.Format(_T("%s_Ana%02d"), pSimulation->baseFileName, analysisNo_);
	csUnkFile=baseFileName;
	csUnkFile.Append(_T(".Unk"));
	csResFile=baseFileName;
	csResFile.Append(_T(".res"));

	return 0;
}
int CAnalysis::importRecovery(FILE *stream)
{
	int i;

	finishRecovery();

	///// import locations in each domain
	for(i=0; i<tlCmptFce; i++)
	{
		if(vpDmProblem[i]->importRecovery(stream))
		{
			::ErrDisplay("Import recovery data fail");
			return 1;
		}
	}
	return 0;
}
int CAnalysis::fileReadVTKData(FILE *stream)
{
	int tlPostDmn, iFaceNo;
	char str[200];
	CString cs;
	fscanf(stream,"%s%d", str, &tlPostDmn);
	ASSERT(tlPostDmn == tlCmptFce);
	cs.Format( _T("%s"), &str );
	cs.MakeUpper();
	if(cs.Compare(_T("TLPOSTFACE:")))
	{
		ASSERT(FALSE); return 1;
	}

	for(int i=0; i<tlCmptFce; i++)
	{
		fscanf(stream,"%s%d", str, &iFaceNo);
		ASSERT(iFaceNo == vpCmptFce[i]->pGeoFce->fceNo_);
		cs.Format( _T("%s"), &str );
		cs.MakeUpper();
		if(cs.Compare(_T("FACENO:")))
		{
			ASSERT(FALSE); return 1;
		}

		if(vpDmProblem[i]->fileReadVTKData(stream))
		{
			::ErrDisplay("Wrong VTK Data");
			finish();
			return 1;
		}
	}

	return 0;
}
////////////////////////////////////////////////////////////////////
void CAnalysis::getSizeInformation()
{
	int i;
	lSizeOfClustreTree=0.0;//simply summed for all domains
	lSizeOfFullMatrix=0.0;//simply summed for all domains
	lSizeOfRkMatrix=0.0;//simply summed for all domains
	lSizeOfSuperMatrix=0.0;//simply summed for all domains
	lSizeOfSystemFullMatrix=0.0;//simply summed for all domains

//	lSizeOfSurfaceInvH=0.0;//simply summed for all surfaces

	for(i=0; i<tlCmptFce; i++)
	{
		lSizeOfClustreTree+=vpDmProblem[i]->pDomainMX->lSizeOfClustreTree;
		lSizeOfFullMatrix+=vpDmProblem[i]->pDomainMX->lSizeOfFullMatrix;
		lSizeOfRkMatrix+=vpDmProblem[i]->pDomainMX->lSizeOfRkMatrix;
		lSizeOfSuperMatrix+=vpDmProblem[i]->pDomainMX->lSizeOfSuperMatrix;
		lSizeOfSystemFullMatrix+=vpDmProblem[i]->pDomainMX->lSizeOfSystemFullMatrix;
	}

	/*if(pSuperAll)
	{
        if(bAssembled && !bSolved)
		{
			lSizeOfAssembledSuperAll=pSuperAll->getSize()/1.0e6;
			lSizeOfAssembledSuperAllFull=pSuperAll->getFullMatrixSize()/1.0e6;
			lSizeOfAssembledSuperAllRank=pSuperAll->getRKMatrixSize()/1.0e6;
		}
        if(bSolved)
		{
			lSizeOfSolvedSuperAll=pSuperAll->getSize()/1.0e6;
			lSizeOfSolvedSuperAllFull=pSuperAll->getFullMatrixSize()/1.0e6;
			lSizeOfSolvedSuperAllRank=pSuperAll->getRKMatrixSize()/1.0e6;
		}
	}
	if(ctrlOutCore.bOutOfCore==TRUE)
	{
		lSizeOfDiscFiles=pOutCoreStruct->getSizeOfAllFiles()/1.0e6;
	}*///
}
void CAnalysis::outputInfoFile()
{
	int i;
	//output initial information
	CString csInfoFile;
	csInfoFile=baseFileName;
	csInfoFile.Append(_T(".inf"));
	tempStream = fopen( (char*)(LPCTSTR)csInfoFile, "w");
	if( tempStream != NULL )
	{
		fprintf(tempStream,"\n\n   The given parameters:\n"); 
		fprintf(tempStream,"The node control parameters:\n");
		fprintf(tempStream,"   The dSupportScale=%4.2g, cDScale=%4.2g.\n",
								ctrlNode.dSupportScale,ctrlNode.cDScale);
		fprintf(tempStream,"   t1Step=%d,  t2Step=%d\n",ctrlNode.t1Step,ctrlNode.t2Step);
		fprintf(tempStream,"   IntgCellType=%d,  InptOffset=%4.2g\n",ctrlNode.glbElmtIntpType,ctrlNode.InptOffset);

		fprintf(tempStream,"\nThe Hierarchical matrix control parameters:\n");
		fprintf(tempStream,"   The maximum rank for Low-rank matrix=%2d.\n",ctrlHMTRX.RKM_k);
		fprintf(tempStream,"   The max number of nodes in a leaf=%4d\n",ctrlHMTRX.MAX_nodes);
		fprintf(tempStream,"   The criteria to judge neighbourhood(eta)=%12.4g\n",ctrlHMTRX.RKM_eta);
		fprintf(tempStream,"   The tolerance for ACA to converge(epson)=%12.4g\n",ctrlHMTRX.ACA_epson);
		fprintf(tempStream,"   The maximum tolerance for ACA to converge(lowest_epson)=%12.4g\n",ctrlHMTRX.ACA_lowest_epson);

		fprintf(tempStream,"\nThe GMRES control parameters:\n");
		fprintf(tempStream,"   The given number of iteration of inner loop=%3d.\n",ctrlGMRES.IT_inner);
		fprintf(tempStream,"   The given number of iteration of outer loop=%3d.\n",ctrlGMRES.IT_outer);
		fprintf(tempStream,"   The desired tolerance of convergence=%10.4g\n",ctrlGMRES.tolerance);

		fprintf(tempStream,"\nThe Accuracy control parameters:\n");
		fprintf(tempStream,"   FULL_SVD_EPS1=%12.4g\n",svdEPS.Full_1);
		fprintf(tempStream,"   FULL_SVD_EPS2=%12.4g\n",svdEPS.Full_2);
		fprintf(tempStream,"   FULL_SVD_EPS3=%12.4g\n\n",svdEPS.Full_3);
		fprintf(tempStream,"   RK_SVD_EPS1=%12.4g\n",svdEPS.RK_1);
		fprintf(tempStream,"   RK_SVD_EPS2=%12.4g\n",svdEPS.RK_2);
		fprintf(tempStream,"   RK_SVD_EPS3=%12.4g\n\n",svdEPS.RK_3);

		fprintf(tempStream,"The total number of SVDs=%12d\n\n",svdEPS.rk_nrsvd);

		fprintf(tempStream,"\nTotal unknown=%12d\n", tlCrvSrcPt);

		//////////////////////matrices size information
		fprintf(tempStream,"\n-----SIZE INFORMATION FOR MATRICES------\n");
		fprintf(tempStream,"The virtual full system matrix size for all domains:%10.4g MB\n", lSizeOfSystemFullMatrix);
		fprintf(tempStream,"The summed super matrices size for all domains:  %10.4g MB\n", lSizeOfSuperMatrix);
		fprintf(tempStream,"   The summed full matrix size for all domains:  %10.4g MB\n", lSizeOfFullMatrix);
		fprintf(tempStream,"   The summed rank matrix size for all domains:  %10.4g MB\n", lSizeOfRkMatrix);
		fprintf(tempStream,"The cluster tree size for all domains:  %10.4g MB\n\n", lSizeOfClustreTree);
		fprintf(tempStream,"The inversed H matrix size for all surfaces:  %10.4g MB\n\n", pCmptStructure->lSizeOfSurfaceInvH);

		if(ctrlOutCore.bOutOfCore==TRUE)
		{
			fprintf(tempStream,"Sizes of Files\n");
			fprintf(tempStream,"The size of DISC FILEs for Assembled SuperAll:  %10.4g MB\n\n", lSizeOfDiscFiles);
			for(i=0; i<tlCmptFce; i++)
			{
				fprintf(tempStream,"\nDomain No.:%4d\n", vpCmptFce[i]->pGeoFce->fceNo_);
				fprintf(tempStream,"     Size for DISC FILEs =%10.4g MB\n", vpDmProblem[i]->pDomainMX->lSizeOfDiscFiles);
			}
		}
		////////////
		for(i=0; i<tlCmptFce; i++)
		{
			fprintf(tempStream,"\nDomain No.=%4d; tlSrcPt=%4d; idxCrvSrcPt=%4d; tlCelPt=%4d; idxCrvCelPt=%4d;\n",
				                 vpCmptFce[i]->pGeoFce->fceNo_, 
								 vpCmptFce[i]->tlCrvSrcPt+vpCmptFce[i]->tlSrfSrcPt, 
								 vpCmptFce[i]->idxSrcPt,
								 vpCmptFce[i]->tlCrvCelPt+vpCmptFce[i]->tlSrfCelPt,
								 vpCmptFce[i]->idxCelPt);

			fprintf(tempStream,"\nEdgeNo.   EdgeKind   cruveID   sideID   idxCrvSrcPt   tlIntgCel   tlSrcPt   tlCelPt\n");
			for(int j=0; j<vpCmptFce[i]->tlCmptEdg; j++)
			{
		        if(vpCmptFce[i]->vpCmptEdg[j])
					fprintf(tempStream,"%5d%9d%10d%7d%15d%10d%8d\n",
										vpCmptFce[i]->vpCmptEdg[j]->pGeoEdg->edgNo_,
										vpCmptFce[i]->vpCmptEdg[j]->pGeoEdg->egKind,
										vpCmptFce[i]->vpCmptEdg[j]->pGeoEdg->crvNo_,
										vpCmptFce[i]->vpCmptEdg[j]->pGeoEdg->sideID,
										vpCmptFce[i]->vpCmptEdg[j]->idxCrvCelPt,
										vpCmptFce[i]->vpCmptEdg[j]->tlCrvSrcPt,
										vpCmptFce[i]->vpCmptEdg[j]->tlCrvIntgCel);
			}
			fprintf(tempStream,"\n\nInformation for Hierarchical Matrices at initializing stage\n");
	       /* if(DmProblemArray[i]->pDomainMX->pCluster)
			{
				fprintf(tempStream,"   The total number of leafs =%5d:\n", DmProblemArray[i]->pDomainMX->pCluster->tlAllLeaf);
				fprintf(tempStream,"      The highest level of leafs =%3d,  the lowest level of leafs =%3d\n",
										DmProblemArray[i]->pDomainMX->pCluster->max_level, DmProblemArray[i]->pDomainMX->pCluster->min_level);
				fprintf(tempStream,"      The max number of nodes in a leaf =%4d,  the min number of nodes in a leaf =%4d\n",
										DmProblemArray[i]->pDomainMX->pCluster->max_nodes, DmProblemArray[i]->pDomainMX->pCluster->min_nodes);
			}*///9.4
		}
		fprintf(tempStream,"\nThe weight and volume fractions of the inclusions are: %12.4g%12.4g", pCmptStructure->Wfraction, pCmptStructure->Vfraction);
		fclose(tempStream);
	}
}
////////////////////////////////////////////////////////////////////
void CAnalysis::output()
{
	CString csResFile;
	csResFile=baseFileName;
	if(csResFile != "")
	{
		csResFile.Append(_T(".res"));
		FILE *ResStream;
		ResStream = fopen( (char*)(LPCTSTR)csResFile, "w");
		if( ResStream != NULL )
		{
			fprintf(ResStream,"====================%s==================\n\n", pSimulation->baseFileName);
			fprintf(ResStream,"////////////////////////The source data/////////////////////\n");
			pCmptStructure->output(ResStream);
			outputSource(ResStream);


			fprintf(ResStream,"////////////////////////The result/////////////////////////////\n");
			outputTimeInfo(ResStream);
			outputResult(ResStream);
			fclose(ResStream);
		}
	}
}
void CAnalysis::outputTimeInfo(FILE *stream)
{
	fprintf(stream,"////////////////////////The time information////////////////\n");
	LONGLONG totalSeconds=0;
	csTime = initBegin.Format( "%A, %B %d, %Y" );
	fprintf(stream,"Beginning on %s\n", csTime);
	csTime = solveEnd.Format( "%A, %B %d, %Y" );
	fprintf(stream,"Ending on %s\n\n", csTime);

	csTime = pCmptStructure->meshBegin.Format( "%H:%M:%S, on %B %d."); //***###modification2
	fprintf(stream,"meshing begin at %s\n", csTime);
	csTime = pCmptStructure->meshEnd.Format( "%H:%M:%S, on %B %d.");
	fprintf(stream,"        end at %s\n", csTime);
	CTimeSpan tSpan = pCmptStructure->meshEnd - pCmptStructure->meshBegin;
	nSeconds=tSpan.GetTotalSeconds();
	totalSeconds+=nSeconds;
	fprintf(stream,"        time used: %d seconds.\n\n", nSeconds);

	csTime = initBegin.Format( "%H:%M:%S, on %B %d.");
	fprintf(stream,"Initializing begin at %s\n", csTime);
	csTime = initEnd.Format( "%H:%M:%S, on %B %d.");
	fprintf(stream,"             end at %s\n", csTime);
	tSpan = initEnd - initBegin;
	nSeconds=tSpan.GetTotalSeconds();
	totalSeconds+=nSeconds;
	fprintf(stream,"             time used: %d seconds.\n\n", nSeconds);

	csTime = getmtxBegin.Format( "%H:%M:%S, on %B %d.");
	fprintf(stream,"Getting matrices begin at %s\n", csTime);
	csTime = getmtxEnd.Format( "%H:%M:%S, on %B %d.");
	fprintf(stream,"                 end at %s\n", csTime);
	tSpan = getmtxEnd - getmtxBegin;
	nSeconds=tSpan.GetTotalSeconds();
	totalSeconds+=nSeconds;
	fprintf(stream,"                 time used: %d seconds.\n\n", nSeconds);

	csTime = assmbBegin.Format( "%H:%M:%S, on %B %d.");
	fprintf(stream,"Assambling begin at %s\n", csTime);
	csTime = assmbEnd.Format( "%H:%M:%S, on %B %d.");
	fprintf(stream,"                 end at %s\n", csTime);
	tSpan = assmbEnd - assmbBegin;
	nSeconds=tSpan.GetTotalSeconds();
	totalSeconds+=nSeconds;
	fprintf(stream,"                 time used: %d seconds.\n\n", nSeconds);

	csTime = solveBegin.Format( "%H:%M:%S, on %B %d.");
	fprintf(stream,"Solving the equations begin at %s\n", csTime);
	csTime = solveEnd.Format( "%H:%M:%S, on %B %d.");
	fprintf(stream,"                      end at %s\n", csTime);
	tSpan = solveEnd - solveBegin;
	nSeconds=tSpan.GetTotalSeconds();
	totalSeconds+=nSeconds;
	fprintf(stream,"                      time used: %d seconds.\n\n", nSeconds);
	fprintf(stream,"Total time used:  %d seconds.\n\n\n", totalSeconds);
}
void CAnalysis::outputResult(FILE *stream)
{
	int i;

	fprintf(stream,"\nThe Total Number of Nodes:%d:\n", tlCrvSrcPt);

	//////////////////////matrices size information
	fprintf(stream,"\n-----SIZE INFORMATION FOR MATRICES------\n");
	fprintf(stream,"The virtual full system matrix size for all domains%10.4g MB\n", lSizeOfSystemFullMatrix);
	fprintf(stream,"The summed super matrices size for all domains:  %10.4g MB\n", lSizeOfSuperMatrix);
	fprintf(stream,"   The summed full matrix size for all domains:  %10.4g MB\n", lSizeOfFullMatrix);
	fprintf(stream,"   The summed rank matrix size for all domains:  %10.4g MB\n", lSizeOfRkMatrix);
	fprintf(stream,"The cluster tree size for all domains:  %10.4g MB\n\n", lSizeOfClustreTree);
	fprintf(stream,"The inversed H matrix size for all surfaces:  %10.4g MB\n\n", pCmptStructure->lSizeOfSurfaceInvH);

	fprintf(stream,"The size of Assembled SuperAll:  %10.4g MB\n", lSizeOfAssembledSuperAll);
	fprintf(stream,"The summed full matrix size of Assembled SuperAll:  %10.4g MB\n", lSizeOfAssembledSuperAllFull);
	fprintf(stream,"The summed rank matrix size of Assembled SuperAll:  %10.4g MB\n\n", lSizeOfAssembledSuperAllRank);

	fprintf(stream,"The size of Solved SuperAll:  %10.4g MB\n", lSizeOfSolvedSuperAll);
	fprintf(stream,"The summed full matrix size of Solved SuperAll:  %10.4g MB\n", lSizeOfSolvedSuperAllFull);
	fprintf(stream,"The summed rank matrix size of Solved SuperAll:  %10.4g MB\n\n", lSizeOfSolvedSuperAllRank);

	////////////////
	if(ctrlOutCore.bOutOfCore==TRUE)
	{
		fprintf(stream,"\nSizes of Files\n");
		fprintf(stream,"The size of DISC FILEs for Assembled SuperAll:  %10.4g MB\n\n", lSizeOfDiscFiles);
		for(i=0; i<tlCmptFce; i++)
		{
			fprintf(stream,"\nDomain No.:%4d\n", vpCmptFce[i]->pGeoFce->fceNo_);
			fprintf(stream,"     Size for DISC FILEs =%10.4g MB\n", vpDmProblem[i]->pDomainMX->lSizeOfDiscFiles);
		}
	}

	for(i=0; i<tlCmptFce; i++)
	{
		fprintf(stream,"\n\nDomain No.=%4d\n", vpCmptFce[i]->pGeoFce->fceNo_);
		fprintf(tempStream,"\nDomain No.=%4d; tlSrcPt=%4d; idxCrvSrcPt=%4d; tlCelPt=%4d; idxCrvCelPt=%4d;\n",
                             vpCmptFce[i]->pGeoFce->fceNo_, 
                             vpCmptFce[i]->tlCrvSrcPt+vpCmptFce[i]->tlSrfSrcPt, 
                             vpCmptFce[i]->idxSrcPt,
                             vpCmptFce[i]->tlCrvCelPt+vpCmptFce[i]->tlSrfCelPt,
                             vpCmptFce[i]->idxCelPt);

		fprintf(tempStream,"\nEdgeNo.   EdgeKind   cruveID   sideID   idxCrvSrcPt   tlIntgCel   tlSrcPt   tlCelPt\n");
		for(int j=0; j<vpCmptFce[i]->tlCmptEdg; j++)
		{
			if(vpCmptFce[i]->vpCmptEdg[j])
				fprintf(tempStream,"%5d%9d%10d%7d%15d%10d%8d\n",
									vpCmptFce[i]->vpCmptEdg[j]->pGeoEdg->edgNo_,
									vpCmptFce[i]->vpCmptEdg[j]->pGeoEdg->egKind,
									vpCmptFce[i]->vpCmptEdg[j]->pGeoEdg->crvNo_,
									vpCmptFce[i]->vpCmptEdg[j]->pGeoEdg->sideID,
									vpCmptFce[i]->vpCmptEdg[j]->idxCrvSrcPt,
									vpCmptFce[i]->vpCmptEdg[j]->tlCrvSrcPt,
									vpCmptFce[i]->vpCmptEdg[j]->tlCrvIntgCel);
		}

		fprintf(stream,"\nInformation for clustering at initializing stage\n");
		/*fprintf(stream,"   The total number of leafs =%5d:\n", DmProblemArray[i]->pDomainMX->pCluster->tlAllLeaf);
		fprintf(stream,"      The highest level of leafs =%3d,  the lowest level of leafs =%3d\n",
								DmProblemArray[i]->pDomainMX->pCluster->max_level, DmProblemArray[i]->pDomainMX->pCluster->min_level);
		fprintf(stream,"      The max number of nodes in a leaf =%4d,  the min number of nodes in a leaf =%4d\n",
		*///9.4						DmProblemArray[i]->pDomainMX->pCluster->max_nodes, DmProblemArray[i]->pDomainMX->pCluster->min_nodes);
	}
	fprintf(stream,"\nThe weight and volume fractions of the inclusions are: %12.4g%12.4g",
				                      pCmptStructure->Wfraction, pCmptStructure->Vfraction);
	fprintf(stream,"\n\n");

	for(i=0; i<tlCmptFce; i++)
	{
		vpDmProblem[i]->outputResults(stream);
	}
}
int  CAnalysis::outputUnknownsFile()
{
	int i,j, tlUnknown;

	if(analysisKind == potentialAnal || analysisKind == acousticsAnal)
	{
		tlUnknown=tlCrvSrcPt;
	}
	else if(analysisKind == planeStrainAnal || analysisKind ==planeStressAnal)
	{
		tlUnknown=tlCrvSrcPt*2;
	}
	else
	{
		::ErrDisplay("Such kind of analysis has not been defined yet");
		::ErrDisplay("outputUnknownsFile failed");
		return 1;
	}

	//	CString csUnkFile;
	//	csUnkFile=baseFileName;
	if(csUnkFile != "")
	{
		//		csUnkFile.Append(_T(".Unk"));
		FILE *UnkStream;
		//		UnkStream = fopen( (char*)(LPCTSTR)csUnkFile, "a+");
		UnkStream = fopen( (char*)(LPCTSTR)csUnkFile, "w");
		if( UnkStream != NULL )
		{
			long ni=tlUnknown/10;
			long nj=tlUnknown-ni*10;
			fprintf(UnkStream,"\n***************   vtF   ***************\n");
			for(i=0; i<ni; i++)
			{
				fprintf(UnkStream,"%7d ", i*10);
				for(j=0; j<10; j++) fprintf(UnkStream,"%11.4g, ", vtFAll[i*10+j]);
				fprintf(UnkStream,"\n");
			}
			fprintf(UnkStream,"%7d ", ni*10);
			for(j=0; j<nj; j++) fprintf(UnkStream,"%11.4g, ", vtFAll[ni*10+j]);
			fprintf(UnkStream,"\n");

			fclose(UnkStream);
		}
	}
	return 0;
}
////////////////////////////////////////////////////////////////////
void CAnalysis::outputResultRecovery(FILE *stream)
{
	for(int k=0; k<tlCmptFce; k++)
	{
		vpDmProblem[k]->outputResults(stream);
	}
}
void CAnalysis::Routine_BEGIN(char *routineName, CTime *timeBegin)
{
	logStream = NULL;
	if(pSimulation->csLogFile != "") logStream = fopen( (char*)(LPCTSTR)pSimulation->csLogFile, "a+");
	if( logStream != NULL )
	{
		tBegin = CTime::GetCurrentTime();
		csTime = tBegin.Format( "%H:%M:%S, on %A, %B %d, %Y" );
		fprintf(logStream, "  ---%s", routineName);
		fprintf(logStream, " starts at %s\n", csTime);
		fprintf(logStream,"    --...\n");
		fclose(logStream);
		*timeBegin=tBegin;
	}
}
void CAnalysis::Subroutine_BEGIN(char *subClassName, int i)
{
	logStream = NULL;
	if(pSimulation->csLogFile != "") logStream = fopen( (char*)(LPCTSTR)pSimulation->csLogFile, "a+");
	if( logStream != NULL )
	{
		tBegin = CTime::GetCurrentTime();
		csTime = tBegin.Format( "%H:%M:%S, on %A, %B %d, %Y" );
		fprintf(logStream, "        %s %4d begins ", subClassName, i);
		fprintf(logStream, "at %s\n", csTime);
		fprintf(logStream,"    --...\n");
		fclose(logStream);
	}
}
void CAnalysis::Routine_END(char *routineName, CTime timeBegin)
{
	logStream = NULL;
	if(pSimulation->csLogFile != "") logStream = fopen( (char*)(LPCTSTR)pSimulation->csLogFile, "a+");
	if( logStream != NULL )
	{
		tEnd = CTime::GetCurrentTime();
		csTime = tEnd.Format( "%H:%M:%S, on %A, %B %d, %Y" );
		fprintf(logStream, "  --- %s", routineName);
		fprintf(logStream, " ends at %s\n", csTime);
		tSpan = tEnd - timeBegin;
		nSeconds=tSpan.GetTotalSeconds();
		fprintf(logStream,"    --Seconds used: %d\n\n", nSeconds);
		fclose(logStream);
	}
}
void CAnalysis::WriteOneNumberToLog(char *str, int num)
{
	logStream=NULL;
	if(pSimulation->csLogFile != "") logStream = fopen( (char*)(LPCTSTR)pSimulation->csLogFile, "a+");
	if( logStream != NULL )
	{
		if(num<10000000) fprintf(logStream, "%s is %6d", str, num);
		else
		{
			double dNum=num;
			fprintf(logStream, "%s is %13.7g", str, dNum);
		}
		fprintf(logStream, "\n\n");
		fclose(logStream);
	}
}
void CAnalysis::WriteTwoNumbersToLog(char *str1, char *str2, int num1, int num2)
{
	logStream=NULL;
	if(pSimulation->csLogFile != "") logStream = fopen( (char*)(LPCTSTR)pSimulation->csLogFile, "a+");
	if( logStream != NULL )
	{
		if(num1<10000000) fprintf(logStream, "%s is %6d", str1, num1);
		else
		{
			double dNum=num1;
			fprintf(logStream, "%s is %13.7g", str1, dNum);
		}
		if(num2<10000000) fprintf(logStream, ";   %s is %6d", str2, num2);
		else
		{
			double dNum=num2;
			fprintf(logStream, ";   %s is %13.7g", str2, dNum);
		}
		fprintf(logStream, "\n\n");
		fclose(logStream);
	}
}
void CAnalysis::WriteThreeNumbersToLog(char *str1, char *str2, char *str3, int num1, int num2, int num3)
{
	logStream=NULL;
	if(pSimulation->csLogFile != "") logStream = fopen( (char*)(LPCTSTR)pSimulation->csLogFile, "a+");
	if( logStream != NULL )
	{
		if(num1<10000000) fprintf(logStream, "%s is %6d", str1, num1);
		else
		{
			double dNum=num1;
			fprintf(logStream, "%s is %13.7g", str1, dNum);
		}
		if(num2<10000000) fprintf(logStream, ";   %s is %6d", str2, num2);
		else
		{
			double dNum=num2;
			fprintf(logStream, ";   %s is %13.7g", str2, dNum);
		}
		if(num3<10000000) fprintf(logStream, ";   %s is %6d", str3, num3);
		else
		{
			double dNum=num3;
			fprintf(logStream, ";   %s is %13.7g", str3, dNum);
		}
		fprintf(logStream, "\n\n");
		fclose(logStream);
	}
}
