#include "stdafx.h"

#include "SpaceDeformer2D.h"
#include "Utils/STL_Macros.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "Utils/MatlabInterface.h"


const MTypeId SpaceDeformer2D::mTypeId(0x6723c);
const MString SpaceDeformer2D::mTypeName("SpaceDeformer2D");

MObject SpaceDeformer2D::mCageAttr;
MObject SpaceDeformer2D::mCoordinateTypeAttr;



SpaceDeformer2D::SpaceDeformer2D() : mIsFirstTime(true)
{

}

SpaceDeformer2D::~SpaceDeformer2D()
{

}

void* SpaceDeformer2D::creator()
{
	return new SpaceDeformer2D();
}

MStatus SpaceDeformer2D::initialize()
{
	MStatus stat;

	MFnTypedAttribute cageAttr;
	mCageAttr = cageAttr.create("cage" ,"cage", MFnData::kMesh, MObject::kNullObj, &stat);
	CHECK_MSTATUS(addAttribute(mCageAttr));
	CHECK_MSTATUS(attributeAffects(mCageAttr, outputGeom));

	MFnEnumAttribute coordinateTypeAttr;
	mCoordinateTypeAttr = coordinateTypeAttr.create("coordinateType" ,"coordinateType", 0, &stat);
	CHECK_MSTATUS(coordinateTypeAttr.setKeyable(true));
	CHECK_MSTATUS(addAttribute(mCoordinateTypeAttr));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Cauchy", 0));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Cauchy Interpolation", 1));
	CHECK_MSTATUS(coordinateTypeAttr.addField("Point to point", 2));
	CHECK_MSTATUS(attributeAffects(mCoordinateTypeAttr, outputGeom));

	return MStatus::kSuccess;
}


MStatus SpaceDeformer2D::deform(MDataBlock& block, MItGeometry& iter, const MMatrix& mat, unsigned int multiIndex)
{
	MStatus stat;

	MDataHandle coordHandle = block.inputValue(mCoordinateTypeAttr, &stat);
	short coordinateType = coordHandle.asShort();

	MDataHandle handle = block.inputValue(mCageAttr, &stat);
	CHECK_MSTATUS_AND_RETURN_IT(stat);
	MObject cageMesh = handle.asMesh();

	MFnMesh cageMeshFn(cageMesh, &stat);
	if(stat != MS::kSuccess)
	{
		return stat;
	}

	updateCage(cageMeshFn);

	if(mIsFirstTime)
	{
		stat = doSetup(iter);
		CHECK_MSTATUS_AND_RETURN_IT(stat);
		mIsFirstTime = false;
	}

	//compute the deformation of all the internal points. This is done by simply multiplying the coordinate matrix by the cage vertices vector
	//gmm::mult(mCoordinates, mCageVertices, mInternalPoints);

	if (coordinateType == 0)
	{

		gmm::mult(mCoordinates, mCageVertices, mInternalPoints);


		//update the new deformed position of all the internal vertices
		for (iter.reset(); !iter.isDone(); iter.next())
		{
			int i = iter.index();

			MPoint pt = iter.position();

			//Complex c(pt[0], pt[1]);

			iter.setPosition(MPoint(mInternalPoints[iter.index()].real(), mInternalPoints[iter.index()].imag(), 0.0));
		}
	}
	
	else if (coordinateType == 1)
	{
		gmm::mult(mOuterCageCoordinates, mCageVertices, mCageVertices);

		gmm::mult(mCoordinates, mCageVertices, mInternalPoints);

		//update the new deformed position of all the internal vertices
		for (iter.reset(); !iter.isDone(); iter.next())
		{
			int i = iter.index();

			MPoint pt = iter.position();

			//Complex c(pt[0], pt[1]);

			iter.setPosition(MPoint(mInternalPoints[iter.index()].real(), mInternalPoints[iter.index()].imag(), 0.0));
		}
	}
	return stat;
}


MStatus SpaceDeformer2D::updateCage(MFnMesh& cageMeshFn)
{
	MStatus stat;

	int numFaces = cageMeshFn.numPolygons(&stat);

	assert(numFaces == 1);

	MIntArray vertexIndices;
	cageMeshFn.getPolygonVertices(0, vertexIndices);
	int numV = vertexIndices.length();
	assert(numV >= 3);

	gmm::clear(mCageVertices);
	gmm::resize(mCageVertices, numV, 1);

	MPointArray vertexArray;
	stat = cageMeshFn.getPoints(vertexArray);

	assert(numV == vertexArray.length());

	for(int i = 0; i < numV; i++)
	{
		MPoint p = vertexArray[i];
		Complex c(p[0], p[1]);
		mCageVertices(i, 0) = c;
	}
	return MS::kSuccess;
}


MStatus SpaceDeformer2D::doSetup(MItGeometry& iter)
{
	MStatus stat;

	MGlobal::displayInfo((MString)"in doSetup ");

	int m = iter.count(&stat); //num internal points
	int n = mCageVertices.nrows(); //num vertices in the cage

	GMMDenseComplexColMatrix orgCage(n,1);
	int i = 0;
	while (i < n)
	{
		orgCage[i]= mCageVertices[i]; // save origanl one
		i++;
	}

	
	 i = 0;
	int prev = 0;
	int next = 0;
	std::complex<double> topVec;
	std::complex<double> botVec;
	std::complex<double> ans;
	std::complex<double> img(0, 1);
	std::complex<double> real(1, 0);

		gmm::clear(mCoordinates);
	gmm::resize(mCoordinates, m, n);

	gmm::clear(mOuterCageCoordinates);
	gmm::resize(mOuterCageCoordinates, n, n);

	gmm::clear(mInternalPoints);
	gmm::resize(mInternalPoints, m, 1);

	while (i < n)
	{
		prev = i - 1;
		next = i + 1;

		if (i == 0)
		{
			prev = n - 1;
		}
		else if (i == n - 1)
		{
			next = 0;
		}

		topVec = (-img) * ((orgCage[next]-orgCage[i]) / abs(orgCage[next] - orgCage[i]));
		botVec = (-img) * ((orgCage[i] - orgCage[prev]) / abs(orgCage[i] - orgCage[prev]));

		ans = (topVec + botVec)/abs(topVec + botVec);

		mCageVertices[i] += 0.01 * (ans);

		i++;
	}
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageAfter", mCageVertices);
	MatlabGMMDataExchange::SetEngineDenseMatrix("cageBefore", orgCage);



	
	int zjp, zjm;
	std::complex<double> ajp;
	std::complex<double> aj;
	std::complex<double> bjp;
	std::complex<double> bjm;
	std::complex<double> bj;


	for ( i = 0; i < n; i++)
	{
		//getting point on cage and putting in z0
		Complex z0 = orgCage[i];

		for (int j = 0; j < n; j++)
		{
			Complex K(0.0, 0.0);

			//This is to update the nodes 
			zjm = j - 1;
			zjp = j + 1;

			if (j == 0)
			{
				zjm = n - 1;
			}
			else if (j == n - 1)
			{
				zjp = 0;
			}

			//Here I update all the vectors so can plug into formula (p stand for + and m for -)
			ajp = (mCageVertices[zjp] - mCageVertices[j]);
			aj = (mCageVertices[j] - mCageVertices[zjm]);

			bjp = (mCageVertices[zjp] - z0);
			bjm = (mCageVertices[zjm] - z0);
			bj = (mCageVertices[j] - z0);

			//Now will plug into formula and update k
			K = (real / (2 * M_PI * (img))) * (((bjp / ajp) * log(bjp / bj)) - ((bjm / aj) * log(bj / bjm)));

			mOuterCageCoordinates(i, j) = K;

						
		}
	}

	int result = MatlabGMMDataExchange::SetEngineDenseMatrix("mOuter", mOuterCageCoordinates); // here is where I send the matrix to matlab
	MatlabInterface::GetEngine().Eval("My_Matrix=inv(mOuter)"); //creates a matrix with name My_Matrix and inverses it
	result = MatlabGMMDataExchange::GetEngineDenseMatrix("My_Matrix", mOuterCageCoordinates);


	for(iter.reset(); !iter.isDone(); iter.next())
	{
		int i = iter.index();
		MPoint pt = iter.position();

		Complex z0(pt[0], pt[1]); //internal point

		mInternalPoints(i, 0) = z0;

		for(int j = 0; j < n; j++)
		{
			Complex K(0.0, 0.0);

			//This is to update the nodes 
			zjm = j - 1;
			zjp = j + 1;

			if (j == 0)
			{
				zjm = n - 1;
			}
			else if (j == n - 1)
			{
				zjp = 0;
			}

			//Here I update all the vectors so can plug into formula (p stand for + and m for -)
			ajp = (mCageVertices[zjp] - mCageVertices[j]);
			aj = (mCageVertices[j] - mCageVertices[zjm]);

			bjp = (mCageVertices[zjp] - z0);
			bjm = (mCageVertices[zjm] - z0);
			bj = (mCageVertices[j] - z0);

			//Now will plug into formula and update k
			K = (real / (2 * M_PI * (img))) * (((bjp / ajp) * log(bjp / bj)) - ((bjm / aj) * log(bj / bjm)));
		
			mCoordinates(i, j) = K;
		}
	}

	return MS::kSuccess;
}
