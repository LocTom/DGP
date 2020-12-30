#include "stdafx.h"

#include "harmonicFlatteningCmd.h"
#include <math.h>
#include "Utils/STL_Macros.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/MatlabInterface.h"
#include "Utils/GMM_Macros.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "Utils/Utilities.h"


#define maxFlag "-max"
#define maxFlagLong "-maxFlag"
#define minFlag "-min"
#define minFlagLong "-minFlagss"


harmonicFlatteningCmd::harmonicFlatteningCmd()
{
}


void* harmonicFlatteningCmd::creator()
{
	return new harmonicFlatteningCmd;
}

MString harmonicFlatteningCmd::commandName()
{
	return "harmonicFlatteningCmd";
}

bool harmonicFlatteningCmd::isUndoable() const
{
	return false;
}

void angleCalc(const MPoint &a, const MPoint &b, const MPoint &c,double &cotan )
{
	double pos1[4];
	double pos2[4];
	double pos3[4];

	double vector1[3];
	double vector2[3];


	a.get(pos1);
	b.get(pos2);
	c.get(pos3);

	//creats the first vector by subtracing two points from another
	vector1[0] = pos2[0] - pos1[0];
	vector1[1] = pos2[1] - pos1[1];
	vector1[2] = pos2[2] - pos1[2];

	//creats the second vector by subtracing two points from another
	vector2[0] = pos3[0] - pos1[0];
	vector2[1] = pos3[1] - pos1[1];
	vector2[2] = pos3[2] - pos1[2];

	MVector va(vector1);
	MVector vb(vector2);

	cotan = (va * vb) /( (va ^ vb).length());



}

void updateLaplacian(const double &alpha,const double &beta,const double &gama, GMMSparseRowMatrix &laplacian, GMMSparseRowMatrix &uvnonboundry, MItMeshVertex &vertex_it2, const int &a, const int &b, const int &c, MIntArray &originalIndex, MFloatArray &u, MFloatArray &v)
{
	int previous;
	vertex_it2.setIndex(a, previous);
	if (!vertex_it2.onBoundary())
	{
		laplacian(originalIndex[a], originalIndex[a]) -= beta+gama;

		vertex_it2.setIndex(b, previous);
		if (!vertex_it2.onBoundary())
		{
			laplacian(originalIndex[a], originalIndex[b]) += gama;
		}
		else
		{
			uvnonboundry(originalIndex[a], 0) -=(gama) * u[b];
			uvnonboundry(originalIndex[a], 1) -= (gama) * v[b];
		}

		vertex_it2.setIndex(c, previous);
		if (!vertex_it2.onBoundary())
		{
		
			laplacian(originalIndex[a], originalIndex[c]) += beta;
		}
		else
		{
			uvnonboundry(originalIndex[a], 0) -= (beta) * u[c];
			uvnonboundry(originalIndex[a], 1) -= (beta) * v[c];
		}
	}
}

void updateLaplacianUniform(GMMSparseRowMatrix& laplacian, GMMSparseRowMatrix& uvnonboundry, MItMeshVertex& vertex_it2, const int& a, const int& b, const int& c, MIntArray& originalIndex, MFloatArray& u, MFloatArray& v)
{
	int previous;
	vertex_it2.setIndex(a, previous);
	if (!vertex_it2.onBoundary())
	{
		laplacian(originalIndex[a], originalIndex[a]) -= 1;

		vertex_it2.setIndex(b, previous);
		if (!vertex_it2.onBoundary())
		{
			laplacian(originalIndex[a], originalIndex[b]) += 0.5;
		}
		else
		{
			uvnonboundry(originalIndex[a], 0) -= 0.5*u[b];
			uvnonboundry(originalIndex[a], 1) -= 0.5*v[b];
		}

		vertex_it2.setIndex(c, previous);
		if (!vertex_it2.onBoundary())
		{

			laplacian(originalIndex[a], originalIndex[c]) += 0.5;
		}
		else
		{
			uvnonboundry(originalIndex[a], 0) -= 0.5*u[c];
			uvnonboundry(originalIndex[a], 1) -= 0.5*v[c];
		}
	}
}


MStatus	harmonicFlatteningCmd::doIt(const MArgList& argList)
{
	MStatus stat = MS::kSuccess;

	//This code is here just as an example of how to use the Matlab interface.
	//You code for inverting a matrix should be written as part of a new Maya command with the name "inverseMatrixCmd”.
	//test Matlab engine
	if (0)
	{



		MatlabInterface::GetEngine().EvalToCout("My_Matrix = [1 2 3; 4 5 6]"); //creates a 2x3 matrix with name My_Matrix
		GMMDenseColMatrix M(2, 4);
		M(0, 0) = 8.0;
		M(1, 2) = -4.0;
		int result = MatlabGMMDataExchange::SetEngineDenseMatrix("M", M);

		GMMDenseColMatrix My_Matrix;
		result = MatlabGMMDataExchange::GetEngineDenseMatrix("My_Matrix", My_Matrix);
		cout << "printing the GMM Matrix: " << My_Matrix << endl;
	}


	MSyntax commandSyntax = syntax();

	MArgDatabase argData(commandSyntax, argList, &stat);
	MCHECKERROR(stat, "Wrong syntax for command " + commandName());

	MSelectionList objectsList;
	stat = argData.getObjects(objectsList);
	MCHECKERROR(stat, "Can't access object list");

	MObject object;
	stat = objectsList.getDependNode(0, object);
	MCHECKERROR(stat, "Can't access object");

	MObject meshObject;
	stat = Maya_Utils::getMe_a_Mesh(object, meshObject);
	MCHECKERROR(stat, "Object is not a mesh");

	MFnMesh meshFn(meshObject, &stat);
	MCHECKERROR(stat, "Can't access mesh");

	int numPolygons = meshFn.numPolygons(&stat); //number of pollygons/faces in mesh
	int numVerticies = meshFn.numVertices(&stat); //number of vertices in mesh
	int numBoundryVertices=0;
	int numEdges = meshFn.numEdges(&stat); //number of edges in mesh
	int euler = numPolygons + numVerticies - numEdges;


	MItMeshPolygon poly(meshObject);
	if (!poly.isPlanar(&stat) || poly.isLamina(&stat) || poly.isHoled(&stat))
	{
		MCHECKERROR(MS::kFailure, "The given polygon shape is either self intersecting, holed or non-planar which are not supported");
	}

	unsigned int temp;
	for (int i = 0; i < numPolygons; i++)
	{
		temp = poly.polygonVertexCount();
		if (3 != temp)
			MCHECKERROR(MS::kFailure, "this is not a triangle mesh!");
		poly.next();
	}

	//this checks the euler number 
	if (euler != 1)
	{
		MCHECKERROR(MS::kFailure, "The mesh is not a Toplogycal Disk");
	}



	MItMeshVertex vertex_it(meshFn.object());
	MItMeshEdge edge_it(meshFn.object());
	MIntArray vertexNeighborList;
	MIntArray edgesNeighborList;
	MIntArray pedgesNeighborList;
	MFloatArray u(numVerticies,0);
	MFloatArray v(numVerticies,0);


	int curIndex;
	int numEdge; 
	double outerLength=0;
	double length;
	int firstVertex = vertex_it.index();
	int firstBoundryVertex;
	int flag = 0;
	int i = 0;
	int previous;
	double rNumerater = 0;
	double r = 0;
	int j = 0;
	int epreviouse;


	// this loop gets counts the number of boundry nodes
	while (!vertex_it.isDone()) 
	{
		if (vertex_it.onBoundary())
		{
			numBoundryVertices++;
		}
		vertex_it.next(); //thsi moves to next node using append funciton
	}


	// this loop gets me the total length of the boundry
	while (!edge_it.isDone()) 
	{
		if (edge_it.onBoundary())
		{
			edge_it.getLength(length);
			outerLength += length;
		}
		edge_it.next(); //thsi moves to next node using append funciton
	}

	MItMeshVertex vertex_itr(meshFn.object());
	vertex_it.setIndex(firstVertex, previous); //brings iterator back to the first node
	int a = vertex_it.index();

	// this loop gets the first node on boundry 
	while (!vertex_itr.onBoundary()) 
	{
		vertex_itr.next(); //thsi moves to next node using append funciton
		a=vertex_itr.index();
	}
	
	u[vertex_itr.index()] = 1;
	v[vertex_itr.index()] = 0;

	firstBoundryVertex = vertex_itr.index();

	int previous2 = firstBoundryVertex;
	int keep = firstBoundryVertex;
	int numedges;
	int trash;
	int b;
	// this loop sets all the boundry nodes to their locataion in UV
	do
	{
		vertex_itr.getConnectedVertices(vertexNeighborList);

		previous2 = vertex_itr.index();

		while (flag != 1) // this loop goes through preivous nodes nighbors and finds one that is on boundry
		{
			vertex_itr.setIndex(vertexNeighborList[i], trash); //moves vertex to next boundry node
			if (vertex_itr.onBoundary())
			{

				if (keep != vertex_itr.index())
				{
					keep = previous2;
					flag++;
				}
			}
			i++;
		}

		flag = 0;
		i = 0;

		vertex_itr.getConnectedEdges(edgesNeighborList);
		a = vertex_itr.index();
		vertex_itr.setIndex(previous2, previous);
		b = vertex_itr.index();
		vertex_itr.getConnectedEdges(pedgesNeighborList);
		vertex_itr.numConnectedEdges(numedges);
		vertex_itr.setIndex(previous, trash); //brings it back

		while (flag != 1) // this loop goes through edges till its finds their connecting edge
		{
			while (j< numedges && flag != 1) {
				if (edgesNeighborList[i] == pedgesNeighborList[j])
				{
					edge_it.setIndex(edgesNeighborList[i], epreviouse);
					edge_it.getLength(length);
					flag++;
				}
				j++;
			}
			j = 0;
			i++;
		}
		i = 0;
		flag = 0;
		rNumerater += length;
		r = rNumerater / outerLength;
		u[vertex_itr.index()] = cos(2 * M_PI * r);
		v[vertex_itr.index()] = sin(2 * M_PI * r);

	} while (vertex_itr.index() != firstBoundryVertex);


	//here I am copying u and v to u2 and v2 so that I have it for the unifrom
	MFloatArray u2(u);
	MFloatArray v2(v);


	MItMeshVertex vertex_it2(meshFn.object()); // making a new vertex iterator
	MItMeshPolygon face_it(meshFn.object());
	MIntArray triVertex;
	MPointArray triPos;
	GMMSparseRowMatrix laplacian( numVerticies - numBoundryVertices, numVerticies - numBoundryVertices);
	GMMSparseRowMatrix uvnonboundry(numVerticies - numBoundryVertices, 2);
	MIntArray originalIndex(numVerticies, 0);
	MIntArray originalIndexReverse(numVerticies-numBoundryVertices, 0);




	double cotalpha;
	double cotbeta;
	double cotgama;
	int counter=0;
	int newIndex;

	i, j = 0;

	firstVertex = vertex_it2.index();


	//this loop counts how many boundry nodes it has before it in order to know in which wich row to put it in i
	while (!vertex_it2.isDone()) 
	{
		if (vertex_it2.onBoundary())
		{
			counter++;
		}
		else
		{
			originalIndex[vertex_it2.index()] = vertex_it2.index() - counter;
			originalIndexReverse[vertex_it2.index() - counter] = vertex_it2.index();
		}
		vertex_it2.next();
	}

	vertex_it2.setIndex(firstVertex, previous); // bring it back to original vertex


	while (!face_it.isDone()) //this loop 
	{
		face_it.getVertices(triVertex);
		face_it.getPoints(triPos);



		angleCalc(triPos[0], triPos[1], triPos[2], cotalpha);
		angleCalc(triPos[1], triPos[0], triPos[2], cotbeta);
		angleCalc(triPos[2], triPos[1], triPos[0], cotgama);

		//this updates the three rows of the three vertices on face
		updateLaplacian(cotalpha, cotbeta, cotgama, laplacian, uvnonboundry, vertex_it2, triVertex[0], triVertex[1], triVertex[2], originalIndex, u, v);
		updateLaplacian(cotbeta, cotgama, cotalpha, laplacian, uvnonboundry, vertex_it2, triVertex[1], triVertex[2], triVertex[0], originalIndex, u, v);
		updateLaplacian(cotgama ,cotalpha, cotbeta, laplacian, uvnonboundry, vertex_it2, triVertex[2], triVertex[0], triVertex[1], originalIndex, u, v);


		face_it.next();
	}
	vertex_it2.setIndex(firstVertex, previous); // bring it back to original vertex


	int result = MatlabGMMDataExchange::SetEngineSparseMatrix("laplacian", laplacian); // here is where I send th
	MatlabGMMDataExchange::SetEngineSparseMatrix("uvnonboundry", uvnonboundry); // here is where I send th
	MatlabInterface::GetEngine().Eval("Solution=solve_linear_system_with_cholesky(-laplacian,-uvnonboundry);"); //creates a 2x3 matrix with name My_Matrix and inverses it
	GMMDenseColMatrix Solution(numVerticies - numBoundryVertices, 2);
	result = MatlabGMMDataExchange::GetEngineDenseMatrix("Solution", Solution);

	while (i < numVerticies - numBoundryVertices)
	{
			u[originalIndexReverse[i]] = Solution(i,0);
			v[originalIndexReverse[i]] = Solution(i,1);
			double b = Solution(i, 0);
			double a = Solution(i, 1);
			i++;
	}

	meshFn.deleteUVSet("Harmonic-Cotangent");
	MString uvSetName = meshFn.createUVSetWithName("Harmonic-Cotangent");
	meshFn.setUVs(u, v, &uvSetName);
	MIntArray uvCounts;
	MIntArray uvIds;
	meshFn.getVertices(uvCounts, uvIds);
	meshFn.assignUVs(uvCounts, uvIds, &uvSetName);







	//Now I amm adding code for Unifrom set

	GMMSparseRowMatrix laplacianUniform(numVerticies - numBoundryVertices, numVerticies - numBoundryVertices);
	GMMSparseRowMatrix uvnonboundryUniform(numVerticies - numBoundryVertices, 2);
	MItMeshPolygon face_it2(meshFn.object());

	MItMeshVertex vertex_it3(meshFn.object());

	while (!face_it2.isDone()) //this loop 
	{
		face_it2.getVertices(triVertex);
		face_it2.getPoints(triPos);



		//this updates the three rows of the three vertices on face
		updateLaplacianUniform(laplacianUniform, uvnonboundryUniform, vertex_it3, triVertex[0], triVertex[1], triVertex[2], originalIndex, u2, v2);
		updateLaplacianUniform(laplacianUniform, uvnonboundryUniform, vertex_it3, triVertex[1], triVertex[2], triVertex[0], originalIndex, u2, v2);
		updateLaplacianUniform(laplacianUniform, uvnonboundryUniform, vertex_it3, triVertex[2], triVertex[0], triVertex[1], originalIndex, u2, v2);


		face_it2.next();
	}
	vertex_it3.setIndex(firstVertex, previous); // bring it back to original vertex


	int result2 = MatlabGMMDataExchange::SetEngineSparseMatrix("laplacianUniform", laplacianUniform); // here is where I send th
	MatlabGMMDataExchange::SetEngineSparseMatrix("uvnonboundryUniform", uvnonboundryUniform); // here is where I send th
	MatlabInterface::GetEngine().Eval("Solution2=solve_linear_system_with_cholesky(-laplacianUniform,-uvnonboundryUniform);"); //creates a 2x3 matrix with name My_Matrix and inverses it
	GMMDenseColMatrix Solution2(numVerticies - numBoundryVertices, 2);
	result2 = MatlabGMMDataExchange::GetEngineDenseMatrix("Solution2", Solution2);
	i = 0;
	while (i < numVerticies - numBoundryVertices)
	{
		u2[originalIndexReverse[i]] = Solution2(i, 0);
		v2[originalIndexReverse[i]] = Solution2(i, 1);
		i++;
	}

	meshFn.deleteUVSet("Harmonic-Uniform");
	MString uvSetName2 = meshFn.createUVSetWithName("Harmonic-Uniform");
	meshFn.setUVs(u2, v2, &uvSetName2);
	MIntArray uvCounts2;
	MIntArray uvIds2;
	meshFn.getVertices(uvCounts2, uvIds2);
	meshFn.assignUVs(uvCounts2, uvIds2, &uvSetName2);



	return MS::kSuccess;}


MSyntax harmonicFlatteningCmd::syntax()
{
	MStatus stat = MS::kSuccess;
	MSyntax commandSyntax;

	// Hint - you need to use here the addFlag method of MSyntax class


	stat = commandSyntax.setObjectType(MSyntax::kSelectionList, 1, 1); //expect exactly one object
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");
	
	// here I am defining flags for min and max and telling the program to be ready to acept doubles
	//commandSyntax.addFlag(minFlag, minFlagLong, MSyntax::kDouble);
	//commandSyntax.addFlag(maxFlag, maxFlagLong, MSyntax::kDouble);

	commandSyntax.useSelectionAsDefault(true);
	return commandSyntax;}
