#include "stdafx.h"

#include "LSCMCmd.h"
#include <math.h>
#include "Utils/STL_Macros.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/MatlabInterface.h"
#include "Utils/GMM_Macros.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "Utils/Utilities.h"
#include <complex>


#define maxFlag "-max"
#define maxFlagLong "-maxFlag"
#define minFlag "-min"
#define minFlagLong "-minFlagss"


LSCMCmd::LSCMCmd()
{
}


void* LSCMCmd::creator()
{
	return new LSCMCmd;
}

MString LSCMCmd::commandName()
{
	return "LSCMCmd";
}

bool LSCMCmd::isUndoable() const
{
	return false;
}

void cotCalc(const MPoint& a, const MPoint& b, const MPoint& c, double& cotan)
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

	cotan = (va * vb) / ((va ^ vb).length());



}

void updateLaplacianLSCM(const double& alpha, const double& beta, const double& gama, GMMSparseComplexRowMatrix& laplacian, GMMSparseComplexRowMatrix& uvnonboundry, MItMeshVertex& vertex_it2, const int& a, const int& b, const int& c, MIntArray& originalIndex, MFloatArray& u, MFloatArray& v, const int& secondBoundry, const int& firstBoundryVertex)
{

	std::complex<double> im(0,1);





	int previous;
	vertex_it2.setIndex(a, previous);
	if (vertex_it2.index() != secondBoundry && vertex_it2.index() != firstBoundryVertex)
	{
		laplacian(originalIndex[a], originalIndex[a]) -= beta + gama;

		vertex_it2.setIndex(b, previous);
		if (vertex_it2.index() != secondBoundry && vertex_it2.index() != firstBoundryVertex)
		{
			laplacian(originalIndex[a], originalIndex[b]) += (gama - im);
		}
		else
		{
			uvnonboundry(originalIndex[a], 0) -= (gama-im)*((double)u[b] + im * (double)v[b]);
			
		}

		vertex_it2.setIndex(c, previous);
		if (vertex_it2.index() != secondBoundry && vertex_it2.index() != firstBoundryVertex)
		{

			laplacian(originalIndex[a], originalIndex[c]) += beta + im;
		}
		else
		{
			uvnonboundry(originalIndex[a], 0) -= (beta+im)*((double)u[c] + im * (double)v[c]);
		
		}
	}
}



MStatus	LSCMCmd::doIt(const MArgList& argList)
{
	MStatus stat = MS::kSuccess;

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
	int numBoundryVertices = 0;
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
	MFloatArray u(numVerticies, 0);
	MFloatArray v(numVerticies, 0);


	int curIndex;
	int numEdge;
	double outerLength = 0;
	double length;
	int firstVertex = vertex_it.index();
	int firstBoundryVertex;
	int flag = 0;
	int i = 0;
	int previous;
	int j = 0;
	int epreviouse;





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


	// this loop gets the first node on boundry 
	while (!vertex_it.onBoundary())
	{
		vertex_it.next(); //thsi moves to next node using append funciton
	}

	u[vertex_it.index()] = 0;
	v[vertex_it.index()] = 0;

	firstBoundryVertex = vertex_it.index();

	MItMeshEdge edge_itr(meshFn.object());
	MItMeshVertex vertex_itr(meshFn.object());
	MIntArray neighborEdges;
	int numedges;
	double halfLength = outerLength / 2;
	double currentLength = 0;
	flag = 0;
	i = 0;
	j = 0;
	int secondBoundry=0;

	int previous2 = firstBoundryVertex;
	int keep = firstBoundryVertex;
	int trash=0;
	int b;
	numedges = 0;
	int a;

	vertex_itr.setIndex(firstBoundryVertex, trash);

	do
	{
		vertex_itr.getConnectedVertices(vertexNeighborList);

		previous2 = vertex_itr.index();

		while (flag != 1) // this loop goes through preivous nodes neighbors and finds one that is on boundry
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
		vertex_itr.setIndex(previous2, previous);
		vertex_itr.getConnectedEdges(pedgesNeighborList);
		vertex_itr.numConnectedEdges(numedges);
		vertex_itr.setIndex(previous, trash); //brings it back


		while (flag != 1) // this loop goes through edges till its finds their connecting edge
		{
			while (j < numedges && flag != 1) {
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
		currentLength += length;
		// this loop finds the second boundry node
	
			
				if (currentLength > halfLength)
				{
					secondBoundry = edge_it.index(0);
					flag++;
				}
			
				b = vertex_itr.index();


	} while (flag == 0);

	
	u[secondBoundry] = 1;
	v[secondBoundry] = 0;
	


	MItMeshVertex vertex_it2(meshFn.object()); // making a new vertex iterator
	MItMeshPolygon face_it(meshFn.object());
	MIntArray triVertex;
	MPointArray triPos;
	GMMSparseComplexRowMatrix laplacianLSCM(numVerticies - 2, numVerticies - 2);
	GMMSparseComplexRowMatrix uvnonboundryLSCM(numVerticies - 2, 1);
	MIntArray originalIndex(numVerticies, 0);
	MIntArray originalIndexReverse(numVerticies - 2, 0);



	double cotalpha;
	double cotbeta;
	double cotgama;
	int counter = 0;
	int newIndex;

	i, j = 0;

	firstVertex = vertex_it2.index();


	//this loop counts how many of the two boundry nodes it has before it in order to know in which wich row to put it in i
	while (!vertex_it2.isDone())
	{
		if (vertex_it2.index()==secondBoundry || vertex_it2.index() == firstBoundryVertex)
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



		cotCalc(triPos[0], triPos[1], triPos[2], cotalpha);
		cotCalc(triPos[1], triPos[0], triPos[2], cotbeta);
		cotCalc(triPos[2], triPos[1], triPos[0], cotgama);

		//this updates the three rows of the three vertices on face
		updateLaplacianLSCM(cotalpha, cotbeta, cotgama, laplacianLSCM, uvnonboundryLSCM, vertex_it2, triVertex[0], triVertex[1], triVertex[2], originalIndex, u, v, secondBoundry, firstBoundryVertex);
		updateLaplacianLSCM(cotbeta, cotgama, cotalpha, laplacianLSCM, uvnonboundryLSCM, vertex_it2, triVertex[1], triVertex[2], triVertex[0], originalIndex, u, v, secondBoundry, firstBoundryVertex);
		updateLaplacianLSCM(cotgama, cotalpha, cotbeta, laplacianLSCM, uvnonboundryLSCM, vertex_it2, triVertex[2], triVertex[0], triVertex[1], originalIndex, u, v, secondBoundry, firstBoundryVertex);


		face_it.next();
	}
	vertex_it2.setIndex(firstVertex, previous); // bring it back to original vertex


	int result = MatlabGMMDataExchange::SetEngineSparseMatrix("laplacianLSCM", laplacianLSCM); // here is where I send th
	MatlabGMMDataExchange::SetEngineSparseMatrix("uvnonboundryLSCM", uvnonboundryLSCM); // here is where I send th
	MatlabInterface::GetEngine().Eval("SolutionLSCM=solve_linear_system_with_cholesky(-laplacianLSCM,-uvnonboundryLSCM);"); //
	GMMDenseComplexColMatrix SolutionLSCM(numVerticies - 2, 1);
	result = MatlabGMMDataExchange::GetEngineDenseMatrix("SolutionLSCM", SolutionLSCM);
	if(result==-1)
		MCHECKERROR(MStatus::kFailure, "cant get matrix");
	i = 0;
	while (i < numVerticies - 2)
	{
		u[originalIndexReverse[i]] = real(SolutionLSCM(i, 0));
		v[originalIndexReverse[i]] = imag(SolutionLSCM(i, 0));
		int a = originalIndexReverse[i];
		double b = real(SolutionLSCM(i, 0));
		double c = imag(SolutionLSCM(i, 0));
		i++;
	}

	int look= u[secondBoundry];
	v[secondBoundry] = 0;

	meshFn.deleteUVSet("LSCM");
	MString uvSetName = meshFn.createUVSetWithName("LSCM");
	stat = meshFn.setUVs(u, v, &uvSetName);
	MCHECKERROR(stat, "cant set");

	MIntArray uvCounts;
	MIntArray uvIds;
	meshFn.getVertices(uvCounts, uvIds);
	stat=meshFn.assignUVs(uvCounts, uvIds, &uvSetName);
	MCHECKERROR(stat, "cant assign");


	return MS::kSuccess;
}


MSyntax LSCMCmd::syntax()
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
	return commandSyntax;
}
