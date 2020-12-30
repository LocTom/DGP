#include "stdafx.h"

#include "colorMeshVerticesCmd.h"
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



colorMeshVerticesCmd::colorMeshVerticesCmd()
{

}

void* colorMeshVerticesCmd::creator()
{
	return new colorMeshVerticesCmd;
}

MString colorMeshVerticesCmd::commandName()
{
	return "colorMeshVerticesCmd";
}

bool colorMeshVerticesCmd::isUndoable() const
{
	return false;
}

void angleUpdater(const MPoint &a, const MPoint &b, const MPoint &c, MDoubleArray &vertexSum, MIntArray &triVertex, int p)
{
	double pos1[4];
	double pos2[4];
	double pos3[4];

	double vector1[3];
	double vector2[3];

	double angle;

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

	MVector vA(vector1);
	MVector vB(vector2);

	angle = vA.angle(vB); //gets angle between two vectors vA and vB
	vertexSum[triVertex[p]] += angle;
}

MStatus	colorMeshVerticesCmd::doIt(const MArgList& argList)
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
		cout << "printing the GMM Matrix: " <<  My_Matrix << endl;
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

	int numVerticesInMeshCreated = 0;
	int numPolygons = meshFn.numPolygons(&stat);

	MItMeshPolygon poly(meshObject);
	if(!poly.isPlanar(&stat) || poly.isLamina(&stat) || poly.isHoled(&stat))
	{
		MCHECKERROR(MS::kFailure, "The given polygon shape is either self intersecting, holed or non-planar which are not supported");
	}

	unsigned int temp; 
	for (int i=0; i<numPolygons; i++)
	{
		temp=poly.polygonVertexCount();
		if ( 3 != temp )
			MCHECKERROR(MS::kFailure, "this is not a triangle mesh!");
		poly.next();
	}


	
	/***************** this part should be changed ****************/
	meshFn.deleteColorSet("NumberColorSet");
	meshFn.deleteColorSet("DegreeColorSet");
	MString s1 = meshFn.createColorSetWithName("NumberColorSet");
	MString s2 = meshFn.createColorSetWithName("DegreeColorSet");
	meshFn.setCurrentColorSetName( s1 );

	MItMeshVertex vertex_it(meshFn.object());
	MIntArray vertexList; 
	MColorArray colors; 

	int curIndex;

	int numEdge; 

	// here I am creating all the colors that the the vertex can be colored with
	MColor red(1, 0, 0);
	MColor blue(0, 0, 1);
	MColor yellow(1, 1, .5);
	MColor green(0, 1, 0);
	MColor magenta(1, 0, 1);
	MColor aqua(0, 1, 1);
	MColor purple(0.5, 0, 1);


	while ( !vertex_it.isDone() ) // this is the loop in which I run through all the nodes and update their color
	{
		curIndex = vertex_it.index(); // gets index of the vertex we are at now
		vertex_it.numConnectedEdges(numEdge); 
		vertexList.append(curIndex);
			if (numEdge <= 3)
			{
				colors.append(purple);
			}
			else if (numEdge == 4)
			{
				colors.append(aqua);
			}
			else if (numEdge == 5)
			{
				colors.append(magenta);
			}
			else if (numEdge == 6)
			{
				colors.append(green);
			}
			else if (numEdge == 7)
			{
				colors.append(yellow);
			}
			else if (numEdge == 8)
			{
				colors.append(blue);
			}
			else
			{
				colors.append(red);
			}
		
		
		vertex_it.next(); //thsi moves to next node using append funciton
	}

	meshFn.setVertexColors ( colors, vertexList);	//takes the array of colors and aray of nodes and applies color to each node
	meshFn.setDisplayColors( true );
	/**************************************************************/
    meshFn.setCurrentColorSetName(s2); // now starting second color set 

	MItMeshPolygon face_it(meshFn.object());
	MDoubleArray vertexSum; 
	MIntArray vertexList2;
	MColorArray colors2; 
	MIntArray triVertex;
	MPointArray triPos;
	MItMeshVertex vertex_itr(meshFn.object());
	MItMeshVertex vertex_itre(meshFn.object());
	float r, g, b;
	double max, min;
	max = -M_PI;
	min = 2 * M_PI;

	int numV;
	numV = meshFn.numVertices();
	

	/*double pos1[4];
	double pos2[4];
	double pos3[4];

	double vector1[3];
	double vector2[3];
	
	double angle; */
	while (numV > 0) // this loop fills the vertex sum array with 0s at the begining
	{
		vertexSum.append(0);
		numV--;
	}

	while (!face_it.isDone()) //this is the loop where I go face by face and find angles for each one of the tree traingles by calling angleUpdater
	{
		face_it.getVertices(triVertex);
		face_it.getPoints(triPos);
		angleUpdater(triPos[0], triPos[1], triPos[2], vertexSum, triVertex,0);
		angleUpdater(triPos[1], triPos[0], triPos[2], vertexSum, triVertex,1);
		angleUpdater(triPos[2], triPos[1], triPos[0], vertexSum, triVertex,2);

		face_it.next();
	}

	while (!vertex_itr.isDone()) //this loop devides nodes into two sections with boundries and without and also updates max and min if need be
	{
		if (vertex_itr.onBoundary())
		{
			vertexSum[vertex_itr.index()] = M_PI - vertexSum[vertex_itr.index()]; //updates list to to pi-sum of angles
			if (vertexSum[vertex_itr.index()] > max) // this if updates max if relavent
			{
				max = vertexSum[vertex_itr.index()];
			}
			if (vertexSum[vertex_itr.index()] < min) //this if updates min if relavent
			{
				min = vertexSum[vertex_itr.index()];
			}
		}

		else 
		{
			vertexSum[vertex_itr.index()] = 2*M_PI - vertexSum[vertex_itr.index()]; //updates list to to 2pi-sum of angles
			if (vertexSum[vertex_itr.index()] > max) // this if updates max if relavent
			{
				max = vertexSum[vertex_itr.index()];
			}
			if (vertexSum[vertex_itr.index()] < min) // this if updates min if relavent
			{
				min = vertexSum[vertex_itr.index()];
			}
		}
	
		vertexList2.append(vertex_itr.index());
		
		vertex_itr.next();
	}

	if (argData.isFlagSet(maxFlag))
	{
		max = argData.flagArgumentDouble(maxFlag, 0);
	}
	if (argData.isFlagSet(minFlag))
	{
		min = argData.flagArgumentDouble(minFlag, 0);
	}

	while (!vertex_itre.isDone())
	{
		mapColor(vertexSum[vertex_itre.index()], r, g, b, min, max);
		colors2.append(MColor(r, g, b)); // sets color to the correct color
		//MGlobal::displayInfo((MString)"The value is " + vertexSum[vertex_itr.index()]);

		vertex_itre.next();
	}


	MGlobal::displayInfo((MString)"The min is " + min + " the max is " + max);
	meshFn.setVertexColors(colors2, vertexList2);
	meshFn.setDisplayColors(true);



	return MS::kSuccess;}

MSyntax colorMeshVerticesCmd::syntax()
{
	MStatus stat = MS::kSuccess;
	MSyntax commandSyntax;

	// Hint - you need to use here the addFlag method of MSyntax class


	stat = commandSyntax.setObjectType(MSyntax::kSelectionList, 1, 1); //expect exactly one object
	MCHECKERRORNORET(stat, "Can't create Syntax object for this command");
	
	// here I am defining flags for min and max and telling the program to be ready to acept doubles
	commandSyntax.addFlag(minFlag, minFlagLong, MSyntax::kDouble);
	commandSyntax.addFlag(maxFlag, maxFlagLong, MSyntax::kDouble);

	commandSyntax.useSelectionAsDefault(true);
	return commandSyntax;}
