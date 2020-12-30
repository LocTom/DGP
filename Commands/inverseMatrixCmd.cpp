#include "stdafx.h"

#include "inverseMatrixCmd.h"
#include <math.h>
#include "Utils/STL_Macros.h"
#include "Utils/Maya_Macros.h"
#include "Utils/Maya_Utils.h"
#include "Utils/MatlabInterface.h"
#include "Utils/GMM_Macros.h"
#include "Utils/MatlabGMMDataExchange.h"
#include "Utils/Utilities.h"


inverseMatrixCmd::inverseMatrixCmd()
{

}

void* inverseMatrixCmd::creator()
{
	return new inverseMatrixCmd;
}

MString inverseMatrixCmd::commandName()
{
	return "inverseMatrixCmd";
}

bool inverseMatrixCmd::isUndoable() const
{
	return false;
}

MStatus	inverseMatrixCmd::doIt(const MArgList& argList)
{
	MStatus stat = MS::kSuccess;

	MSyntax commandSyntax = syntax();

	MArgDatabase argData(commandSyntax, argList, &stat);
	MCHECKERROR(stat, "Wrong syntax for command " + commandName());

	GMMDenseColMatrix M(3, 3); // defines 3 by 3 matrix
	int x = 0;
	int y = 0;
	int counter = 0;
	double var = 0;

	// this neseted loop gets arguments and makes them into a 3 by 3 matirx
	while (x < 3)
	{
		while (y < 3) 
		{
			argData.getCommandArgument(counter, var);
			M(x, y) = var;
			y++;
			counter++;
		}
		x++;
		y = 0;
	}

	int result = MatlabGMMDataExchange::SetEngineDenseMatrix("M", M); // here is where I send the 3by3 matrix to matlab
	MatlabInterface::GetEngine().EvalToCout("My_Matrix=inv(M)"); //creates a 2x3 matrix with name My_Matrix and inverses it

	GMMDenseColMatrix My_Matrix;
	 result = MatlabGMMDataExchange::GetEngineDenseMatrix("My_Matrix", My_Matrix); 


	cout << "printing the GMM Matrix: " << My_Matrix << endl; //prints the matrix in the output
	





	return MS::kSuccess;
}

MSyntax inverseMatrixCmd::syntax()
{
	MStatus stat = MS::kSuccess;
	MSyntax commandSyntax;

	// Hint - you need to use here the addFlag method of MSyntax class
	int counter = 0;
	while (counter < 9) {
		commandSyntax.addArg(MSyntax::kDouble);
		counter++;
	}
	

	commandSyntax.useSelectionAsDefault(true);
	return commandSyntax;
}
