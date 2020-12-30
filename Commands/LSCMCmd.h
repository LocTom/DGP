#pragma once



class LSCMCmd : public MPxCommand
{

public:

	LSCMCmd();
	virtual MStatus	doIt(const MArgList& argList);
	static void* creator();
	static MSyntax syntax();
	static MString commandName();
	virtual bool isUndoable() const;

};
