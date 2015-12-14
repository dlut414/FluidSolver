// PassiveScalarProblemDll.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "PassiveScalarProblemDll.h"

namespace SIM {

	static FSPtr objPtr;

	void PassiveScalarProblemDll::Initialize() {
		objPtr = new FS();
		objPtr->init();
	}

	void PassiveScalarProblemDll::Run() {
		objPtr->stepGL();
	}

	NPtr PassiveScalarProblemDll::Position() {
		return NPtr(objPtr->position());
	}

	NPtr PassiveScalarProblemDll::Scalar() {
		return NPtr(objPtr->scalar());
	}

}