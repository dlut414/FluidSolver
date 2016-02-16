// ThermalFlowProblemDll.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "ThermalFlowProblemDll.h"
#include "FractionalStep_KM.h"
#include <PreInformation.h>

namespace SIM {

	typedef FractionalStep_KM<Parameters::DataType, Parameters::Dimension, Parameters::Order> FS;
	typedef FS* FSPtr;

	static FSPtr objPtr;

	void ThermalFlowProblemDll::Initialize() {
		objPtr = new FS();
		objPtr->init();
	}

	void ThermalFlowProblemDll::Run() {
		objPtr->stepGL();
	}

	int ThermalFlowProblemDll::Number() {
		return objPtr->part->np;
	}
	NPtr ThermalFlowProblemDll::Type() {
		return NPtr(objPtr->type());
	}
	NPtr ThermalFlowProblemDll::Position() {
		return NPtr(objPtr->position());
	}
	NPtr ThermalFlowProblemDll::Scalar() {
		return NPtr(objPtr->scalar());
	}

	void ThermalFlowProblemDll::SaveData() {
		objPtr->saveData();
	}
	void ThermalFlowProblemDll::SensorOut() {
		objPtr->sensorOut();
	}

}


