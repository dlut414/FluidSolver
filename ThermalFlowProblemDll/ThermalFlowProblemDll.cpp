/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Published under CC BY-NC
*/
// ThermalFlowProblemDll_.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "ThermalFlowProblemDll.h"
#include "FractionalStep_KM.h"
#include <PreInformation.h>

namespace SIM {

	typedef FractionalStep_KM<Parameters::DataType, Parameters::Dimension, Parameters::Order> FS;
	typedef FS* FSPtr;

	static FSPtr objPtr;

	void ThermalFlowProblemDll2D::Initialize() {
		objPtr = new FS();
		objPtr->init();
	}

	void ThermalFlowProblemDll2D::Run() {
		objPtr->stepGL();
	}

	int ThermalFlowProblemDll2D::Number() {
		return objPtr->part->np;
	}
	NPtr ThermalFlowProblemDll2D::Type() {
		return NPtr(objPtr->type());
	}
	NPtr ThermalFlowProblemDll2D::PositionX() {
		return NPtr(objPtr->positionX());
	}
	NPtr ThermalFlowProblemDll2D::PositionY() {
		return NPtr(objPtr->positionY());
	}
	NPtr ThermalFlowProblemDll2D::Scalar() {
		return NPtr(objPtr->scalar());
	}

	void ThermalFlowProblemDll2D::SaveData() {
		objPtr->saveData();
	}
	void ThermalFlowProblemDll2D::SensorOut() {
		objPtr->sensorOut();
	}

}


