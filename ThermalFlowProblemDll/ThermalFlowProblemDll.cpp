/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Published under CC BY-NC
*/
// ThermalFlowProblemDll_.cpp : Defines the exported functions for the DLL application.
/// specify the schemes by defining different solver classes

#include "stdafx.h"
#include "ThermalFlowProblemDll.h"
#include "FractionalStep_KM_A.h"
#include "FractionalStep_DD.h"
#include "FractionalStep_KM_DC.h"
#include "FractionalStep_KM_DC_A.h"
#include "FractionalStep_X.h"
#include <PreInformation.h>

namespace SIM {

	typedef FractionalStep_X<Parameters::DataType, Parameters::Dimension, Parameters::Order> FS;
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
		return NPtr(objPtr->PositionX());
	}
	NPtr ThermalFlowProblemDll2D::PositionY() {
		return NPtr(objPtr->PositionY());
	}
	NPtr ThermalFlowProblemDll2D::VelocityX() {
		return NPtr(objPtr->VelocityX());
	}
	NPtr ThermalFlowProblemDll2D::VelocityY() {
		return NPtr(objPtr->VelocityY());
	}
	NPtr ThermalFlowProblemDll2D::Pressure() {
		return NPtr(objPtr->Pressure());
	}
	NPtr ThermalFlowProblemDll2D::Temperature() {
		return NPtr(objPtr->Temperature());
	}
	NPtr ThermalFlowProblemDll2D::Divergence() {
		return NPtr(objPtr->Divergence());
	}
	void ThermalFlowProblemDll2D::SaveData() {
		objPtr->saveData();
	}
	void ThermalFlowProblemDll2D::SensorOut() {
		objPtr->sensorOut();
	}

}


