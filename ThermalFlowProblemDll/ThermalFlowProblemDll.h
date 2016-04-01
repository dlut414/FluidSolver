/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Published under CC BY-NC
*/
//ThermalFlowProblemDll.h
///defination of class PassiveScalarProblemDll
#pragma once
#ifdef THERMALFLOWPROBLEMDLL_EXPORTS
#define THERMALFLOWPROBLEMDLL_API __declspec(dllexport)
#else
#define THERMALFLOWPROBLEMDLL_API __declspec(dllimport)
#endif

namespace SIM {

	typedef void* NPtr;

	class ThermalFlowProblemDll2D {
	public:
		static THERMALFLOWPROBLEMDLL_API void Initialize();
		static THERMALFLOWPROBLEMDLL_API void Run();
		static THERMALFLOWPROBLEMDLL_API int Number();
		static THERMALFLOWPROBLEMDLL_API NPtr Type();
		static THERMALFLOWPROBLEMDLL_API NPtr PositionX();
		static THERMALFLOWPROBLEMDLL_API NPtr PositionY();
		static THERMALFLOWPROBLEMDLL_API NPtr VelocityX();
		static THERMALFLOWPROBLEMDLL_API NPtr VelocityY();
		static THERMALFLOWPROBLEMDLL_API NPtr Pressure();
		static THERMALFLOWPROBLEMDLL_API NPtr Temperature();
		static THERMALFLOWPROBLEMDLL_API NPtr Divergence();
		static THERMALFLOWPROBLEMDLL_API void SaveData();
		static THERMALFLOWPROBLEMDLL_API void SensorOut();
	};

}