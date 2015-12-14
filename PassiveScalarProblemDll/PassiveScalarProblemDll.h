/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Published under CC BY-NC
*/
//PassiveScalarProblemDll.h
///defination of class PassiveScalarProblemDll
#pragma once
#ifdef PASSIVESCALARPROBLEMDLL_EXPORTS
#define PASSIVESCALARPROBLEMDLL_API __declspec(dllexport)
#else
#define PASSIVESCALARPROBLEMDLL_API __declspec(dllimport)
#endif

#include <FractionalStep_x.h>

namespace SIM {

	struct Parameters {
		typedef double DataType;
		typedef DataType* DataTypePtr;
		enum { Dim = 2, };
		enum { Order = 2, };
	};
	typedef FractionalStep_x<Parameters::DataType, Parameters::Dim, Parameters::Order> FS;
	typedef FS* FSPtr;
	typedef void* NPtr;

	class PassiveScalarProblemDll {
	public:
		static PASSIVESCALARPROBLEMDLL_API void Initialize();
		static PASSIVESCALARPROBLEMDLL_API void Run();
		static PASSIVESCALARPROBLEMDLL_API NPtr Position();
		static PASSIVESCALARPROBLEMDLL_API NPtr Scalar();
	};

}