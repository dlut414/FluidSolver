/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Published under CC BY-NC
*/
//VisualizationDll.h
///defination of class VisualizationDll
#pragma once
#ifdef VISUALIZATIONDLL_EXPORTS
#define VISUALIZATIONDLL_API __declspec(dllexport)
#else
#define VISUALIZATIONDLL_API __declspec(dllimport)
#endif

#include <DrawParticle.h>

namespace VIS {

	template <typename R>
	class VisualizationDll {
		typedef DrawParticle<R> DP;
		typedef DrawParticle<R>* DPPtr;
	public:
		static VISUALIZATIONDLL_API void Initialize(int argc, char** argv);
		static VISUALIZATIONDLL_API void Run(const int& dim, const int& num, void* tp, void* pos, void* s);
		static VISUALIZATIONDLL_API void Finalize();

	private:
		static void fps();
		static void onMouse(int, int, int, int);
		static void onMotion(int, int);
		static void onMouseWheel(int, int, int, int);
		static void onReshape(int, int);
		static void onKeyboard(unsigned char, int, int);
		static void onDisplay();
		static void callBack();

	private:
		static Controller control;
		static DPPtr drawer;
		static int dimension;
		static int number;
		static void* type; 
		static void* position;
		static void* scalar;
	};

	template <typename R> Controller VisualizationDll<R>::control;
	template <typename R> DrawParticle<R>* VisualizationDll<R>::drawer;
	template <typename R> int VisualizationDll<R>::dimension;
	template <typename R> int VisualizationDll<R>::number;
	template <typename R> void* VisualizationDll<R>::type;
	template <typename R> void* VisualizationDll<R>::position;
	template <typename R> void* VisualizationDll<R>::scalar;

}