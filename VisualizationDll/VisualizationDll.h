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

namespace VIS {

	template <typename R>
	class VisualizationDll {
		typedef DrawParticle<R>* DPPtr;
	public:
		static VISUALIZATIONDLL_API void Initialize(int argc, char** argv);
		static VISUALIZATIONDLL_API void Run();
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
	};

}