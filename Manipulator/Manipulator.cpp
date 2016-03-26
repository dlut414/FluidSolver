// Manipulator.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
///Manipulator.cpp main loop
#pragma once
#include <Header.h>
#include <Bitmap.h>
#include <Controller.h>
#include "../VisualizationDll/VisualizationDll.h"
#include "../PassiveScalarProblemDll/PassiveScalarProblemDll.h"
#include "../ThermalFlowProblemDll/ThermalFlowProblemDll.h"
#include <PreInformation.h>

typedef VIS::VisualizationDll Visualization;
//typedef SIM::PassiveScalarProblemDll Simulation;
typedef SIM::ThermalFlowProblemDll2D Simulation;

static VIS::Controller control;

static void Render() {
	//Visualization::Run(&control, Parameters::Dimension, Simulation::Number(), Simulation::Type(), Simulation::Position(), Simulation::Scalar());
	Visualization::Run(&control, Simulation::Number(), Simulation::Type(), Simulation::PositionX(), Simulation::PositionY(), Simulation::Temperature());
}

static void callBack() {
	if (control.b_save) {
		Simulation::SaveData();
		control.b_save = false;
	}
	if (control.b_sens) {
		Simulation::SensorOut();
		control.b_sens = false;
	}
	if (control.b_bmp) {
		static Bitmap bm;
		static int i = 0;
		char name[256];
		sprintf_s(name, "./out/bm%04d.bmp", i++);
		bm.SaveAsBMP(name);
		control.b_bmp = false;
	}
	if (!control.b_stop) {
		Simulation::Run();
		static int count = 0;
		if (count++ % 50 == 0) {
			control.b_bmp = true;
			control.b_sens = true;
		}
		control.b_dirty = true;
	}
}
static void fps() {

}
static void onMouse(int button, int s, int x, int y) {
	control.clickMouse(button, s, x, y);
	if (button == GLUT_LEFT_BUTTON) {
		const int pickID = Visualization::IntersectColorPick(&control, Simulation::Number(), x, y);
		std::cout << pickID << std::endl;
	}
}
static void onMotion(int x, int y) {
	control.moveMouse(x, y);
	glutPostRedisplay();
}
static void onMouseWheel(int button, int dir, int x, int y) {
	control.rollMouse(button, dir, x, y);
}
static void onReshape(int width, int height) {
	glViewport(0, 0, width, width);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluPerspective(-90.0f, float(control.u_width) / float(control.u_height), 1.0f, 100.0f);
	control.reshapeWindow();
}
static void onKeyboard(unsigned char key, int a, int b) {
	glutPostRedisplay();
	control.pressKey(key, a, b);
	callBack();
}
static void onDisplay() {
	glm::mat4 modelMatrix = glm::translate(glm::mat4(1.0f), control.m_pan)
		* glm::toMat4(control.m_rotation)
		* glm::scale(glm::mat4(1.0f), control.m_scale);

	control.m_modelMat = modelMatrix;
	control.m_viewMat = control.m_camera.GetViewMatrix();
	control.m_viewModelMat = control.m_camera.GetViewMatrix() * modelMatrix;
	control.m_projectionMat = control.m_camera.GetProjectionMatrix();
	control.m_projectionMatInv = glm::inverse(control.m_projectionMat);
	control.m_mvp = control.m_projectionMat * control.m_viewModelMat;
	control.m_mvpInv = glm::inverse(control.m_mvp);

	Render();

	glutSwapBuffers();
	glutReportErrors();

	callBack();

	if (control.b_dirty) {
		glutPostRedisplay();
		control.b_dirty = false;
	}
	if (control.b_leave) {
		glutLeaveMainLoop();
	}
}

static void Initialize(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(control.u_width, control.u_height);
	glutCreateWindow("RTRenderer");
	glutMouseFunc(onMouse);
	glutMotionFunc(onMotion);
	glutMouseWheelFunc(onMouseWheel);
	glutReshapeFunc(onReshape);
	glutKeyboardFunc(onKeyboard);
	glutDisplayFunc(onDisplay);

	Simulation::Initialize();
	Visualization::Initialize();
}

static void Run() {
	glutMainLoop();
}

static void Finalize() {

}


int _tmain(int argc, _TCHAR* argv[]) {
	CreateDirectoryA(std::string(".\\out").c_str(), NULL);
	Initialize(argc, (char**)argv);
	Run();
	Finalize();
	return 0;
}

