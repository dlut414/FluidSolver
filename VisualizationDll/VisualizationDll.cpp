// Visualization.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "VisualizationDll.h"
#include <Header.h>
#include <Controller.h>
#include <DrawParticle.h>

namespace VIS {

	template <typename R>
	void VisualizationDll<R>::Initialize(int argc, char** argv) {
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH);
		glutInitWindowPosition(0, 0);
		glutInitWindowSize(control.u_width, control.u_height);
		glutCreateWindow("RTRenderer");
		//glLightModeli           (GL_LIGHT_MODEL_TWO_SIDE, 1);
		glutMouseFunc(onMouse);
		glutMotionFunc(onMotion);
		glutMouseWheelFunc(onMouseWheel);
		glutReshapeFunc(onReshape);
		glutKeyboardFunc(onKeyboard);
		glutDisplayFunc(onDisplay);

		glEnable(GL_TEXTURE_1D);
		glEnable(GL_TEXTURE_2D);
		glEnable(GL_TEXTURE_3D);
		glEnable(GL_CULL_FACE);
		//glDisable               (GL_CULL_FACE);
		glFrontFace(GL_CCW);
		glEnable(GL_POINT_SPRITE_ARB);
		glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LESS);
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER, 0.f);
		//glEnable                (GL_BLEND);
		//glBlendFunc             (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glClearColor(1.f, 1.f, 1.f, 0.f);
		glewInit();
		drawer = new DP(&control);
	}

	template <typename R> 
	void VisualizationDll<R>::Run(const int& dim, const int& num, void* tp, void* pos, void* s) {
		dimension = dim; number = num; type = tp; position = pos; scalar = s;
		glutMainLoop();
	}

	template <typename R> 
	void VisualizationDll<R>::Finalize() {

	}

	template <typename R> void VisualizationDll<R>::fps() {

	}
	template <typename R> void VisualizationDll<R>::onMouse(int button, int s, int x, int y) {
		control.clickMouse(button, s, x, y);
	}
	template <typename R> void VisualizationDll<R>::onMotion(int x, int y) {
		control.moveMouse(x, y);
		glutPostRedisplay();
	}
	template <typename R> void VisualizationDll<R>::onMouseWheel(int button, int dir, int x, int y) {
		control.rollMouse(button, dir, x, y);
	}
	template <typename R> void VisualizationDll<R>::onReshape(int width, int height) {
		glViewport(0, 0, width, width);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluPerspective(-90.0f, float(control.u_width) / float(control.u_height), 1.0f, 100.0f);
		control.reshapeWindow();
	}
	template <typename R> void VisualizationDll<R>::onKeyboard(unsigned char key, int a, int b) {
		glutPostRedisplay();
		control.pressKey(key, a, b);
		callBack();
	}
	template <typename R> void VisualizationDll<R>::onDisplay() {
		glm::mat4 modelMatrix = glm::translate(glm::mat4(1.0f), control.m_pan)
			* glm::toMat4(control.m_rotation)
			* glm::scale(glm::mat4(1.0f), control.m_scale);

		control.m_viewModelMat = control.m_camera.GetViewMatrix() * modelMatrix;
		control.m_projectionMat = control.m_camera.GetProjectionMatrix();
		control.m_projectionMatInv = glm::inverse(control.m_projectionMat);
		control.m_mvp = control.m_projectionMat * control.m_viewModelMat;
		control.m_mvpInv = glm::inverse(control.m_mvp);

		drawer->draw(dimension, number, type, position, scalar);

		glutSwapBuffers();
		glutReportErrors();

		callBack();

		fps();

		if (control.b_dirty) {
			glutPostRedisplay();
			control.b_dirty = false;
		}
		if (control.b_leave) {
			glutLeaveMainLoop();
		}
	}
	template <typename R> void VisualizationDll<R>::callBack() {

	}

	template class VisualizationDll<float>;
	template class VisualizationDll<double>;

}


