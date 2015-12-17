/*
 * LICENCE
 * copyright 2014 ~ ****
 * Some rights reserved.
 * Author: HUFANGYUAN
 * Released under CC BY-NC
*/
#pragma once
#include "Header.h"
#include "Draw_.h"
#include "Shader.h"
#include <Controller.h>
#include <BBox.h>

namespace VIS {

	template <typename R>
	class DrawParticle : public Draw_ {
		template <typename R_ = R>	struct DataType {};
		template <>					struct DataType<float>	{ enum { value = GL_FLOAT, }; };
		template <>					struct DataType<double>	{ enum { value = GL_DOUBLE, }; };
	public:
		DrawParticle() : Draw_() { init(); }
		~DrawParticle() {}

		template <typename I, typename V, typename T>
		void draw(const Controller* const controlPtr, const int& dim, const int& num, const I* const type, const V* const vert, const T* const s1) const {
			///clear framebuffer
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) printf("fbo_def not ready\n");
			glViewport(0, 0, controlPtr->u_width, controlPtr->u_height);
			glFrontFace(GL_CCW);
			glClearDepth(1);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			
			///use program0
			glBindFramebuffer(GL_FRAMEBUFFER, 0);
			if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) printf("0 not ready\n");
			glViewport(0, 0, controlPtr->u_width, controlPtr->u_height);
			glFrontFace(GL_CCW);
			glDepthFunc(GL_LESS);

			glClearDepth(1);
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glUseProgram(shaderObj.programID[0]);
			
			glUniform1i(shaderObj.intID[0], controlPtr->i_visFlag);
			glUniform1f(shaderObj.floatID[0], controlPtr->f_visRange);
			glUniformMatrix4fv(shaderObj.matrixID[0], 1, GL_FALSE, &(controlPtr->m_mvp[0][0]));
			glUniformMatrix4fv(shaderObj.matrixID[1], 1, GL_FALSE, &(controlPtr->m_mvpInv[0][0]));
			glUniformMatrix4fv(shaderObj.matrixID[2], 1, GL_FALSE, &(controlPtr->m_modelMat[0][0]));
			glUniformMatrix4fv(shaderObj.matrixID[3], 1, GL_FALSE, &(controlPtr->m_viewMat[0][0]));
			glUniformMatrix4fv(shaderObj.matrixID[4], 1, GL_FALSE, &(controlPtr->m_projectionMat[0][0]));

			glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
			glBufferData(GL_ARRAY_BUFFER, num*sizeof(int), type, GL_STATIC_DRAW);
			glVertexAttribPointer(0, 1, GL_FLOAT, GL_FALSE, 0, (void*)0);

			glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
			glBufferData(GL_ARRAY_BUFFER, dim*num*sizeof(R), vert, GL_STATIC_DRAW);
			glVertexAttribPointer(1, dim, DataType<>::value, GL_FALSE, 0, (void*)0);

			glBindBuffer(GL_ARRAY_BUFFER, vbo[2]);
			glBufferData(GL_ARRAY_BUFFER, num*sizeof(R), s1, GL_STATIC_DRAW);
			glVertexAttribPointer(2, 1, DataType<>::value, GL_FALSE, 0, (void*)0);

			glEnableVertexAttribArray(0);
			glEnableVertexAttribArray(1);
			glEnableVertexAttribArray(2);

			glDrawArrays(GL_POINTS, 0, (GLsizei)num);

			glDisableVertexAttribArray(0);
			glDisableVertexAttribArray(1);
			glDisableVertexAttribArray(2);
		}

	public:
		GLfloat r_bg_chessboard[36];

		GLuint bgTex;
		GLuint frontFaceTex;
		GLuint dirTex;
		GLuint dbo[2];

		GLuint volNormTex;
		GLuint volDataTex;
		GLuint cellInFluidTex;

		Shader shaderObj;

		BBox<float> box;
		GLuint ix, iy, iz;
		GLuint icx, icy, icz;

		void* volNorm;
		void* volData;
		void* cellInFluid;

	private:
		void init() {
			Draw_::initVAO(Draw_::vao);
			Draw_::initVBO(Draw_::vbo[0]);
			Draw_::initVBO(Draw_::vbo[1]);
			Draw_::initVBO(Draw_::vbo[2]);
			Draw_::initVBO(Draw_::vbo[3]);
			Draw_::initVBO(Draw_::vbo[4]);
			Draw_::initVBO(Draw_::vbo[5]);
			initShader();
		}
		void initShader() {
			shaderObj.programID.push_back(shaderObj.LoadShader());
			//shaderObj.programID.push_back(shaderObj.LoadShader("../VisualizationDll/shader0/vertex.glsl", "../VisualizationDll/shader0/fragment.glsl"));
			shaderObj.matrixID.push_back( glGetUniformLocation(shaderObj.programID[0], "vMvp") );
			shaderObj.matrixID.push_back( glGetUniformLocation(shaderObj.programID[0], "fMvpInv") );
			shaderObj.matrixID.push_back( glGetUniformLocation(shaderObj.programID[0], "vModelMat") );
			shaderObj.matrixID.push_back( glGetUniformLocation(shaderObj.programID[0], "vViewMat") );
			shaderObj.matrixID.push_back( glGetUniformLocation(shaderObj.programID[0], "vProjectionMat") );

			shaderObj.intID.push_back(glGetUniformLocation(shaderObj.programID[0], "flag"));
			shaderObj.floatID.push_back(glGetUniformLocation(shaderObj.programID[0], "range"));
		}
	};

}
