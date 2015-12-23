/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
#define GLSL(version, shader)  "#version " #version "\n" #shader

const GLchar* vertex = GLSL(330 core,

uniform mat4 vMvp;
uniform mat4 vModelMat;
uniform mat4 vViewMat;
uniform mat4 vProjectionMat;

layout(location = 0) in int vType;
layout(location = 1) in vec3 vPos;
layout(location = 2) in float v3;
layout(location = 3) in float v4;

flat out int fType;
out vec4 fPos;
out float f3;
out float f4;

void main() {
	vec3 pos = vec3(vPos.x, vPos.y, 0.2f*v3);
	fType = vType;
	fPos = vMvp * vec4(pos, 1.0f);

	f3 = v3;
	f4 = v4;

	gl_Position = fPos;
	gl_PointSize = 4.;// +scal / 5.e10;
}

);