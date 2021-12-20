#version 430 core 

in float CLIMATE_VALS_VAR;
out vec4 out_Color;
uniform float tMin;
uniform float tMax;


void main(){
	float scalar = (CLIMATE_VALS_VAR - tMin) / (tMax - tMin);
	out_Color = vec4(vec3(scalar), 1.0);
}