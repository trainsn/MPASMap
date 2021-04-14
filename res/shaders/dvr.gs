#version 430 core
layout (points) in;
layout (triangle_strip, max_vertices = 3) out;

in V2G{
	flat int vertex_id;
}v2g[];

out G2F{
	flat int triangle_id;
	smooth vec3 o_pos;
}g2f;

//3D transformation matrices
uniform mat4 uMVMatrix;
uniform mat4 uPMatrix;
//connectivity
uniform samplerBuffer latCell;
uniform samplerBuffer lonCell;
uniform isamplerBuffer cellsOnVertex;

uniform float GLOBLE_RADIUS;

void main(void){
	//only consider the outter most surface mesh 
	int vertexId = v2g[0].vertex_id;
	//query three neighbor cell id 
	ivec3 cellId3 = texelFetch(cellsOnVertex, vertexId).xyz;
	int cellIds[3];
	cellIds[0] = cellId3.x;
	cellIds[1] = cellId3.y;
	cellIds[2] = cellId3.z;

	for (int i = 0; i < 3; i++){
		int cell_id = cellIds[i];
		
		float lat = texelFetch(latCell, cell_id - 1).x;
		float lon = texelFetch(lonCell, cell_id - 1).x; // indexCell: cell_id - 1 	
		vec4 xyzw = vec4(
						GLOBLE_RADIUS * cos(lat) * cos(lon),
						GLOBLE_RADIUS * cos(lat) * sin(lon),
						GLOBLE_RADIUS * sin(lat),
						1.0);

		gl_Position = uPMatrix * uMVMatrix * xyzw;
		g2f.o_pos = xyzw.xyz;
		g2f.triangle_id = vertexId;
		EmitVertex();
	}
	EndPrimitive();
}