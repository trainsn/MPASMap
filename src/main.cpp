#include <glad/glad.h>
#include <GLFW/glfw3.h>
#define _USE_MATH_DEFINES

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/transform2.hpp>

#include <stdlib.h>
#include <cfloat>
#include <netcdf.h>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>
#include "def.h"
#include <assert.h>

#include "shader.h"
#define STBI_MSC_SECURE_CRT
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
using namespace std;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);

// settings
const unsigned int SCR_WIDTH = 1024;
const unsigned int SCR_HEIGHT = 512;

size_t nCells, nEdges, nVertices, nVertLevels, maxEdges, vertexDegree, Time;
vector<double> latVertex, lonVertex, xVertex, yVertex, zVertex;
vector<double> xyzCell, latCell, lonCell; // xCell, yCell, zCell;
vector<int> indexToVertexID, indexToCellID;
vector<int> verticesOnEdge, cellsOnEdge, cellsOnVertex, edgesOnVertex, verticesOnCell;
vector<double> temperature, salinity, thickness;

map<int, int> vertexIndex, cellIndex;

const double max_rho = 6371229.0;
const double layerThickness = 20000.0;
const double eps = 1e-5;

// Base color used for the fog, and clear-to colors.
glm::vec3 base_color(0.0 / 255.0, 0.0 / 255.0, 0.0 / 255.0);

unsigned int VBO, VAO;
unsigned int latCellBuf, latCellTex, lonCellBuf, lonCellTex;
unsigned int latVertexBuf, latVertexTex, lonVertexBuf, lonVertexTex;
unsigned int cellsOnVertexBuf, cellsOnVertexTex;
unsigned int temperatureBuf, temperatureTex;
unsigned int salinityBuf, salinityTex;

void loadMeshFromNetCDF(const string& filename) {
	int ncid;
	int dimid_cells, dimid_edges, dimid_vertices, dimid_vertLevels, dimid_maxEdges,
		dimid_vertexDegree, dimid_Time;
	int varid_latVertex, varid_lonVertex, varid_xVertex, varid_yVertex, varid_zVertex,
		varid_latCell, varid_lonCell, varid_xCell, varid_yCell, varid_zCell,
		varid_verticesOnEdge, varid_cellsOnVertex,
		varid_indexToVertexID, varid_indexToCellID,
		varid_nEdgesOnCell, varid_cellsOncell,
		varid_temperature, varid_salinity, varid_thickness;

	NC_SAFE_CALL(nc_open(filename.c_str(), NC_NOWRITE, &ncid));

	NC_SAFE_CALL(nc_inq_dimid(ncid, "nCells", &dimid_cells));
	NC_SAFE_CALL(nc_inq_dimid(ncid, "nEdges", &dimid_edges));
	NC_SAFE_CALL(nc_inq_dimid(ncid, "nVertices", &dimid_vertices));
	NC_SAFE_CALL(nc_inq_dimid(ncid, "nVertLevels", &dimid_vertLevels));
	NC_SAFE_CALL(nc_inq_dimid(ncid, "maxEdges", &dimid_maxEdges));
	NC_SAFE_CALL(nc_inq_dimid(ncid, "vertexDegree", &dimid_vertexDegree));
	NC_SAFE_CALL(nc_inq_dimid(ncid, "Time", &dimid_Time));

	NC_SAFE_CALL(nc_inq_dimlen(ncid, dimid_cells, &nCells));
	NC_SAFE_CALL(nc_inq_dimlen(ncid, dimid_edges, &nEdges));
	NC_SAFE_CALL(nc_inq_dimlen(ncid, dimid_vertices, &nVertices));
	NC_SAFE_CALL(nc_inq_dimlen(ncid, dimid_vertLevels, &nVertLevels));
	NC_SAFE_CALL(nc_inq_dimlen(ncid, dimid_maxEdges, &maxEdges));
	NC_SAFE_CALL(nc_inq_dimlen(ncid, dimid_vertexDegree, &vertexDegree));
	NC_SAFE_CALL(nc_inq_dimlen(ncid, dimid_Time, &Time));

	NC_SAFE_CALL(nc_inq_varid(ncid, "indexToVertexID", &varid_indexToVertexID));
	NC_SAFE_CALL(nc_inq_varid(ncid, "indexToCellID", &varid_indexToCellID));
	NC_SAFE_CALL(nc_inq_varid(ncid, "latCell", &varid_latCell));
	NC_SAFE_CALL(nc_inq_varid(ncid, "lonCell", &varid_lonCell));
	NC_SAFE_CALL(nc_inq_varid(ncid, "xCell", &varid_xCell));
	NC_SAFE_CALL(nc_inq_varid(ncid, "yCell", &varid_yCell));
	NC_SAFE_CALL(nc_inq_varid(ncid, "zCell", &varid_zCell));
	NC_SAFE_CALL(nc_inq_varid(ncid, "latVertex", &varid_latVertex));
	NC_SAFE_CALL(nc_inq_varid(ncid, "lonVertex", &varid_lonVertex));
	NC_SAFE_CALL(nc_inq_varid(ncid, "xVertex", &varid_xVertex));
	NC_SAFE_CALL(nc_inq_varid(ncid, "yVertex", &varid_yVertex));
	NC_SAFE_CALL(nc_inq_varid(ncid, "zVertex", &varid_zVertex));
	NC_SAFE_CALL(nc_inq_varid(ncid, "verticesOnEdge", &varid_verticesOnEdge));
	NC_SAFE_CALL(nc_inq_varid(ncid, "cellsOnVertex", &varid_cellsOnVertex));
	NC_SAFE_CALL(nc_inq_varid(ncid, "nEdgesOnCell", &varid_nEdgesOnCell));
	NC_SAFE_CALL(nc_inq_varid(ncid, "cellsOnCell", &varid_cellsOncell));
	NC_SAFE_CALL(nc_inq_varid(ncid, "temperature", &varid_temperature));
	NC_SAFE_CALL(nc_inq_varid(ncid, "salinity", &varid_salinity));
	NC_SAFE_CALL(nc_inq_varid(ncid, "layerThickness", &varid_thickness));

	const size_t start_cells[1] = { 0 }, size_cells[1] = { nCells };

	latCell.resize(nCells);
	lonCell.resize(nCells);
	indexToCellID.resize(nCells);
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_latCell, start_cells, size_cells, &latCell[0]));
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_lonCell, start_cells, size_cells, &lonCell[0]));
	NC_SAFE_CALL(nc_get_vara_int(ncid, varid_indexToCellID, start_cells, size_cells, &indexToCellID[0]));
	for (int i = 0; i < nCells; i++) {
		cellIndex[indexToCellID[i]] = i;
		// fprintf(stderr, "%d, %d\n", i, indexToCellID[i]);
	}

	std::vector<double> coord_cells;
	coord_cells.resize(nCells);
	xyzCell.resize(nCells * 3);
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_xCell, start_cells, size_cells, &coord_cells[0]));
	for (int i = 0; i < nCells; i++)
		xyzCell[i * 3] = coord_cells[i];
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_yCell, start_cells, size_cells, &coord_cells[0]));
	for (int i = 0; i < nCells; i++)
		xyzCell[i * 3 + 1] = coord_cells[i];
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_zCell, start_cells, size_cells, &coord_cells[0]));
	for (int i = 0; i < nCells; i++)
		xyzCell[i * 3 + 2] = coord_cells[i];

	//for (int i = 0; i < nCells; i++) {
	//	double x = max_rho * cos(latCell[i]) * cos(lonCell[i]);
	//	double y = max_rho * cos(latCell[i]) * sin(lonCell[i]);
	//	double z = max_rho * sin(latCell[i]);
	//	assert(abs(x - xyzCell[i * 3]) < eps && abs(y - xyzCell[i * 3 + 1]) < eps && abs(z - xyzCell[i * 3 + 2])< eps);
	//}

	const size_t start_vertices[1] = { 0 }, size_vertices[1] = { nVertices };
	latVertex.resize(nVertices);
	lonVertex.resize(nVertices);
	xVertex.resize(nVertices);
	yVertex.resize(nVertices);
	zVertex.resize(nVertices);
	indexToVertexID.resize(nVertices);

	NC_SAFE_CALL(nc_get_vara_int(ncid, varid_indexToVertexID, start_vertices, size_vertices, &indexToVertexID[0]));
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_latVertex, start_vertices, size_vertices, &latVertex[0]));
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_lonVertex, start_vertices, size_vertices, &lonVertex[0]));
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_xVertex, start_vertices, size_vertices, &xVertex[0]));
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_yVertex, start_vertices, size_vertices, &yVertex[0]));
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_zVertex, start_vertices, size_vertices, &zVertex[0]));

	//for (int i = 0; i < nVertices; i++) {
	//	double x = max_rho * cos(latVertex[i]) * cos(lonVertex[i]);
	//	double y = max_rho * cos(latVertex[i]) * sin(lonVertex[i]);
	//	double z = max_rho * sin(latVertex[i]);
	//	assert(abs(x - xVertex[i]) < eps && abs(y - yVertex[i]) < eps && abs(z - zVertex[i]) < eps);
	//}

	for (int i = 0; i < nVertices; i++) {
		vertexIndex[indexToVertexID[i]] = i;
		// fprintf(stderr, "%d, %d\n", i, indexToVertexID[i]);
	}

	const size_t start_edges2[2] = { 0, 0 }, size_edges2[2] = { nEdges, 2 };
	verticesOnEdge.resize(nEdges * 2);

	NC_SAFE_CALL(nc_get_vara_int(ncid, varid_verticesOnEdge, start_edges2, size_edges2, &verticesOnEdge[0]));

	//for (int i=0; i<nEdges; i++) 
	//   fprintf(stderr, "%d, %d\n", verticesOnEdge[i*2], verticesOnEdge[i*2+1]);

	const size_t start_vertex_cell[2] = { 0, 0 }, size_vertex_cell[2] = { nVertices, 3 };
	cellsOnVertex.resize(nVertices * 3);

	NC_SAFE_CALL(nc_get_vara_int(ncid, varid_cellsOnVertex, start_vertex_cell, size_vertex_cell, &cellsOnVertex[0]));

	const size_t start_time_cell_vertLevel[3] = { 0, 0, 0 }, size_time_cell_vertLevel[3] = { Time, nCells, nVertLevels };
	temperature.resize(Time * nCells * nVertLevels);
	salinity.resize(Time * nCells * nVertLevels);
	thickness.resize(Time * nCells * nVertLevels);

	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_temperature, start_time_cell_vertLevel, size_time_cell_vertLevel, &temperature[0]));
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_salinity, start_time_cell_vertLevel, size_time_cell_vertLevel, &salinity[0]));
	NC_SAFE_CALL(nc_get_vara_double(ncid, varid_thickness, start_time_cell_vertLevel, size_time_cell_vertLevel, &thickness[0]));

	NC_SAFE_CALL(nc_close(ncid));

	fprintf(stderr, "%zu, %zu, %zu, %zu\n", nCells, nEdges, nVertices, nVertLevels);
}

void initBuffers() {
	// set up vertex data (and buffer(s)) and configure vertex attributes
	// ------------------------------------------------------------------
	glGenBuffers(1, &VBO);
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, nVertices * sizeof(int), &indexToVertexID[0], GL_STATIC_DRAW);
	//glBufferData(GL_ARRAY_BUFFER, nCells * sizeof(int), &indexToCellID[0], GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribIPointer(0, 1, GL_INT, sizeof(int), 0);
	glBindVertexArray(0);
}

void initTextures() {
	//// Coordinates of cells
	glGenBuffers(1, &latCellBuf);
	glBindBuffer(GL_TEXTURE_BUFFER, latCellBuf);
	vector<float> latCellFloat(latCell.begin(), latCell.end());
	glBufferData(GL_TEXTURE_BUFFER, nCells * sizeof(float), &latCellFloat[0], GL_STATIC_DRAW);
	glGenTextures(1, &latCellTex);

	glGenBuffers(1, &lonCellBuf);
	glBindBuffer(GL_TEXTURE_BUFFER, lonCellBuf);
	vector<float> lonCellFloat(lonCell.begin(), lonCell.end());
	glBufferData(GL_TEXTURE_BUFFER, nCells * sizeof(float), &lonCellFloat[0], GL_STATIC_DRAW);
	glGenTextures(1, &lonCellTex);

	glGenBuffers(1, &cellsOnVertexBuf);
	glBindBuffer(GL_TEXTURE_BUFFER, cellsOnVertexBuf);
	glBufferData(GL_TEXTURE_BUFFER, nVertices * 3 * sizeof(int), &cellsOnVertex[0], GL_STATIC_DRAW);
	glGenTextures(1, &cellsOnVertexTex);

	// Coordinates of vertices 
	glGenBuffers(1, &latVertexBuf);
	glBindBuffer(GL_TEXTURE_BUFFER, latVertexBuf);
	vector<float> latVertexFloat(latVertex.begin(), latVertex.end());
	glBufferData(GL_TEXTURE_BUFFER, nVertices * sizeof(float), &latVertexFloat[0], GL_STATIC_DRAW);
	glGenTextures(1, &latVertexTex);

	glGenBuffers(1, &lonVertexBuf);
	glBindBuffer(GL_TEXTURE_BUFFER, lonVertexBuf);
	vector<float> lonVertexFloat(lonVertex.begin(), lonVertex.end());
	glBufferData(GL_TEXTURE_BUFFER, nVertices * sizeof(float), &lonVertexFloat[0], GL_STATIC_DRAW);
	glGenTextures(1, &lonVertexTex);

	// temperature 
	glGenBuffers(1, &temperatureBuf);
	glBindBuffer(GL_TEXTURE_BUFFER, temperatureBuf);
	vector<float> temperatureFloat(temperature.begin(), temperature.end());
	glBufferData(GL_TEXTURE_BUFFER, Time * nCells * nVertLevels * sizeof(float), &temperatureFloat[0], GL_STATIC_DRAW);
	glGenTextures(1, &temperatureTex);

	// salinity 
	glGenBuffers(1, &salinityBuf);
	glBindBuffer(GL_TEXTURE_BUFFER, salinityBuf);
	vector<float> salinityFloat(salinity.begin(), salinity.end());
	glBufferData(GL_TEXTURE_BUFFER, Time * nCells * nVertLevels * sizeof(float), &salinityFloat[0], GL_STATIC_DRAW);
	glGenTextures(1, &salinityTex);
}

int main(int argc, char **argv)
{
    char filename[1024];
	sprintf(filename, argv[1]);
	fprintf(stderr, "%s\n", argv[1]); 
	
	string filename_s = filename;
	int pos_last_dot = filename_s.rfind(".");
	string input_name = filename_s.substr(0, pos_last_dot);
    int pos_first_dash = filename_s.find("_");
	string fileid = filename_s.substr(0, pos_first_dash);
	
	char input_path[1024];
	sprintf(input_path, "/fs/project/PAS0027/MPAS1/Inter/%s", filename);
	//loadMeshFromNetCDF("D:\\OSU\\Grade1\\in-situ\\6.0\\output.nc");
	//loadMeshFromNetCDF("D:\\OSU\\Grade1\\in-situ\\MPAS-server\\Results\\0070_4.88364_578.19012_0.51473_227.95909_ght0.2_epoch420.nc");
	loadMeshFromNetCDF(input_path);

	// glfw: initialize and configure
	// ------------------------------
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this statement to fix compilation on OS X
#endif

	// glfw window creation
	// --------------------
	GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "MPASMap", NULL, NULL);
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// glad: load all OpenGL function pointers
	// ---------------------------------------
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to initialize GLAD" << std::endl;
		return -1;
	}

	// configure global opengl state
	// -----------------------------
	glEnable(GL_DEPTH_TEST);
	

	// build and compile shaders
	// -------------------------
	//Shader shader("triangle_mesh.vs", "dvr.fs", "dvr.gs");
	//Shader shader("triangle_center.vs", "triangle_center.fs");
	//Shader shader("hexagon_center.vs", "hexagon_center.fs");
	Shader shader("../res/shaders/triangle_mesh.vs", "../res/shaders/triangle_mesh.fs", "../res/shaders/triangle_mesh.gs");

	initTextures();
	initBuffers();

	// render loop
	// -----------
	// while (!glfwWindowShouldClose(window))
	for (int layer_id = 0; layer_id < nVertLevels; layer_id++)
	{
		// render
		// ------
		glClearColor(base_color[0], base_color[1], base_color[2], 1.0); // Set the WebGL background color.
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// draw points
		shader.use();
		
		/*glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_BUFFER, latVertexTex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, latVertexBuf);
		shader.setInt("latVertex", 0);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_BUFFER, lonVertexTex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, lonVertexBuf);
		shader.setInt("lonVertex", 1);*/

		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_BUFFER, latCellTex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, latCellBuf);
		shader.setInt("latCell", 0);

		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_BUFFER, lonCellTex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, lonCellBuf);
		shader.setInt("lonCell", 1);

		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_BUFFER, cellsOnVertexTex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32I, cellsOnVertexBuf);
		shader.setInt("cellsOnVertex", 2);

		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_BUFFER, temperatureTex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, temperatureBuf);
		shader.setInt("temperature", 3);

		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_BUFFER, salinityTex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_R32F, salinityBuf);
		shader.setInt("salinity", 4);

		shader.setInt("TOTAL_LAYERS", nVertLevels);
		shader.setInt("layer_id", layer_id);
		float tMin = -1.93f, tMax = 30.35f;
		/*for (int i = 0; i < temperature.size(); i += nVertLevels) {
			if (temperature[i] < tMin)
				tMin = temperature[i];
			if (temperature[i] > tMax)
				tMax = temperature[i];
		}*/
		shader.setFloat("tMin", tMin);
		shader.setFloat("tMax", tMax);

		glBindVertexArray(VAO);
		glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
		glDrawArrays(GL_POINTS, 0, nVertices);

		stbi_flip_vertically_on_write(1);
		char imagepath[1024];
		sprintf(imagepath, "/fs/project/PAS0027/MPAS1/Inter/%s/layer%d.png", fileid.c_str(), layer_id);
		float* pBuffer = new float[SCR_WIDTH * SCR_HEIGHT * 4];
		unsigned char* pImage = new unsigned char[SCR_WIDTH * SCR_HEIGHT * 3];
		glReadBuffer(GL_BACK);
		glReadPixels(0, 0, SCR_WIDTH, SCR_HEIGHT, GL_RGBA, GL_FLOAT, pBuffer);
		for (unsigned int j = 0; j < SCR_HEIGHT; j++) {
			for (unsigned int k = 0; k < SCR_WIDTH; k++) {
				int index = j * SCR_WIDTH + k;
				pImage[index * 3 + 0] = GLubyte(min(pBuffer[index * 4 + 0] * 255, 255.0f));
				pImage[index * 3 + 1] = GLubyte(min(pBuffer[index * 4 + 1] * 255, 255.0f));
				pImage[index * 3 + 2] = GLubyte(min(pBuffer[index * 4 + 2] * 255, 255.0f));
			}
		}
		stbi_write_png(imagepath, SCR_WIDTH, SCR_HEIGHT, 3, pImage, SCR_WIDTH * 3);
		delete pBuffer;
		delete pImage;

		// glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
		// -------------------------------------------------------------------------------
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	// optional: de-allocate all resources once they've outlived their purpose:
	// ------------------------------------------------------------------------
	glDeleteVertexArrays(1, &VAO);
	glDeleteBuffers(1, &VBO);

	glfwTerminate();
	return 0;
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	// make sure the viewport matches the new window dimensions; note that width and 
	// height will be significantly larger than specified on retina displays.
	glViewport(0, 0, width, height);
}
