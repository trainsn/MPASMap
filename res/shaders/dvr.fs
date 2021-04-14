#version 430 core

out vec4 out_Color;

uniform int nLayers;
uniform int maxEdges;
uniform float threshold;

in G2F{
	flat int triangle_id;
	smooth vec3 o_pos;
}

uniform mat4 inv_mvm;
uniform float stepsize;

//connectivity
uniform samplerBuffer latCell;
uniform samplerBuffer lonCell;
uniform isamplerBuffer cellsOnVertex;
uniform isamplerBuffer edgesOnVertex;
uniform isamplerBuffer cellsOnEdge;
uniform isamplerBuffer verticesOnEdge;
uniform isamplerBuffer verticesOnCell;
uniform isamplerBuffer nEdgesOnCell;
uniform isamplerBuffer maxLevelCell;
uniform samplerBuffer cellVal;
//symbolic definition of dual triangle mesh
#define CLIMATE_VALS_VAR  cellVal
#define TRIANGLE_TO_EDGES_VAR edgesOnVertex
#define EDGE_CORNERS_VAR cellsOnEdge
#define CORNERS_LAT_VAR latCell
#define CORNERS_LON_VAR lonCell
#define FACEID_TO_EDGEID faceId_to_edgeId
#define EDGE_TO_TRIANGLES_VAR verticesOnEdge
#define CORNER_TO_TRIANGLES_VAR verticesOnCell
#define CORNER_TO_TRIANGLES_DIMSIZES nEdgesOnCell

#define FLOAT_MAX  3.402823466e+38F
#define FLOAT_MIN -2.402823466e+36F
uniform float GLOBLE_RADIUS;
uniform float THICKNESS;
uniform int TOTAL_LAYERS;

int	d_mpas_faceCorners[24] = {
    0, 1, 2,  3, 4, 5,//top 0 and bottom 1
    4, 2, 1,  4, 5, 2,//front 2,3
    5, 0, 2,  5, 3, 0,//right 4,5
    0, 3, 1,  3, 4, 1//left 6,7
};

struct Ray(){
	vec3 o;
	vec3 d;
};

struct HitRec{
	float t;	//  t value along the hitted face 
	int hitFaceid;	// hitted face id 
	int nextlayerId;
} 

struct MPASPrism {
	uint m_prismId;
	int m_iLayer;
	vec3 vtxCoordTop[3]; // 3 top vertex coordinates
	vec3 vtxCoordBottom[3]; // 3 botton vertex coordinates 
	int m_idxEdge[3];
	int idxVtx[3]; // triangle vertex index is equvalent to hexagon cell index 
};

#define FLOAT_ERROR 1.0E-8

bool rayIntersectsTriangle(vec3 p, vec3 d,
	vec3 v0, vec3 v1, vec3 v2, inout float t)
{
	vec3 e1, e2, h, s, q;
	float a, f, u, v;
	//float error = 1.0e-4;//0.005f;
	e1 = v1 - v0;
	e2 = v2 - v0;
	//crossProduct(h, d, e2);
	h = cross(d, e2);
	a = dot(e1, h);//innerProduct(e1, h);

	if (a > -FLOAT_ERROR && a < FLOAT_ERROR)
		return(false);

	f = 1.0f / a;
	s = p - v0;//_vector3d(s, p, v0);
	u = f * dot(s, h);//(innerProduct(s, h));

	if (u < -FLOAT_ERROR || u >(1.0f + FLOAT_ERROR))
		return(false);

	q = cross(s, e1);//crossProduct(q, s, e1);
	v = f * dot(d, q);//innerProduct(d, q);

	if (v < -FLOAT_ERROR || u + v >(1.0f + FLOAT_ERROR))
		return(false);

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	t = f * dot(e2, q);//innerProduct(e2, q);

	if (t > FLOAT_ERROR)//ray intersection
		return(true);
	else // this means that there is a line intersection
		// but not a ray intersection
		return (false);
}

bool ReloadVtxInfo(in int triangle_id, in int iLayer, inout MPASPrism prism) {
	prism.m_prismId = triangle_id;
	prism.m_iLayer = iLayer;

	prism.idxVtx[0] = -1;
	prism.idxVtx[1] = -1;
	prism.idxVtx[2] = -1;

	// load first edge 
	// index of first edge of triangle
	ivec3 idxEdges = texelFetch(TRIANGLE_TO_EDGES_VAR, triangle_id).xyz;
	int idxEdge = idxEdges.x; // TRIANGLE_TO_EDGES_VAR[m_prismId * 3 + 0];
	prism.m_idxEdge[0] = idxEdge;
	// index of start corner of this edge
	ivec2 cornerIdxs = texelFetch(EDGE_CORNERS_VAR, idxEdge).xy;
	int iS = cornerIdxs.x; //EDGE_CORNERS_VAR[idxEdge * 2];	
	// index of end corner of this edge 
	int iE = cornerIdxs.y; //EDGE_CORNERS_VAR[idxEdge * 2 + 1];
	prism.idxVtx[0] = iS;
	prism.idxVtx[1] = iE;
	int edge1E = iE;

	// load second edge 
	int idxEdge = idxEdges.y; // TRIANGLE_TO_EDGES_VAR[m_prismId * 3 + 1];
	prism.m_idxEdge[2] = idxEdge;
	// index of start corner of second edge
	ivec2 cornerIdxs = texelFetch(EDGE_CORNERS_VAR, idxEdge).xy;
	int iS = cornerIdxs.x; //EDGE_CORNERS_VAR[idxEdge * 2];	
	// index of end corner of this edge 
	int iE = cornerIdxs.y; //EDGE_CORNERS_VAR[idxEdge * 2 + 1];

	bool normalCase = (edge1E != iS) && (edge1E != iE); // the second edge connects corner 0 and corner 2
	if (iS != prism.idxVtx[0] && iS != prism.idxVtx[1]) {	// find the index of the third corner.
		prism.idxVtx[2] = iS;
	}
	else {
		prism.idxVtx[2] = iE;
	}

	// index of third edge.
	int idxEdge = idxEdges.z; // TRIANGLE_TO_EDGES_VAR[m_prismId * 3 + 1];
	prism.m_idxEdge[2] = idxEdge;
	if (!normalCase) {
		// swap m_idxEdge[1] and m_idexEdge[2]
		prism.m_idxEdge[1] ^= prism.m_idxEdge[2];
		prism.m_idxEdge[2] ^= prism.m_idxEdge[1];
		prism.m_idxEdge[1] ^= prism.m_idxEdge[2];
	}

	// load vertex info based on edge's info
	float lon[3], lat[3]; // longtitude and latitude of three corners
	float maxR = GLOBAL_RADIUS - THICKNESS * iLayer; // Radius of top triangle in current cell
	float minR = GLOBAL_RADIUS - THICKNESS; // Radius of bottom triangle in current cell

	for (int i = 0; i < 3; i++) {	// for each corner index of the triangle (specified by prismId)
		int idxCorner = prism.idxVtx[i];	//TRIANGLE_TO_CORNERS_VAR[m_prismId*3+i];
		lat[i] = texelFetch(CORNERS_LAT_VAR, idxCorner - 1).r;	// CORNERS_LAT_VAR[idxCorner]
		lon[i] = texelFetch(CORNERS_LON_VAR, idxCorner - 1).r;	// CORNERS_LON_VAR[idxCorner]
		prism.vtxCoordTop[i] = vec3(maxR*cos(lat[i])*cos(lon[i]), maxR*cos(lat[i])*sin(lon[i]), maxR*sin(lat[i]));
		prism.vtxCoordBottom[i] = vec3(minR*cos(lat[i])*cos(lon[i]), minR*cos(lat[i])*sin(lon[i]), minR*sin(lat[i]));
	}

	return true; 
}

// Return adjacent mesh cell (which is triangle in the remeshed MPAS mesh) id 
// which shares with current triangle (specified by curTriangleId) the edge belongs to the face (denoted by faceId)
int getAdjacentCellId(inout MPASPrism prism, int faceId) {
	if (prism.m_iLayer == TOTAL_LAYERS - 2 && faceId == 1) {
		return -1; // we reached the deepest layers, no more layers beyond current layer
	} 
	else if (faceId == 0 || faceId == 1) {
		return int(prism.m_prismId);	// currentTriangleId
	}

	const int faceId_to_edgeId[8] = { -1, -1, 2, 2, 1, 1, 0, 0 };
	int edgeId = FACEID_TO_EDGEID[faceId];
	int idxEdge = prism.m_idxEdge[edgeId];
	ivec2 nextTriangleIds = texelFetch(EDGES_TO_TRIANGLES_VAR, idxEdge).xy;	// EDGE_TO_TRIANGLES_VAR[idxEdge * 2];
	if (nextTriangleIds.x == prism.m_prismId){ // curTriangleId
		return nextTriangleIds.y;
	}
	return nextTriangleIds.x;
}

#define FACEID_TO_CORNERID d_gcrm_faceCorners

int rayPrismIntersection(inout MPASPrism prism, in Ray r, inout HitRec tInRec,
	inout HitRec tOutRec, inout int nextCellId) {
	nextCellId = -1;	// assume no netxt prism to shot into
	int nHit = 0;
	int nFaces = 8;
	tOutRec.hitFace = -1; // initialize to tOutRec
	tOutRec.t = -1.0f;
	double min_t = FLOAT_MAX, max_t = -1.0f;
	vec3 vtxCoord[6];
	vtxCoord[0] = prism.vtxCoordTop[0];
	vtxCoord[1] = prism.vtxCoordTop[1];
	vtxCoord[2] = prism.vtxCoordTop[2];
	vtxCoord[3] = prism.vtxCoordBottom[0];
	vtxCoord[4] = prism.vtxCoordBottom[1];
	vtxCoord[5] = prism.vtxCoordBottom[2];

	for (int idxface = 0; idxface < nFaces; idxface++) {	// 8 faces
		vec3 v0 = vtxCoord[FACEID_TO_CORNERID[idxFace * 3]];
		vec3 v1 = vtxCoord[FACEID_TO_CORNERID[idxFace * 3 + 1]];
		vec3 v2 = vtxCoord[FACEID_TO_CORNERID[idxFace * 3 + 2]];

		float t = 0.0;
		vec3 rayO = vec3(r.o);
		vec3 rayD = vec3(r.d);
		bool bhit = rayIntersectsTriangle(rayO, rayD, v0, v1, v2, t);

		if (bhit) {
			nHit++;
			
			if (min_t > t) {
				min_t = t;
				tInRec.t = float(t);
				tInRec.hitFaceid = idxface; 
			}
			if (max_t < t) {
				max_t = t;
				tOutRec.t = float(t);
				tOutRec.t = idxface;
				if (idxface == 1) {
					tOutRec.nextlayerId = prism.m_iLayer + 1;	// the next prism to be traversed is in the lower layer 
				}
				else if (idxface == 0) {
					tOutRec.nextlayerId = prism.m_iLayer - 1;	// the next prism to be traversed is in the upper layer
				} 
				else {
					tOutRec.nextlayerId = prism.m_iLayer;
				}
			}
		}
	}

	if (nHit == 2) {
		nextCellId = getAdjacentCellId(prism, tOutRec.hitFaceid);
	}
	else {	// specical case when ray hit on the edge 
		nextCellId = -1;
	}

	return nHit;
}

void GetAs0(inout MPASPrism prism, inout float A[12])
{
	float x1 = (prism.vtxCoordTop[0].x);
	float y1 = (prism.vtxCoordTop[0].y);
	float z1 = (prism.vtxCoordTop[0].z);

	float x2 = (prism.vtxCoordTop[1].x);
	float y2 = (prism.vtxCoordTop[1].y);
	float z2 = (prism.vtxCoordTop[1].z);

	float x3 = (prism.vtxCoordTop[2].x);
	float y3 = (prism.vtxCoordTop[2].y);
	float z3 = (prism.vtxCoordTop[2].z);

	A[1] = (+x3 * y1 - x1 * y3);
	A[2] = (+x3 * z1 - x1 * z3);
	A[3] = (+y3 * z1 - y1 * z3);
	A[4] = (-x2 * y1 + x1 * y2);
	A[5] = (-x2 * z1 + x1 * z2);
	A[6] = (-y2 * z1 + y1 * z2);
	A[7] = (-x2 * y1 + x3 * y1 + x1 * y2 - x3 * y2 - x1 * y3 + x2 * y3);//x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1 + y3);//(-x2*y1 + x3*y1 + x1*y2 - x3*y2 - x1*y3 + x2*y3);
	A[8] = (-x2 * z1 + x3 * z1 + x1 * z2 - x3 * z2 - x1 * z3 + x2 * z3);//x3 * (z1 - z2) + x1 * (z2 - z3) + x2 * (-z1 + z3);//(-x2*z1 + x3*z1 + x1*z2 - x3*z2 - x1*z3 + x2*z3);
	A[9] = (-y2 * z1 + y3 * z1 + y1 * z2 - y3 * z2 - y1 * z3 + y2 * z3);//y3 * (z1 - z2) + y1 * (z2 - z3) + y2 * (-z1 + z3);//(-y2*z1 + y3*z1 + y1*z2 - y3*z2 - y1*z3 + y2*z3);
	A[10] = (x3*y2*z1);
	A[11] = (-x2 * y3*z1 - x3 * y1*z2 + x1 * y3*z2 + x2 * y1*z3 - x1 * y2*z3);
}

void getScalarValue(inout MPASPrism prism, inout vec3 scalars[2]) {
	for (int iFace = 0; iFace < 2; iFace++) {
		int layerId = prism.m_iLayer + iFace;
		scalars[iFace].x = texelFetch(CLIMATE_VALS_VAR, prism.idxVtx[0] * nLayers + layerId).r; //CLIMATE_VALS_VAR[idxVtx[0] * nLayers + (m_iLayer + iFace)];
		scalars[iFace].y = texelFetch(CLIMATE_VALS_VAR, prism.idxVtx[1] * nLayers + layerId).r; //CLIMATE_VALS_VAR[idxVtx[0] * nLayers + (m_iLayer + iFace)];
		scalars[iFace].z = texelFetch(CLIMATE_VALS_VAR, prism.idxVtx[2] * nLayers + layerId).r; //CLIMATE_VALS_VAR[idxVtx[0] * nLayers + (m_iLayer + iFace)];
	}
}

void GetUV(in vec3 O, in vec3 Q, inout float A[12],
	inout float u, inout float v) {
	vec3 QO = (Q - O);//*Factor;
	float denominator = (A[9] * QO.x - A[8] * QO.y + A[7] * QO.z);
	u = (A[3] * QO.x - A[2] * QO.y + A[1] * QO.z) / denominator;
	v = (A[6] * QO.x - A[5] * QO.y + A[4] * QO.z) / denominator;
}

float GetInterpolationValue2(inout MPASPrism prism, in const float u, in const float v,
	in const vec3 Q, inout vec3 fT, inout vec3 fB) {
	vec3 baryCoord = vec3(1.0 - u - v, u, v);
	vec3 m1 = baryCoord.x * prism.vtxCoordTop[0] + baryCoord.y * prism.vtxCoordTop[1] + baryCoord.z * prism.vtxCoordTop[2];
	vec3 m2 = baryCoord.x * prism.vtxCoordBottom[0] + baryCoord.y * prism.vtxCoordBottom[1] + baryCoord.z * prism.vtxCoordBottom[2];

	float scalar_m1 = dot(baryCoord, fT);
	float scalar_m2 = dot(baryCoord, fB);
	float t3 = length(Q - m2) / length(m1 - m2);
	float lerpedVal = mix(scalar_m2, scalar_m1, t3);	//lerp()

	return lerpedVal;
}

vec3 GetNormal(const vec3 position, const float A[12],
	const float B[8], const float OP1, const float OP4) {
	float delx = 0.0f, dely = 0.0f, delz = 0.0f;//partial derivative of scalar value at sampling point.
	float inv_denom = 1.0f / (A[9] * position.x - A[8] * position.y + A[7] * position.z);
	float C0 = B[0] * OP1*(A[9] * position.x - A[8] * position.y + A[7] * position.z)*(B[1] * position.x + B[3] * position.y + B[2] * position.z);
	//B10 OP1 (A9 qx			- A8 qy				+ A7 qz)		 (B11 qx			+ B13 qy		 + B12 qz)
	float temp = (A[10] + A[11])*OP4 + OP1 * (A[9] * position.x - A[8] * position.y + A[7] * position.z);
	//(A10 + A11) OP4	 + OP1 (A9 qx			- A8 qy				+ A7 qz)
	float C1 = B[6] - B[0] * (B[4] - B[6])*temp;
	//B16 - B10 (B14 - B16)
	float C2 = B[7] - B[0] * (B[5] - B[7])*temp;
	//B17 - B10 (B15 - B17)

	delx = (inv_denom*inv_denom)*
		(A[9] * C0 +
			C1 * ((A[3] * A[8] - A[2] * A[9])*position.y + (-A[3] * A[7] + A[1] * A[9])*position.z) +
			C2 * ((A[6] * A[8] - A[5] * A[9])*position.y + (-A[6] * A[7] + A[4] * A[9])*position.z)
			);
	//A9 C0 +
	//		C1 ((A3 A8		- A2 A9) qy				+ (-A3 A7		+ A1 A9) qz) +
	//		C2 ((A6 A8		- A5 A9) qy				+ (-A6 A7		+ A4 A9) qz)
	dely = (inv_denom*inv_denom)*
		(-A[8] * C0 +
			C1 * ((-A[3] * A[8] + A[2] * A[9])*position.x + (A[2] * A[7] - A[1] * A[8])*position.z) +
			C2 * ((-A[6] * A[8] + A[5] * A[9])*position.x + (A[5] * A[7] - A[4] * A[8])*position.z)
			);
	//-A8 C0 +
	//		C1 ((-A3 A8			+ A2 A9) qx			+ (A2 A7		- A1 A8) qz) +
	//		C2 ((-A6 A8			+ A5 A9) qx			+ (A5 A7		- A4 A8) qz)
	delz = (inv_denom*inv_denom)*
		(A[7] * C0 +
			C1 * ((A[3] * A[7] - A[1] * A[9])*position.x + (-A[2] * A[7] + A[1] * A[8])*position.y) +
			C2 * ((A[6] * A[7] - A[4] * A[9])*position.x + (-A[5] * A[7] + A[4] * A[8])*position.y)
			);
	//A7 C0 +
	// C1 ((A3 A7		- A1 A9) qx				+ (-A2 A7		+ A1 A8) qy) +
	// C2 ((A6 A7		- A4 A9) qx				+ (-A5 A7		+ A4 A8) qy)
	return normalize(vec3(delx, dely, delz));
}

void ComputeVerticesNormalTop(inout MPASPrism prism, inout vec3 vtxNormals[3]) {	
	// for each one of the three corners of top face of current prism
	// compute their shared normals respectively
	for (iCorners = 0; iCorners < 3; iCorners++) {
		int idxCorner = prism.idxVtx[iCorner];	//TRIANGLE_TO_CORNERS_VAR[m_prismtId*3+iCorner]
		// for each nNeighbors prisms (including current prism and nNeighbors-1 neighbor prisms)
		int nNeighbors = texelFetch(CORNER_TO_TRIANGLES_DIMSIZES, idxCorner);
		vec3 avgNormals = vec3(0.0f);
		for (int iPrism = 0; iPrism < nNeighbors; iPrism++) {	//weighted normal 
			int id = texelFetch(CORNERS_TO_TRIANGLE_VAR, idxCorner * maxEdges + iPrism).x;	// CORNER_TO_TRIANGLES_VAR[idxCorner * nNeighbors + iPrism]

			MPASPrism curPrismHitted;	// (id, m_iLayer)
			ReloadVtxInfo(id, prism.m_iLayer, curPrismHitted);
			// find the vertex in curPrism whose global index == idxCorner
			int i = 0;
			for (; i < 3; i++) {
				if (curPrismHitted.idxVtx[i] == idxCorner)
					break;
			}
			float A[12];
			GetAs0(curPrismHitted, A);
			vec3 fTB[2];
			GetScalarValue(curPrismHitted, fTB);
			float OP1, OP4, B[8];
			GetBs(curPrismHitted, A, fTB, B, OP1, OP4);
			avgNormal += GetNormal(curPrismHitted.vtxCoordTop[i], A, B, OP1, OP4);

			if (prism.m_iLayer > 0) {
				MPASPrism curPrismHitted;	// (id, m_iLayer)
				ReloadVtxInfo(id, prism.m_iLayer - 1, curPrismHitted); 
				float A[12];
				GetAs0(curPrismHitted, A);
				vec3 fTB[2];
				GetScalarValue(curPrismHitted, fTB);
				float OP1, OP4, B[8];
				GetBs(curPrismHitted, A, fTB, B, OP1, OP4);
				avgNormal += GetNormal(curPrismHitted.vtxCoordBottom[i], A, B, OP1, OP4);
			}
		}
		vtxNormals[iCorner] = normalize(avgNormal);
	}
}

void ComputeVerticesNormalBottom(inout MPASPrism prism, inout vec3 vtxNormals[3]) {  
	// for each one of the three corners of top face of current prism
	// compute their shared normals respectively
	for (iCorners = 0; iCorners < 3; iCorners++) {
		int idxCorner = prism.idxVtx[iCorner];	//TRIANGLE_TO_CORNERS_VAR[m_prismtId*3+iCorner]
		// for each nNeighbors prisms (including current prism and nNeighbors-1 neighbor prisms)
		int nNeighbors = texelFetch(CORNER_TO_TRIANGLES_DIMSIZES, idxCorner);
		vec3 avgNormals = vec3(0.0f);
		for (int iPrism = 0; iPrism < nNeighbors; iPrism++) {	//weighted normal 
			int id = texelFetch(CORNERS_TO_TRIANGLE_VAR, idxCorner * maxEdges + iPrism).x;	// CORNER_TO_TRIANGLES_VAR[idxCorner * nNeighbors + iPrism]
			MPASPrism curPrismHitted;	// (id, m_iLayer)
			ReloadVtxInfo(id, prism.m_iLayer, curPrismHitted);
			// find the vertex in curPrism whose global index == idxCorner
			int i = 0;
			for (; i < 3; i++) {
				if (curPrismHitted.idxVtx[i] == idxCorner)
					break;
			}
			float A[12];
			GetAs0(curPrismHitted, A);
			vec3 fTB[2];
			GetScalarValue(curPrismHitted, fTB);
			float OP1, OP4, B[8];
			GetBs(curPrismHitted, A, fTB, B, OP1, OP4);
			avgNormal += GetNormal(curPrismHitted.vtxCoordBottom[i], A, B, OP1, OP4);
		}
		vtxNormals[iCorner] = normalize(avgNormal);
	}
}

void main(){
	// perform ray-casting through triangle mesh
	vec3 o_eye = (inv_mvm * vec4(0, 0, 0, 1.0)).xyz;
	Ray ray;
	ray.o = o_eye;
	ray.d = normalize(g2f.o_pos - o_eye);
	int triangle_id = g2f.triangle_id;

	HitRec tInHitRecord, tOutHitRecord, tmpInRec, tmpOutRec;
	tInHitRecord.hitFaceid = -1;
	tInHitRecord.t = FLOAT_MAX;
	tInHitRecord.nextlayerId = -1;

	tOutHitRecord.hitFaceid = -1;
	tOutHitRecord.t = FLOAT_MIN;
	tOutHitRecord.nextlayerId = -1;
	
	tmpInRec.t = FLOAT_MAX;
	tmpOutRec.t = FLOAT_MIN;
	tmpOutRec.nextlayerId = -1;

	int nHit = 0;
	int nextCell = -1;
	int tmpNextCell = -1;
	uint curPrismHittedId = triangle_id;

	MPASPrism curPrismHitted;
	curPrismHitted.idxVtx[0] = -1;
	curPrismHitted.idxVtx[1] = -1;
	curPrismHitted.idxVtx[2] = -1;
	ReloadVtxInfo(int(curPrismHittedId), 0, curPrismHitted);

	int tmpNhit = rayPrismIntersection(curPrismHitted, ray, tmpInRec, tmpOutRec, tmpNextCellId);

	if (tmpNHit > 0) {
		nHit = tmpNHit;
		nextCellId = tmpNextCellId;
		curPrismHittedId = triangle_id;
		tInHitRecord = (tInHitRecord.t > tmpInRec.t) ? tmpInRec : tInHitRecord;
		tOutHitRecord = (tOutHitRecord.t < tmpOutRec.t) ? tmpOutRec : tOutHitRecord;
	}

	// compute normal on three vertices of the top face of current prism
	vec3 vtxNormalsTop[3];
	ComputeVerticesNormalTop(curPrismHitted, vtxNormals);

	do {	// loop through 3d grid
		if (tInHitRecord.t < 0.0f)
			tInHitRecord.t = 0.0f;

		float t = tInHitRecord.t;
		float A[12];
		GetAs0(curPrismHitted, A);
		vec3 fTB[2];
		GetScalarValue(curPrismHitted, fTB);

		position = ray.o + ray.d * t;
		float u, v;
		GetUV(vec3(0.0f), position, A, u, v);
		float scalar_last = GetInterpolationValue2(curPrismHitted, u, v, position, fTB[0], fTB[1]);
		vec3 normals_last = 
		t += stepsize;
		
		float scalar;
		while (t <= tOutHitRecord.t) {
			position = ray.o + ray.d * t;
			float u, v;
			GetUV(vec3(0.0f), position, A, u, v);
			scalar = GetInterpolationValue2(curPrismHitted, u, v, position, fTB[0], fTB[1]);
			
			if ((scalar_last - threshold) * (scalar - threshold) < 0) {
				float offset = (scalar - threshold) / (scalar - scalar_last);
			}

			scalar_last = scalar;
			t += stepsize;
		}

		t = tOutHitRecord.t; 
		position = ray.o + ray.d * t;
		GetUV(vec3(0.0f), position, A, u, v);
		scalar = GetInterpolationValue2(curPrismHitted, u, v, position, fTB[0], fTB[1]);
		if ((scalar_last - threshold) * (scalar - threshold) < 0) {
			float offset = (scalar - threshold) / (scalar - scalar_last);
		}
		
		// intersect ray with next cell 
		if (nextCellId > -1) {	// more prism to intersect
			curPrismHittedId = nextCellId;
			tInHitRecord = tOutHitRecord;
			if (tInHitRecord.nextlayerId <= -1){	// short from inner part to the outside  	
				break;
			}
			ReloadVtxInfo(int(curPrismHittedId), tInHitRecord.nextlayerId, curPrismHitted);	// update curPrismHitted to nextCell

			nHit = rayPrismIntersection(curPrismHitted, ray, tInHitRecord, tOutHitRecord, nextCellId); 
		}
		else {
			break;
		}
	} while (nHit > 0 && tOutHitRecord.nextlayerId != -1)
} 

