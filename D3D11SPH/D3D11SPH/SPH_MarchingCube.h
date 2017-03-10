#ifndef SPH_MARCHING_CUBE_H_
#define SPH_MARCHING_CUBE_H_

#include "Control.h"
#include "SPH_System.h" 
#include<exception>
using namespace std;
//1.Figure out density of every vertex
//2.Compare with the threshold value p =  1050 this is only a try
//

class Grid;
class SPHSystem;

class MarchingCube{

	SPHBox m_spacebox;

	int m_iXAxisGridNumber;
	int m_iYAxisGridNumber;
	int m_iZAxisGridNumber;

	float m_gridlength;
	int m_gridscount;

	float m_thresholdvalue;
	float m_CubeLength;
	float m_CubeRange;

	vector<SPHCube> m_Cubes; //每个单元格对应的顶点

	bool isRender;

public:

	vector<int> m_testindexList;

	vector<SPHTriangle> m_trianglePos;

	void Polygonise(Grid & grid,float isovalue,SPHSystem & sphSystem);

	void GridValueFigure(Grid & grid,SPHSystem & sphSystem);

	XMFLOAT3 Interpolation(float isovalue,XMFLOAT3 p1,XMFLOAT3 p2,float v1,float v2);

	float Thresholdvalue() const { return m_thresholdvalue; }
	void Thresholdvalue(float val) { m_thresholdvalue = val; }

	float CubeLength() const { return m_CubeLength; }
	void CubeLength(float val);

	float CubeRange() const { return m_CubeRange; }
	void CubeRange(float val);

	bool IsRender() const { return isRender; }
	void IsRender(bool val) { isRender = val; }

	void Reset()
	{
		if(!m_testindexList.empty())
			m_testindexList.clear();

		if(!m_trianglePos.empty())
			m_trianglePos.clear();

	};


	void     SetMarchingCube(SPHBox limitbox, float length);
	void     Init();
	void     returnToDefault();

	XMFLOAT3 figureCellPos(int index);

	void     FigureCellVertexPos(int index,XMFLOAT3 * VertexPos);

	void     InitAllVertexPos();

	int Gridscount() const { return m_gridscount; }
};



#endif