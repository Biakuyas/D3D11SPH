#ifndef GRID_H_
#define GRID_H_
#include <vector>
#include"SPH_System.h"
#include"Control.h"
#include "SPH_MarchingCube.h"
#include "SPH_UI.h"

using namespace std;

class SPH_UIControl;
//互相包含，前置声明
class SPHSystem;

class Grid{
private:
	float m_gridlength;
	int m_gridscount;
	SPHBox m_spacebox;
	int m_iXAxisGridNumber;
	int m_iYAxisGridNumber;
	int m_iZAxisGridNumber;
	
public:

	vector<int> m_girdunit[5000];
	void InitGrid(SPHBox spacebox,float particleradius);
	void InsertParticles(SPHSystem & sphsystem);
	int findCell(XMFLOAT3 & vPos);
	int findCell(float x, float y, float z);
	void findNeighborCells(XMFLOAT3 & vPos,int * neighborcells);
	friend class MarchingCube;

	~Grid();

};

#endif