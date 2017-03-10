#include"Grid.h"

void Grid::InitGrid(SPHBox spacebox,float particleradius)
{
	m_spacebox = spacebox;

	m_iXAxisGridNumber = (spacebox.max_x - spacebox.min_x) / ( 2 * particleradius);
	m_iYAxisGridNumber = (spacebox.max_y - spacebox.min_y) / ( 2 * particleradius);
	m_iZAxisGridNumber = (spacebox.max_z - spacebox.min_z) / ( 2 * particleradius);

	m_gridscount =  m_iXAxisGridNumber * m_iYAxisGridNumber * m_iZAxisGridNumber;

	m_gridlength = 2 * particleradius; // 边长为粒子核半径 * 2


	//InitAllVertexPos();

	//m_girdunit.resize(m_gridscount);


}

void Grid::InsertParticles(SPHSystem & sphsystem)
{
	//更新grid表格中的信息时先清空
	for(int i = 0; i < m_gridscount; i++ )
	{
		if(!m_girdunit[i].empty())
			m_girdunit[i].clear();
	}

// 	if(!m_girdunit.empty())
// 		m_girdunit.clear();

	bool isexception = false;
	for(int i = 0; i < sphsystem.m_iParticleCount; i++)
	{
		float pos_x = sphsystem.m_ParticlePool[i].m_vPosRealScale.x;
		float pos_y = sphsystem.m_ParticlePool[i].m_vPosRealScale.y;
		float pos_z = sphsystem.m_ParticlePool[i].m_vPosRealScale.z;

		int index = findCell(pos_x,pos_y,pos_z);

		//if(-1 == index)
		//	continue;

		if(index < 0 || index >= m_gridscount)
		{
			MessageBoxA(0,"Particle exceed the border. Return to Default state automatically",0,0);
			isexception = true;
			break;
		}
		m_girdunit[index].push_back(i); // 将粒子插入grid对应的链表

	}

	if(isexception)
	{
		SPH_UIControl::ReturnToDefault();
	}


}

int Grid::findCell(XMFLOAT3 & vPos)
{
    float x = vPos.x;
	float y = vPos.y;
	float z = vPos.z;
	 
	if(x < m_spacebox.min_x || x > m_spacebox.max_x || y < m_spacebox.min_y || y > m_spacebox.max_y || z < m_spacebox.min_z || z > m_spacebox.max_z)
	{
	//	MessageBoxA(0,"Over Grid",0,0);
		//return -1;
	}


	int index_x,index_y,index_z;

	index_x = (x - m_spacebox.min_x) / m_gridlength;
	index_y = (y - m_spacebox.min_y) / m_gridlength;
	index_z = (z - m_spacebox.min_z) / m_gridlength;

	int index = index_z * m_iYAxisGridNumber * m_iXAxisGridNumber + index_y * m_iXAxisGridNumber + index_x;

	return index;

}

int Grid::findCell(float x, float y, float z)
{
	if(x < m_spacebox.min_x || x > m_spacebox.max_x || y < m_spacebox.min_y || y > m_spacebox.max_y || z < m_spacebox.min_z || z > m_spacebox.max_z)
	{
	//	MessageBoxA(0,"Over Grid",0,0);
	//	return -1;
	}

	int index_x,index_y,index_z;

	index_x = (x - m_spacebox.min_x) / m_gridlength;
	index_y = (y - m_spacebox.min_y) / m_gridlength;
	index_z = (z - m_spacebox.min_z) / m_gridlength;

	int index = index_z * m_iYAxisGridNumber * m_iXAxisGridNumber + index_y * m_iXAxisGridNumber + index_x;

	return index;
}

void Grid::findNeighborCells(XMFLOAT3 & vPos,int * neighborcells)
{
	float x = vPos.x;
	float y = vPos.y;
	float z = vPos.z;

	int index = findCell(vPos);

	XMFLOAT3 neighborPos[8];
	neighborPos[0] = XMFLOAT3(x - m_gridlength/2.0f,y - m_gridlength/2.0f,z - m_gridlength/2.0f);
	neighborPos[1] = XMFLOAT3(x - m_gridlength/2.0f,y - m_gridlength/2.0f,z + m_gridlength/2.0f);
	neighborPos[2] = XMFLOAT3(x - m_gridlength/2.0f,y + m_gridlength/2.0f,z - m_gridlength/2.0f);
	neighborPos[3] = XMFLOAT3(x - m_gridlength/2.0f,y + m_gridlength/2.0f,z + m_gridlength/2.0f);
	neighborPos[4] = XMFLOAT3(x + m_gridlength/2.0f,y - m_gridlength/2.0f,z - m_gridlength/2.0f);
	neighborPos[5] = XMFLOAT3(x + m_gridlength/2.0f,y - m_gridlength/2.0f,z + m_gridlength/2.0f);
	neighborPos[6] = XMFLOAT3(x + m_gridlength/2.0f,y + m_gridlength/2.0f,z - m_gridlength/2.0f);
	neighborPos[7] = XMFLOAT3(x + m_gridlength/2.0f,y + m_gridlength/2.0f,z + m_gridlength/2.0f);


	neighborcells[0] = index;
	for(int i = 1;i < 9; i++)
	{
		neighborcells[i] = findCell(neighborPos[i - 1]);

		// 有一定概率重复
		for(int j = 0; j < i; j++)
			if(neighborcells[i] == neighborcells[j])
				neighborcells[i] = -1;
	}

}

Grid::~Grid()
{

}

/*XMFLOAT3 Grid::figureCellPos(int index)
{
	int index_x,index_y,index_z;

	index_z = index / (m_iYAxisGridNumber * m_iXAxisGridNumber) ;
	index_y = (index - index_z * (m_iYAxisGridNumber * m_iXAxisGridNumber)) / m_iXAxisGridNumber ;
	index_x = index - index_z * (m_iYAxisGridNumber * m_iXAxisGridNumber) - index_y * m_iXAxisGridNumber;

	float x = index_x * m_gridlength + m_spacebox.min_x + m_gridlength/2.0f;
	float y = index_y * m_gridlength + m_spacebox.min_y + m_gridlength/2.0f;
	float z = index_z * m_gridlength + m_spacebox.min_z + m_gridlength/2.0f;

	return XMFLOAT3(x,y,z); //返回方格的中心位置



}

void Grid::FigureCellVertexPos(int index,XMFLOAT3 * VertexPos)
{

	XMFLOAT3 vCenter = figureCellPos(index);

	VertexPos[0] = XMFLOAT3(vCenter.x - m_gridlength/2.0f,vCenter.y - m_gridlength/2.0f,vCenter.z + m_gridlength/2.0f);
	VertexPos[1] = XMFLOAT3(vCenter.x + m_gridlength/2.0f,vCenter.y - m_gridlength/2.0f,vCenter.z + m_gridlength/2.0f);
	VertexPos[2] = XMFLOAT3(vCenter.x + m_gridlength/2.0f,vCenter.y - m_gridlength/2.0f,vCenter.z - m_gridlength/2.0f);
	VertexPos[3] = XMFLOAT3(vCenter.x - m_gridlength/2.0f,vCenter.y - m_gridlength/2.0f,vCenter.z - m_gridlength/2.0f);
	VertexPos[4] = XMFLOAT3(vCenter.x - m_gridlength/2.0f,vCenter.y + m_gridlength/2.0f,vCenter.z + m_gridlength/2.0f);
	VertexPos[5] = XMFLOAT3(vCenter.x + m_gridlength/2.0f,vCenter.y + m_gridlength/2.0f,vCenter.z + m_gridlength/2.0f);
	VertexPos[6] = XMFLOAT3(vCenter.x + m_gridlength/2.0f,vCenter.y + m_gridlength/2.0f,vCenter.z - m_gridlength/2.0f);
	VertexPos[7] = XMFLOAT3(vCenter.x - m_gridlength/2.0f,vCenter.y + m_gridlength/2.0f,vCenter.z - m_gridlength/2.0f);
}

void Grid::InitAllVertexPos()
{
	for(int i = 0;i < 4096; i++)
	{
		m_vertex[i].resize(8);

		XMFLOAT3 VertexPos[8];

		FigureCellVertexPos(i,VertexPos);

		for(int j = 0;j < 8;j++)

		m_vertex[i][j].pos = VertexPos[j];

	}
}*/


