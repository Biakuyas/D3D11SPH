#ifndef SPH_SYSTEM_H_
#define SPH_SYSTEM_H_
#include"Particle.h"
#include"VertexIndexControl.h"

const int Particle_number = 8000;
const float PI = 3.141592f;
const float WallThickness = 0.008f;
const float WallLength = 8;

struct Box
{
	float min_x,min_y,min_z;
	float max_x,max_y,max_z;

	Box(float minx,float miny,float minz,float maxx,float maxy, float maxz)
	{
		min_x = minx;
		min_y = miny;
		min_z = minz;
		max_x = maxx;
		max_y = maxy;
		max_z = maxz;
	}
	Box()
	{
		min_x = 0;
		min_y = 0;
		min_z = 0;
		max_x = 0;
		max_y = 0;
		max_z = 0;
	}

};

class SPHSystem
{
private:
	float m_unitScale;
	float m_viscosity;
	float m_restDensity; 
	float m_pointMass;
	float m_smoothRadius;
	float m_gasConstantK;
	XMFLOAT3 m_gravity;
	Box m_limitBox;

	float m_boundaryCoefficient;
	float m_boundaryDampening ;

	float m_kernelPoly6;
	float m_kernelSpiky;
	float m_kernelViscosity;
	float Distance(XMFLOAT3 a, XMFLOAT3 b){ return pow((a.x - b.x) *(a.x - b.x) + (a.y - b.y)*(a.y - b.y)+ (a.z - b.z)*(a.z - b.z),0.5f);}
	float DistanceSquare(XMFLOAT3 a, XMFLOAT3 b){ return (a.x - b.x) *(a.x - b.x) + (a.y - b.y)*(a.y - b.y)+ (a.z - b.z)*(a.z - b.z);}

public:
	Particle m_ParticlePool[Particle_number];
	SPHSystem();
	~SPHSystem(){};

	int m_iParticleCount;
	void SPH_DensityPressure();
	void SPH_AccelerationFigure();
	void SPH_ParticleUpdate(const float deltatime);
	void SPH_Init();

	float DotF3(XMFLOAT3 a ,XMFLOAT3 b){return a.x * b.x + a.y * b.y + a.z * b.z;}
	XMFLOAT3 MinusF3(XMFLOAT3 a ,XMFLOAT3 b){return XMFLOAT3(a.x - b.x,a.y - b.y,a.z - b.z);}
	XMFLOAT3 PlusF3(XMFLOAT3 a ,XMFLOAT3 b){return XMFLOAT3(a.x + b.x,a.y + b.y,a.z + b.z);}
	XMFLOAT3 MultiplyFS(XMFLOAT3 a ,float fb){return XMFLOAT3(a.x * fb,a.y * fb,a.z * fb);}
	XMFLOAT3 DivideFS(XMFLOAT3 a ,float fb){ if(fb == 0) MessageBoxA(0,"DivideError",0,0);return XMFLOAT3(a.x / fb,a.y / fb,a.z / fb);}
};




#endif