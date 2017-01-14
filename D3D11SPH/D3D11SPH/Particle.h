#include<windows.h>
#include<xnamath.h>
#ifndef PARTICLE_H_
#define PARTICLE_H_

class Particle
{
private:

public:

	XMFLOAT3 m_vPosRender;
	XMFLOAT3 m_vPosRealScale;
	XMFLOAT3 m_vAcceleration;
	XMFLOAT3 m_vVelocity;
	float m_fDensity;
	float m_fPessure;
	float m_fViscosity;


	Particle();
	~Particle(){};
};






#endif