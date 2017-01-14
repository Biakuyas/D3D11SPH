#include"SPH_System.h"

SPHSystem::SPHSystem()
{
	m_unitScale         = 0.004f;			// 尺寸单位
	m_viscosity			= 1.0f;				// 粘度
	m_restDensity		= 1000.f;			// 密度
	m_pointMass			= 0.0004f;			// 粒子质量
	m_gasConstantK		= 1.0f;				// 理想气体方程常量
	m_smoothRadius		= 0.01f;			// 光滑核半径
//	m_iParticleCount      = 0;
	m_gravity = XMFLOAT3(0,-9.8f,0);

	//容器
	m_limitBox = Box(-WallLength * m_unitScale,0 * m_unitScale,-WallLength * m_unitScale,WallLength * m_unitScale,2 * WallLength * m_unitScale,WallLength * m_unitScale);

//	m_limitBox = Box(0,0,0,0 * m_unitScale,2 * 0 * m_unitScale,100 * m_unitScale);
	m_boundaryCoefficient = 10000.f;
	m_boundaryDampening = 256.f;

	//Poly6 Kernel
	m_kernelPoly6 = 315.0f/(64.0f * PI * pow(m_smoothRadius, 9));
	//Spiky Kernel
	m_kernelSpiky = -45.0f/(PI * pow(m_smoothRadius, 6));
	//Viscosity Kernel
	m_kernelViscosity = 45.0f/(PI * pow(m_smoothRadius, 6));

	


}

void SPHSystem::SPH_Init()
{
	m_iParticleCount      = 0;

	float pointDistance	= pow(m_pointMass/m_restDensity, 1.f/3.f); //粒子间距
//	_addFluidVolume(initFluidBox, pointDistance/m_unitScale);

	for(float i = m_limitBox.min_x; i <= m_limitBox.max_x; i += pointDistance)
		for(float j = m_limitBox.min_y; j <= m_limitBox.max_y; j += pointDistance)
			for(float p = m_limitBox.min_z; p <= m_limitBox.max_z; p += pointDistance)
			{
				
				m_ParticlePool[m_iParticleCount].m_vPosRender = XMFLOAT3(i / m_unitScale,j / m_unitScale,p / m_unitScale);
				m_ParticlePool[m_iParticleCount].m_vPosRealScale = XMFLOAT3(i,j + 10 *m_unitScale ,p);
				m_ParticlePool[m_iParticleCount].m_fViscosity = 0;
				m_ParticlePool[m_iParticleCount].m_fDensity = 0;
				m_ParticlePool[m_iParticleCount].m_fPessure = 0;
				m_ParticlePool[m_iParticleCount].m_vAcceleration = XMFLOAT3(0,0,0);
				m_ParticlePool[m_iParticleCount].m_vVelocity = XMFLOAT3(0,0,0);
				//if(m_iParticleCount > 10000)
				//	MessageBoxA(0,"Over Pool",0,0);
				m_iParticleCount++;
			}

	/*for(int i = 0; i < Particle_number / 100; i++)
		for(int j = 0; j < 10; j++)
			for(int p = 0; p < 10; p++)
			{
				m_ParticlePool[i * 100 + j * 10 + p].m_vPosRender = XMFLOAT3(p * 2, i * 2 + 10, j * 2);
				m_ParticlePool[i * 100 + j * 10 + p].m_vPosRealScale = XMFLOAT3(p * 2 * m_unitScale, (i * 2 + 10) * m_unitScale, j * 2 * m_unitScale);
				m_ParticlePool[i * 100 + j * 10 + p].m_fViscosity = 0;
				m_ParticlePool[i * 100 + j * 10 + p].m_fDensity = 0;
				m_ParticlePool[i * 100 + j * 10 + p].m_fPessure = 0;
				m_ParticlePool[i * 100 + j * 10 + p].m_vAcceleration = XMFLOAT3(0,0,0);
				m_ParticlePool[i * 100 + j * 10 + p].m_vVelocity = XMFLOAT3(0,0,0);

			}*/


}

void SPHSystem::SPH_DensityPressure()
{
	//这次遍历的时候可以记录一下周围的点，下次就不用遍历来找了，不过我估计会用更加优化的方法来做
	/*for(int i = 0;i < Particle_number; i++)
	{
		float sum = 0;

		for(int j = 0; j < Particle_number; j++)
		{
			float h2 = m_smoothRadius * m_smoothRadius;
			float r2 = DistanceSquare(m_ParticlePool[i].m_vPosRealScale,m_ParticlePool[j].m_vPosRealScale);

			if(r2 > h2)
				continue;

			sum += pow((h2 - r2),3);
	
		}
		m_ParticlePool[i].m_fDensity = m_kernelPoly6 * sum * m_pointMass; //因为所有粒子的重量都一样所以这里可以提出来乘
		m_ParticlePool[i].m_fPessure = (m_ParticlePool[i].m_fDensity - m_restDensity) * m_gasConstantK;
	}*/
	for(int i = 0;i < m_iParticleCount; i++)
	{
		float sum = 0;

		for(int j = 0; j < m_iParticleCount; j++)
		{
			float h2 = m_smoothRadius * m_smoothRadius;
			float r2 = DistanceSquare(m_ParticlePool[i].m_vPosRealScale,m_ParticlePool[j].m_vPosRealScale);

			if(r2 > h2)
				continue;

			sum += pow((h2 - r2),3);

		}
		m_ParticlePool[i].m_fDensity = m_kernelPoly6 * sum * m_pointMass; //因为所有粒子的重量都一样所以这里可以提出来乘
		m_ParticlePool[i].m_fPessure = (m_ParticlePool[i].m_fDensity - m_restDensity) * m_gasConstantK;
	}


}

void SPHSystem::SPH_AccelerationFigure()
{
/*
	for(int i = 0;i < Particle_number; i++)
	{
		XMFLOAT3 accel_sum = XMFLOAT3(0,0,0);

		for(int j = 0; j < Particle_number; j++)
		{
			float h2 = m_smoothRadius * m_smoothRadius;
			float r2 = DistanceSquare(m_ParticlePool[i].m_vPosRealScale,m_ParticlePool[j].m_vPosRealScale);

			if(r2 > h2)
				continue;

			float r = Distance(m_ParticlePool[i].m_vPosRealScale,m_ParticlePool[j].m_vPosRealScale);


			float h2_r2 = h2 - r2;
			float h_r = m_smoothRadius - r;
			XMFLOAT3 ri_rj = MinusF3(m_ParticlePool[i].m_vPosRealScale,m_ParticlePool[j].m_vPosRealScale);

			float pscalar = -m_pointMass * m_kernelSpiky * h_r * h_r * (m_ParticlePool[i].m_fPessure + m_ParticlePool[j].m_fPessure) / (2.0f * m_ParticlePool[i].m_fDensity * m_ParticlePool[j].m_fDensity);
			
			if(r != 0)
			accel_sum = PlusF3(accel_sum,DivideFS(MultiplyFS(ri_rj,pscalar),r) ); //(ri_rj * pscalar / r )


			float vscalar = m_kernelViscosity * m_viscosity * h_r * m_pointMass / (m_ParticlePool[i].m_fDensity * m_ParticlePool[j].m_fDensity);
			accel_sum = PlusF3(accel_sum,MultiplyFS(MinusF3(m_ParticlePool[j].m_vVelocity,m_ParticlePool[i].m_vVelocity),vscalar) );  //(vj - vi * vscalar)

		}

		m_ParticlePool[i].m_vAcceleration = accel_sum;
	}
*/

	for(int i = 0;i < m_iParticleCount; i++)
	{
		XMFLOAT3 accel_sum = XMFLOAT3(0,0,0);

		for(int j = 0; j < m_iParticleCount; j++)
		{
			float h2 = m_smoothRadius * m_smoothRadius;
			float r2 = DistanceSquare(m_ParticlePool[i].m_vPosRealScale,m_ParticlePool[j].m_vPosRealScale);

			if(r2 > h2)
				continue;

			float r = Distance(m_ParticlePool[i].m_vPosRealScale,m_ParticlePool[j].m_vPosRealScale);


			float h2_r2 = h2 - r2;
			float h_r = m_smoothRadius - r;
			XMFLOAT3 ri_rj = MinusF3(m_ParticlePool[i].m_vPosRealScale,m_ParticlePool[j].m_vPosRealScale);

			float pscalar = -m_pointMass * m_kernelSpiky * h_r * h_r * (m_ParticlePool[i].m_fPessure + m_ParticlePool[j].m_fPessure) / (2.0f * m_ParticlePool[i].m_fDensity * m_ParticlePool[j].m_fDensity);

			if(r != 0)
				accel_sum = PlusF3(accel_sum,DivideFS(MultiplyFS(ri_rj,pscalar),r) ); //(ri_rj * pscalar / r )


			float vscalar = m_kernelViscosity * m_viscosity * h_r * m_pointMass / (m_ParticlePool[i].m_fDensity * m_ParticlePool[j].m_fDensity);
			accel_sum = PlusF3(accel_sum,MultiplyFS(MinusF3(m_ParticlePool[j].m_vVelocity,m_ParticlePool[i].m_vVelocity),vscalar) );  //(vj - vi * vscalar)

		}

		m_ParticlePool[i].m_vAcceleration = accel_sum;
	}
}

void SPHSystem::SPH_ParticleUpdate(const float deltatime)
{
	for(int i = 0;i < m_iParticleCount;i++)
	{

	//	float testdeltatime = deltatime * 0.05f;//在1/1000时间状态下进行模拟，不然有点看不清
		float testdeltatime = 0.003f;
		//加上重力和弹力更新加速度
		m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,m_gravity);

		//一个上口不封的盒子

		//X-Wall min 
		float fbounceParameter = WallThickness - (m_ParticlePool[i].m_vPosRealScale.x - m_limitBox.min_x);
		//到达墙壁范围内
		if(fbounceParameter > 0) 
		{
			//所受到的墙的弹力减去阻力
			XMFLOAT3 vec(1,0,0);
			float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

			XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

			m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

		}

		//X-Wall max 
		fbounceParameter = WallThickness - (m_limitBox.max_x - m_ParticlePool[i].m_vPosRealScale.x);
		//到达墙壁范围内
		if(fbounceParameter > 0) 
		{
			//所受到的墙的弹力减去阻力
			XMFLOAT3 vec(-1,0,0);
			float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

			XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

			m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

		}

		//Y-Wall min 
		fbounceParameter = WallThickness - (m_ParticlePool[i].m_vPosRealScale.y - m_limitBox.min_y);
		//到达墙壁范围内
		if(fbounceParameter > 0) 
		{
			//所受到的墙的弹力减去阻力
			XMFLOAT3 vec(0,1,0);
			float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

			XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

			m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

		}

		//Z-Wall min 
		fbounceParameter = WallThickness - (m_ParticlePool[i].m_vPosRealScale.z - m_limitBox.min_z);
		//到达墙壁范围内
		if(fbounceParameter > 0) 
		{
			//所受到的墙的弹力减去阻力
			XMFLOAT3 vec(0,0,1);
			float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

			XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

			m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

		}

		//Z-Wall max 
		fbounceParameter = WallThickness - (m_limitBox.max_z - m_ParticlePool[i].m_vPosRealScale.z);
		//到达墙壁范围内
		if(fbounceParameter > 0) 
		{
			//所受到的墙的弹力减去阻力
			XMFLOAT3 vec(0,0,-1);
			float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

			XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

			m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

		}


		XMFLOAT3 vAcceleration = MultiplyFS(m_ParticlePool[i].m_vAcceleration,testdeltatime); //at，增加的速度

		XMFLOAT3 vAccelerationAverage = DivideFS(vAcceleration,2.0f); // at/2.0f

		XMFLOAT3 vVelocityAverage = PlusF3(m_ParticlePool[i].m_vVelocity,vAccelerationAverage);// (v + at/2) 

		XMFLOAT3 vDisplacement = MultiplyFS(vVelocityAverage,testdeltatime); // vt + 1/2 * at^2

		m_ParticlePool[i].m_vVelocity = PlusF3(m_ParticlePool[i].m_vVelocity,vAcceleration); // v2 = v1 + at

		m_ParticlePool[i].m_vPosRealScale = PlusF3(m_ParticlePool[i].m_vPosRealScale,vDisplacement);

		m_ParticlePool[i].m_vPosRender = DivideFS(m_ParticlePool[i].m_vPosRealScale,m_unitScale);
	}
}






