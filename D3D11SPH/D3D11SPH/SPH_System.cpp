#include"SPH_System.h"

SPHSystem::SPHSystem()
{
	m_unitScale         = 0.004f;	// �ߴ絥λ
	m_viscosity			= 1.0f;				// ճ��
	m_restDensity		= 1000.f;			// �ܶ�
	m_pointMass			= 0.0004f;			// ��������
	m_gasConstantK		= 1.0f;				// �������巽�̳���
	m_smoothRadius		= 0.01f;			// �⻬�˰뾶
	m_wallLength        = 0.048f;
	m_updateinteval     = 0.003f;

	m_iParticleCount      = 0;
	m_iParticleMaxNumber = 500;
// 	m_XRecord = -10  * m_unitScale;
// 	m_ZRecord = -5 * m_unitScale;
	m_isBarrierActive = false;
//	m_iParticleCount      = 0;
	m_gravity = XMFLOAT3(0,-9.8f,0);

	//����
	m_limitBox = SPHBox(-m_wallLength ,0 ,-m_wallLength ,m_wallLength ,2 * m_wallLength ,m_wallLength);

//	m_limitBox = Box(0,0,0,0 * m_unitScale,2 * 0 * m_unitScale,100 * m_unitScale);
	m_boundaryCoefficient = 10000.f;
	m_boundaryDampening = 256.f;

	//Poly6 Kernel
	m_kernelPoly6 = 315.0f/(64.0f * PI * pow(m_smoothRadius, 9));
	//Spiky Kernel
	m_kernelSpiky = -45.0f/(PI * pow(m_smoothRadius, 6));
	//Viscosity Kernel
	m_kernelViscosity = 45.0f/(PI * pow(m_smoothRadius, 6));

	
	m_SPHRenderState = eBillPartilce;

}

void SPHSystem::SPH_CreateParticlePool()
{


	float pointDistance	= pow(m_pointMass/m_restDensity, 1.f/3.f); //���Ӽ��
//	_addFluidVolume(initFluidBox, pointDistance/m_unitScale);

	/*for(float i = m_limitBox.min_x; i <= m_limitBox.max_x; i += pointDistance)
		for(float j = m_limitBox.min_y; j <= m_limitBox.max_y; j += pointDistance)
			for(float p = m_limitBox.min_z; p <= m_limitBox.max_z; p += pointDistance)
			{
				
				m_ParticlePool[m_iParticleCount].m_vPosRender = XMFLOAT3(i / m_unitScale,j / m_unitScale,p / m_unitScale);
				m_ParticlePool[m_iParticleCount].m_vPosRealScale = XMFLOAT3(i,j + (10 * m_unitScale) ,p);
				m_ParticlePool[m_iParticleCount].m_fViscosity = 0;
				m_ParticlePool[m_iParticleCount].m_fDensity = 0;
				m_ParticlePool[m_iParticleCount].m_fPessure = 0;
				m_ParticlePool[m_iParticleCount].m_vAcceleration = XMFLOAT3(0,0,0);
				m_ParticlePool[m_iParticleCount].m_vVelocity = XMFLOAT3(0,0,0);
				//if(m_iParticleCount > 10000)
				//	MessageBoxA(0,"Over Pool",0,0);
				m_iParticleCount++;
			}*/

	//int ParticleNunmber_LmtJ = m_iParticleMaxNumber /((20 * m_unitScale / pointDistance) *(15 * m_unitScale / pointDistance))
	

	for(float j = WallThickness / 2.0 + m_wallLength / 2.0f; m_iParticleCount < m_iParticleMaxNumber  ; j += pointDistance)
	  for(float i = -10  * m_unitScale; i <= 10 * m_unitScale && m_iParticleCount < m_iParticleMaxNumber; i += pointDistance)
		for(float p = -5 * m_unitScale; p <= 10 * m_unitScale && m_iParticleCount < m_iParticleMaxNumber; p += pointDistance)
			{

// 				m_ParticlePool[m_iParticleCount].m_vPosRender = XMFLOAT3(i / m_unitScale,j / m_unitScale,p / m_unitScale);
// 				m_ParticlePool[m_iParticleCount].m_vPosRealScale = XMFLOAT3(i,j ,p);
// 				m_ParticlePool[m_iParticleCount].m_fViscosity = 0;
// 				m_ParticlePool[m_iParticleCount].m_fDensity = 0;
// 				m_ParticlePool[m_iParticleCount].m_fPessure = 0;
// 				m_ParticlePool[m_iParticleCount].m_vAcceleration = XMFLOAT3(0,0,0);
// 				m_ParticlePool[m_iParticleCount].m_vVelocity = XMFLOAT3(0,0,0);
				//if(m_iParticleCount > 10000)
				//	MessageBoxA(0,"Over Pool",0,0);

				Particle particle;

				particle.m_vPosRender = XMFLOAT3(i / m_unitScale,j / m_unitScale,p / m_unitScale);
				particle.m_vPosRealScale = XMFLOAT3(i,j ,p);
				particle.m_fViscosity = 0;
				particle.m_fDensity = 0;
				particle.m_fPessure = 0;
				particle.m_vAcceleration = XMFLOAT3(0,0,0);
				particle.m_vVelocity = XMFLOAT3(0,0,0);

				m_ParticlePool.push_back(particle);

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

void SPHSystem::SPH_DensityPressure(Grid & grids)
{
	//��α�����ʱ����Լ�¼һ����Χ�ĵ㣬�´ξͲ��ñ��������ˣ������ҹ��ƻ��ø����Ż��ķ�������
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
		m_ParticlePool[i].m_fDensity = m_kernelPoly6 * sum * m_pointMass; //��Ϊ�������ӵ�������һ��������������������
		m_ParticlePool[i].m_fPessure = (m_ParticlePool[i].m_fDensity - m_restDensity) * m_gasConstantK;
	}*/
	/*for(int i = 0;i < m_iParticleCount; i++)
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
		m_ParticlePool[i].m_fDensity = m_kernelPoly6 * sum * m_pointMass; //��Ϊ�������ӵ�������һ��������������������
		m_ParticlePool[i].m_fPessure = (m_ParticlePool[i].m_fDensity - m_restDensity) * m_gasConstantK;
	}*/
	float h2 = m_smoothRadius * m_smoothRadius;

	//vector<int>::iterator it;

	for(int i = 0;i < m_iParticleCount; i++)
	{
		float sum = 0;

		int neighborcells[9]; // ���2*2*2 + ;
		grids.findNeighborCells(m_ParticlePool[i].m_vPosRealScale,neighborcells);

		//ÿ��Ҫ����ڽӱ����¼���
		if(!m_ParticlePool[i].m_neightbortable.empty())

		m_ParticlePool[i].m_neightbortable.clear();

		for(int j = 0; j < 9; j++)
		{
			if(-1 == neighborcells[j])
				continue;

			
	
			for(int p = 0; p < grids.m_girdunit[ neighborcells[j] ].size(); p++)
			{
				int fParticleIndex = grids.m_girdunit[ neighborcells[j] ][p];
				float r2 = DistanceSquare(m_ParticlePool[i].m_vPosRealScale,m_ParticlePool[fParticleIndex].m_vPosRealScale);

				if(r2 > h2)
					continue;

				m_ParticlePool[i].m_neightbortable.push_back(fParticleIndex);

				sum += pow((h2 - r2),3);

			}

		}
		m_ParticlePool[i].m_fDensity = m_kernelPoly6 * sum * m_pointMass; //��Ϊ�������ӵ�������һ��������������������
		m_ParticlePool[i].m_fPessure = (m_ParticlePool[i].m_fDensity - m_restDensity) * m_gasConstantK;
			/*if(m_ParticlePool[i].m_fDensity < 1050)
		{
			XMFLOAT3 a = m_ParticlePool[i].m_vPosRealScale;
			int b;
			b = 3;
		}*/
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

	/* for(int i = 0;i < m_iParticleCount; i++)
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
	}*/
	
	float h2 = m_smoothRadius * m_smoothRadius;

	for(int i = 0;i < m_iParticleCount; i++)
	{
		XMFLOAT3 accel_sum = XMFLOAT3(0,0,0);

 		for(int j = 0; j < m_ParticlePool[i].m_neightbortable.size(); j++)
 		{
			int index = m_ParticlePool[i].m_neightbortable[j];

			float r2 = DistanceSquare(m_ParticlePool[i].m_vPosRealScale,m_ParticlePool[index].m_vPosRealScale);

			if(r2 > h2)
				continue;

			float r = Distance(m_ParticlePool[i].m_vPosRealScale,m_ParticlePool[index].m_vPosRealScale);


			float h2_r2 = h2 - r2;
			float h_r = m_smoothRadius - r;
			XMFLOAT3 ri_rj = MinusF3(m_ParticlePool[i].m_vPosRealScale,m_ParticlePool[index].m_vPosRealScale);

 			float pscalar = -m_pointMass * m_kernelSpiky * h_r * h_r * (m_ParticlePool[i].m_fPessure + m_ParticlePool[index].m_fPessure) / (2.0f * m_ParticlePool[i].m_fDensity * m_ParticlePool[index].m_fDensity);

			if(r != 0)
				accel_sum = PlusF3(accel_sum,DivideFS(MultiplyFS(ri_rj,pscalar),r) ); //(ri_rj * pscalar / r )


			float vscalar = m_kernelViscosity * m_viscosity * h_r * m_pointMass / (m_ParticlePool[i].m_fDensity * m_ParticlePool[index].m_fDensity);
			accel_sum = PlusF3(accel_sum,MultiplyFS(MinusF3(m_ParticlePool[index].m_vVelocity,m_ParticlePool[i].m_vVelocity),vscalar) );  //(vj - vi * vscalar)
 
		}
		m_ParticlePool[i].m_vAcceleration = accel_sum;
	}
}

void SPHSystem::SPH_ParticleUpdate(const float deltatime)
{
	for(int i = 0;i < m_iParticleCount;i++)
	{

	//	float testdeltatime = deltatime * 0.05f;//��1/1000ʱ��״̬�½���ģ�⣬��Ȼ�е㿴����
		float testdeltatime = m_updateinteval;
		//���������͵������¼��ٶ�
		m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,m_gravity);

		//һ���Ͽڲ���ĺ���

		//X-Wall min 
		float fbounceParameter = WallThickness - (m_ParticlePool[i].m_vPosRealScale.x - m_limitBox.min_x);
		//����ǽ�ڷ�Χ��
		if(fbounceParameter > 0) 
		{
			//���ܵ���ǽ�ĵ�����ȥ����
			XMFLOAT3 vec(1,0,0);
			float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

			XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

			m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

		}

		//X-Wall max 
		fbounceParameter = WallThickness - (m_limitBox.max_x - m_ParticlePool[i].m_vPosRealScale.x);
		//����ǽ�ڷ�Χ��
		if(fbounceParameter > 0) 
		{
			//���ܵ���ǽ�ĵ�����ȥ����
			XMFLOAT3 vec(-1,0,0);
			float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

			XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

			m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

		}

		//Y-Wall min 
		fbounceParameter =  WallThickness - (m_ParticlePool[i].m_vPosRealScale.y - m_limitBox.min_y);
		//����ǽ�ڷ�Χ��
		if(fbounceParameter > 0) 
		{
			//���ܵ���ǽ�ĵ�����ȥ����
			XMFLOAT3 vec(0,1,0);
			float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

			XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

			m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

		}

		//Z-Wall min 
		fbounceParameter = WallThickness - (m_ParticlePool[i].m_vPosRealScale.z - m_limitBox.min_z);
		//����ǽ�ڷ�Χ��
		if(fbounceParameter > 0) 
		{
			//���ܵ���ǽ�ĵ�����ȥ����
			XMFLOAT3 vec(0,0,1);
			float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

			XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

			m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

		}

		//Z-Wall max 
		fbounceParameter = WallThickness - (m_limitBox.max_z - m_ParticlePool[i].m_vPosRealScale.z);
		//����ǽ�ڷ�Χ��
		if(fbounceParameter > 0) 
		{
			//���ܵ���ǽ�ĵ�����ȥ����
			XMFLOAT3 vec(0,0,-1);
			float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

			XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

			m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

		}

		//���ϰ�ʱ�����
		if(m_isBarrierActive)
		{

			SPHBox m_Barrier = SPHBox(-m_wallLength / 3.0f ,0 ,-m_wallLength / 3.0f,m_wallLength / 3.0f ,2 * m_wallLength / 6.0f,m_wallLength / 3.0f);
			float m_boundaryDampening = 400;
			float m_boundaryCoefficient = 40000;
			float boundarylimit = 0.004f; //�������limit��Ŀ����Ϊ�˾�������ˮ������ǽ��ͬʱ������ײЧ������Ϊ��ȡ���ǵ������ܵĽ�������������п��ܷ�����
			//��������ȶ���Y > Z > X�����ܻᷢ��ײ��ǽ�ڽ��紦�������򲻶Ե�����
			//Y_max��ǽ��
			if(m_ParticlePool[i].m_vPosRealScale.x < m_Barrier.max_x && m_ParticlePool[i].m_vPosRealScale.x  > m_Barrier.min_x
				&& m_ParticlePool[i].m_vPosRealScale.z < m_Barrier.max_z && m_ParticlePool[i].m_vPosRealScale.z  > m_Barrier.min_z)
			{
				float fbounceParameter = m_Barrier.max_y - m_ParticlePool[i].m_vPosRealScale.y;
				//����ǽ�ڷ�Χ��
				if(fbounceParameter > 0 && fbounceParameter < boundarylimit) 
				{
					//���ܵ���ǽ�ĵ�����ȥ����
					XMFLOAT3 vec(0,1,0);
					float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

					XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

					m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

				}
			}

			//X_MAX��ǽ��
			if(m_ParticlePool[i].m_vPosRealScale.z < m_Barrier.max_z - boundarylimit && m_ParticlePool[i].m_vPosRealScale.z  > m_Barrier.min_z + boundarylimit
				&& m_ParticlePool[i].m_vPosRealScale.y < m_Barrier.max_y - boundarylimit && m_ParticlePool[i].m_vPosRealScale.x > 0)
			{

				float fbounceParameter = m_Barrier.max_x - m_ParticlePool[i].m_vPosRealScale.x;
				//����ǽ�ڷ�Χ��
				if(fbounceParameter > 0 && fbounceParameter < boundarylimit) 
				{
					//���ܵ���ǽ�ĵ�����ȥ����
					XMFLOAT3 vec(1,0,0);
					float fbounceForce = m_boundaryCoefficient  * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

					XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

					m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

				}


			}

			//X_MIN��ǽ��
			if(m_ParticlePool[i].m_vPosRealScale.z < m_Barrier.max_z - boundarylimit && m_ParticlePool[i].m_vPosRealScale.z  > m_Barrier.min_z + boundarylimit
				&& m_ParticlePool[i].m_vPosRealScale.y < m_Barrier.max_y - boundarylimit && m_ParticlePool[i].m_vPosRealScale.x < 0)
			{

				float fbounceParameter =  m_ParticlePool[i].m_vPosRealScale.x - m_Barrier.min_x;
				//����ǽ�ڷ�Χ��
				if(fbounceParameter > 0 && fbounceParameter < boundarylimit)  
				{
					//���ܵ���ǽ�ĵ�����ȥ����
					XMFLOAT3 vec(-1,0,0);
					float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

					XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

					m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

				}
// 				if(m_ParticlePool[i].m_vPosRealScale.x >= m_Barrier.min_x)
// 				{
// 					//m_ParticlePool[i].m_vPosRealScale.x = m_Barrier.min_x;
// 					m_ParticlePool[i].m_vVelocity.x = 0;
// 				}
			}

			//Z_Max��ǽ��
			if(m_ParticlePool[i].m_vPosRealScale.x < m_Barrier.max_x && m_ParticlePool[i].m_vPosRealScale.x  > m_Barrier.min_x
				&& m_ParticlePool[i].m_vPosRealScale.y < m_Barrier.max_y - boundarylimit && m_ParticlePool[i].m_vPosRealScale.z > 0)
			{
				float fbounceParameter = m_Barrier.max_z - m_ParticlePool[i].m_vPosRealScale.z;
				//����ǽ�ڷ�Χ��
				if(fbounceParameter > 0 && fbounceParameter < boundarylimit) 
				{
					//���ܵ���ǽ�ĵ�����ȥ����
					XMFLOAT3 vec(0,0,1);
					float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

					XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

					m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

				}
// 				if(m_ParticlePool[i].m_vPosRealScale.z <= m_Barrier.max_z )
// 				{
// 					//m_ParticlePool[i].m_vPosRealScale.z = m_Barrier.max_z;
// 					m_ParticlePool[i].m_vVelocity.z = 0;
// 				}
			}

			//Z_Min��ǽ��
			if(m_ParticlePool[i].m_vPosRealScale.x < m_Barrier.max_x && m_ParticlePool[i].m_vPosRealScale.x  > m_Barrier.min_x
				&& m_ParticlePool[i].m_vPosRealScale.y < m_Barrier.max_y - boundarylimit && m_ParticlePool[i].m_vPosRealScale.z < 0)
			{

				float fbounceParameter =  m_ParticlePool[i].m_vPosRealScale.z - m_Barrier.min_z;
				//����ǽ�ڷ�Χ��
				if(fbounceParameter > 0 && fbounceParameter < boundarylimit) 
				{
					//���ܵ���ǽ�ĵ�����ȥ����
					XMFLOAT3 vec(0,0,-1);
					float fbounceForce = m_boundaryCoefficient * fbounceParameter - m_boundaryDampening * DotF3(vec,m_ParticlePool[i].m_vVelocity);

					XMFLOAT3 vBounceAcc = MultiplyFS(vec,fbounceForce);

					m_ParticlePool[i].m_vAcceleration = PlusF3(m_ParticlePool[i].m_vAcceleration,vBounceAcc);

				}

// 				if(m_ParticlePool[i].m_vPosRealScale.z > m_Barrier.min_z )
// 				{
// 	//				m_ParticlePool[i].m_vPosRealScale.z = m_Barrier.min_z;
// 					m_ParticlePool[i].m_vVelocity.z = 0;
// 				}
			}





		}


		XMFLOAT3 vAcceleration = MultiplyFS(m_ParticlePool[i].m_vAcceleration,testdeltatime); //at�����ӵ��ٶ�

		XMFLOAT3 vAccelerationAverage = DivideFS(vAcceleration,2.0f); // at/2.0f

		XMFLOAT3 vVelocityAverage = PlusF3(m_ParticlePool[i].m_vVelocity,vAccelerationAverage);// (v + at/2) 

		XMFLOAT3 vDisplacement = MultiplyFS(vVelocityAverage,testdeltatime); // vt + 1/2 * at^2

		m_ParticlePool[i].m_vVelocity = PlusF3(m_ParticlePool[i].m_vVelocity,vAcceleration); // v2 = v1 + at

		m_ParticlePool[i].m_vPosRealScale = PlusF3(m_ParticlePool[i].m_vPosRealScale,vDisplacement);

		m_ParticlePool[i].m_vPosRender = DivideFS(m_ParticlePool[i].m_vPosRealScale,m_unitScale);
	}
}

void SPHSystem::SPH_Reset()
{
	m_iParticleCount      = 0;
	if(!m_ParticlePool.empty())
		m_ParticlePool.clear();
}

void SPHSystem::SPH_Init()
{

	m_unitScale         = 0.004f;	// �ߴ絥λ
	m_viscosity			= 1.0f;				// ճ��
	m_restDensity		= 1000.f;			// �ܶ�
	m_pointMass			= 0.0004f;			// ��������
	m_gasConstantK		= 1.0f;				// �������巽�̳���
	m_smoothRadius		= 0.01f;			// �⻬�˰뾶
	m_wallLength        = 0.048f;
	m_updateinteval     = 0.003f;

	m_iParticleCount      = 0;
	m_iParticleMaxNumber = 500;
	// 	m_XRecord = -10  * m_unitScale;
	// 	m_ZRecord = -5 * m_unitScale;

	//	m_iParticleCount      = 0;
	m_gravity = XMFLOAT3(0,-9.8f,0);

	//����
	m_limitBox = SPHBox(-m_wallLength ,0 ,-m_wallLength ,m_wallLength ,2 * m_wallLength ,m_wallLength);

	//	m_limitBox = Box(0,0,0,0 * m_unitScale,2 * 0 * m_unitScale,100 * m_unitScale);
	m_boundaryCoefficient = 10000.f;
	m_boundaryDampening = 256.f;
}






