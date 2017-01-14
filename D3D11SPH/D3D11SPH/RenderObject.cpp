#include "RenderObject.h"

HRESULT RenderObject::ObjectInit()
{
	//VertexIndexInit();

	VertexIndexControl::VertexIndexInit(m_pd3dDevice);

	HRESULT hr;

	hr = FXControl::Init(m_pd3dDevice);
	if(FAILED(hr))
		return E_FAIL;


	hr = VertexIndexControl::SetLayout(m_pd3dDevice);
	if(FAILED(hr))
		return E_FAIL;


	hr = StateControl::InitAllState(m_pd3dDevice);
	if(FAILED(hr))
		return E_FAIL;


	// Read Texture From File
	 hr = D3DX11CreateShaderResourceViewFromFile(m_pd3dDevice,L"Textures/brick01.dds",0,0,&m_d3d11SRVFloor,0);

	if(FAILED(hr))
	{
		MessageBoxA(0,"Load Texture Failed",0,0);
		return E_FAIL;
	}

	hr = D3DX11CreateShaderResourceViewFromFile(m_pd3dDevice,L"Textures/checkboard.dds",0,0,&m_d3d11SRVWall,0);

	if(FAILED(hr))
	{
		MessageBoxA(0,"Load Texture Failed",0,0);
		return E_FAIL;
	}


	hr = D3DX11CreateShaderResourceViewFromFile(m_pd3dDevice,L"Textures/ice.dds",0,0,&m_d3d11SRVMirror,0);

	if(FAILED(hr))
	{
		MessageBoxA(0,"Load Texture Failed",0,0);
		return E_FAIL;
	}

	m_SPHSystem.SPH_Init();

	return S_OK;



}


void RenderObject::Release()
{
	MyD3D11app::Release();

}

RenderObject::RenderObject():MyD3D11app()
{
	//MyD3D11app::MyD3D11app();

	m_fRadius = 100;

	m_CameraPosx = 0;
	m_CameraPosy = 0;
	m_CameraPosz = -m_fRadius;

	m_MousePosLastx = 0;
	m_MousePosLasty = 0;


	m_Phi = 0;
	m_Thea = 0;

	XMMATRIX IdenticalMatrix = XMMatrixIdentity();
	m_matWorldTranslate = IdenticalMatrix;
	m_matWorldRotate = IdenticalMatrix;
	m_matWorldScailing = IdenticalMatrix;

	m_v3WorldCoordinate = XMFLOAT3(0,0,0);
	m_v3Rotation = XMFLOAT3(0,0,0);
	m_v3Scailing = XMFLOAT3(1,1,1);

	m_vLightPos = XMFLOAT3(-2.5,10,-2.5);

	m_v3WorldCoordinate_Skull = XMFLOAT3(0,0,-4);
	m_v3Rotation_Skull = XMFLOAT3(0,90,0);
	m_v3Scailing_Skull = XMFLOAT3(0.5,0.5,0.5);

}

RenderObject::~RenderObject()
{
	Release();

}

int RenderObject::Init(HINSTANCE hInstance,int nShowCmd)
{
	if(-1 == MyD3D11app::Init(hInstance,nShowCmd))
		return -1;
	ObjectInit();

	return 0;

}

void RenderObject::UpdateScene(const float dt)
{
//matWorld
	m_matWorldTranslate = XMMatrixTranslation(m_v3WorldCoordinate.x,m_v3WorldCoordinate.y,m_v3WorldCoordinate.z);
	m_matWorldRotate =  XMMatrixRotationX(m_v3Rotation.x) * XMMatrixRotationY(m_v3Rotation.y);
	m_matWorldScailing = XMMatrixScaling(m_v3Scailing.x,m_v3Scailing.y,m_v3Scailing.z);

	//m_matWorld = m_matWorldScailing * m_matWorldRotate * m_matWorldTranslate;


// 	XMMATRIX matWorldTranslate_Skull = XMMatrixTranslation(m_v3WorldCoordinate_Skull.x,m_v3WorldCoordinate_Skull.y,m_v3WorldCoordinate_Skull.z);
// 	XMMATRIX matWorldRotate_Skull =  XMMatrixRotationX(m_v3Rotation_Skull.x) * XMMatrixRotationY(m_v3Rotation_Skull.y);
// 	XMMATRIX matWorldScailing_Skull = XMMatrixScaling(m_v3Scailing_Skull.x,m_v3Scailing_Skull.y,m_v3Scailing_Skull.z);
// 
// 	m_matSkullWorld = matWorldScailing_Skull * matWorldRotate_Skull * matWorldTranslate_Skull;

//matView
	if(m_Phi > PI - 0.1)
		m_Phi = PI - 0.1;
	if(m_Phi < 0.1f)
		m_Phi = 0.1f;

	float x = m_fRadius*sinf(m_Phi)*cosf(m_Thea);
	float z = m_fRadius*sinf(m_Phi)*sinf(m_Thea);
	float y = m_fRadius*cosf(m_Phi);

	m_CameraPosx = x;
	m_CameraPosy = y;
	m_CameraPosz = z;

	XMVECTOR pos = XMVectorSet(m_CameraPosx,m_CameraPosy,m_CameraPosz,1.0);
	XMVECTOR target = XMVectorSet(0,0,0,1.0);
	XMVECTOR up = XMVectorSet(0,1.0f,0,1.0);

	m_matView = XMMatrixLookAtLH(pos,target,up);

	//切线空间
	XMMATRIX matWorldInverseTranspose;
	matWorldInverseTranspose = MatrixInverseTranspose(m_matWorld);



//matWorldViewProj
	//m_matWorldViewProj = m_matWorld * m_matView * m_matProj;


	//摄像机位置
	XMFLOAT3 EyePos(m_CameraPosx,m_CameraPosy,m_CameraPosz);
	//平行光方向
	XMFLOAT3 LightDirection(m_v3WorldCoordinate.x - m_vLightPos.x, m_v3WorldCoordinate.y - m_vLightPos.y, m_v3WorldCoordinate.z - m_vLightPos.z);


	//setFX variable
	//FXControl::m_fxMatWorldViewProj->SetMatrix((float*)&m_matWorldViewProj);
	//FXControl::m_fxMatWorld->SetMatrix((float*)&m_matWorld);
	//FXControl::m_fxMatWorldInverseTranspose->SetMatrix((float*)&matWorldInverseTranspose);
	FXControl::m_fxEyePos->SetRawValue((float*)&EyePos,0,sizeof(XMFLOAT3));
	FXControl::m_fxLightDir->SetRawValue((float*)&LightDirection,0,sizeof(LightDirection));

	m_SPHSystem.SPH_DensityPressure();
	m_SPHSystem.SPH_AccelerationFigure();
	m_SPHSystem.SPH_ParticleUpdate(dt);


}

HRESULT RenderObject::Resize(float window_width,float window_height)
{
	if(FAILED(MyD3D11app::Resize(window_width,window_height)))
		return E_FAIL;

	m_matProj = XMMatrixPerspectiveFovLH(PI/4.0f,(float)window_width/(float)window_height,1.0f,1000.0f);
}

void RenderObject::DrawScene()
{
	//Start to Draw
	m_pd3dDevContext->ClearRenderTargetView(m_pRenderTargetView,(const float*)&Colors::LightBlue);
	m_pd3dDevContext->ClearDepthStencilView(m_pDepthStencilView,D3D11_CLEAR_DEPTH | D3D11_CLEAR_STENCIL,1.0,0);

	m_pd3dDevContext->IASetInputLayout(VertexIndexControl::m_InputLayout);
	m_pd3dDevContext->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_POINTLIST);
//Common Variable
	UINT iVertexStrides = sizeof(Vertex);
	UINT ioffset = 0;
	float fBlendFactors[] = {0,0,0,0};
	XMMATRIX matWorldInverseTranspose;
//Technique
	D3DX11_TECHNIQUE_DESC techDesc;
	FXControl::m_fxTech->GetDesc(&techDesc);


//Different Object Render
// Set Buffer
	m_pd3dDevContext->IASetVertexBuffers(0,1,&VertexIndexControl::m_VB,&iVertexStrides,&ioffset);
	//m_pd3dDevContext->IASetIndexBuffer(VertexIndexControl::m_IB,DXGI_FORMAT_R32_UINT,0);
	// Render 500 Particle
	for(int i = 0;i < m_SPHSystem.m_iParticleCount; i++)
	{
		//Set State
		m_pd3dDevContext->RSSetState(StateControl::m_BackCullRs);
		//Set Matrix
		//matWorldInverseTranspose = MatrixInverseTranspose(m_matWorld);
		m_matWorld = XMMatrixTranslation(m_SPHSystem.m_ParticlePool[i].m_vPosRender.x,m_SPHSystem.m_ParticlePool[i].m_vPosRender.y,m_SPHSystem.m_ParticlePool[i].m_vPosRender.z);
		m_matWorldViewProj = m_matWorld * m_matView * m_matProj;
		//Set FX Variable
		FXControl::m_fxMatWorldViewProj->SetMatrix((float*)&m_matWorldViewProj);
		FXControl::m_fxMatWorld->SetMatrix((float*)&m_matWorld);
		FXControl::m_fxMatWorldInverseTranspose->SetMatrix((float*)&matWorldInverseTranspose);

		//Particle
		//Set Texture
		//FXControl::m_fxDiffuseMap->SetResource(m_d3d11SRVFloor);
		for(int p = 0; p < techDesc.Passes; p++ )
		{
			FXControl::m_fxTech->GetPassByIndex(p)->Apply(0,m_pd3dDevContext);
			m_pd3dDevContext->Draw(1,i);
		}
	}

	/*//Wall
	//Set Texture
	FXControl::m_fxDiffuseMap->SetResource(m_d3d11SRVWall);
	for(int p = 0; p < techDesc.Passes; p++ )
	{
		FXControl::m_fxTech->GetPassByIndex(p)->Apply(0,m_pd3dDevContext);
		m_pd3dDevContext->Draw(18,6);
	}

	//Skull
//SetBuffer
	m_pd3dDevContext->IASetVertexBuffers(0,1,&VertexIndexControl::m_VBSkull,&iVertexStrides,&ioffset);
//Set Matrix
	m_matWorldViewProj = m_matSkullWorld * m_matView * m_matProj;
	matWorldInverseTranspose = MatrixInverseTranspose(m_matSkullWorld);
// Set FXVariable
	FXControl::m_fxMatWorldViewProj->SetMatrix((float*)&m_matWorldViewProj);
	FXControl::m_fxMatWorld->SetMatrix((float*)&m_matSkullWorld);
	FXControl::m_fxMatWorldInverseTranspose->SetMatrix((float*)&matWorldInverseTranspose);

	for(int p = 0; p < techDesc.Passes; p++ )
	{
		FXControl::m_fxTech->GetPassByIndex(p)->Apply(0,m_pd3dDevContext);
		m_pd3dDevContext->DrawIndexed(VertexIndexControl::m_iIndexCount,0,0);
	}


	// Mirror Stencil
//SetBuffer
	m_pd3dDevContext->IASetVertexBuffers(0,1,&VertexIndexControl::m_VB,&iVertexStrides,&ioffset);
//Set State
	m_pd3dDevContext->OMSetBlendState(StateControl::m_NoColorWriteBS,fBlendFactors,0xffffffff);
	m_pd3dDevContext->OMSetDepthStencilState(StateControl::m_DSSMirror,1);
//SetMatrix
	m_matWorldViewProj = m_matWorld * m_matView * m_matProj;
//Set FX Variable
	FXControl::m_fxMatWorldViewProj->SetMatrix((float*)&m_matWorldViewProj);
	FXControl::m_fxMatWorld->SetMatrix((float*)&m_matWorld);

	for(int p = 0; p < techDesc.Passes; p++ )
	{
		FXControl::m_fxTech->GetPassByIndex(p)->Apply(0,m_pd3dDevContext);
		m_pd3dDevContext->Draw(6,24);
	}
//Set Default State
	m_pd3dDevContext->OMSetBlendState(0,fBlendFactors,0xffffffff);
	m_pd3dDevContext->OMSetDepthStencilState(0,0);


	//reflected Skull
//Set Buffer
	m_pd3dDevContext->IASetVertexBuffers(0,1,&VertexIndexControl::m_VBSkull,&iVertexStrides,&ioffset);
//Set State
	m_pd3dDevContext->OMSetDepthStencilState(StateControl::m_DSSSkullReflection,1);
	m_pd3dDevContext->RSSetState(StateControl::m_BackCullRsClockWise);
//Set Matrix
	XMVECTOR mirrorPlane = XMVectorSet(0,0,1.0,0);
	XMMATRIX matReflect = XMMatrixReflect(mirrorPlane);
	XMMATRIX matSkullReflectWorld = m_matSkullWorld * matReflect;

	m_matWorldViewProj = matSkullReflectWorld * m_matView * m_matProj;

//Set FXVariable
	FXControl::m_fxMatWorldViewProj->SetMatrix((float*)&m_matWorldViewProj);
	FXControl::m_fxMatWorld->SetMatrix((float*)&matSkullReflectWorld);

	for(int p = 0; p < techDesc.Passes; p++ )
	{
		FXControl::m_fxTech->GetPassByIndex(p)->Apply(0,m_pd3dDevContext);
		m_pd3dDevContext->DrawIndexed(VertexIndexControl::m_iIndexCount,0,0);
	}
// Set Default State
	m_pd3dDevContext->OMSetDepthStencilState(0,0);
	m_pd3dDevContext->RSSetState(StateControl::m_BackCullRs);


	//Mirror
//Set Buffer
	m_pd3dDevContext->IASetVertexBuffers(0,1,&VertexIndexControl::m_VB,&iVertexStrides,&ioffset);
//Set State
	m_pd3dDevContext->OMSetBlendState(StateControl::m_TransparentBS,fBlendFactors,0xffffffff);
//Set Matrix
	m_matWorldViewProj = m_matWorld * m_matView * m_matProj;
	matWorldInverseTranspose = MatrixInverseTranspose(m_matWorld);
//Set FXVariable
	FXControl::m_fxMatWorldViewProj->SetMatrix((float*)&m_matWorldViewProj);
	FXControl::m_fxMatWorld->SetMatrix((float*)&m_matWorld);
	FXControl::m_fxMatWorldInverseTranspose->SetMatrix((float*)&matWorldInverseTranspose);


// Set Different Technique
	//transparent technique
	D3DX11_TECHNIQUE_DESC transparentTech;
	FXControl::m_fxTransparentTech->GetDesc(&transparentTech);

// Set Texture
	FXControl::m_fxDiffuseMap->SetResource(m_d3d11SRVMirror);
	for(int p = 0; p < transparentTech.Passes; p++ )
	{
		FXControl::m_fxTransparentTech->GetPassByIndex(p)->Apply(0,m_pd3dDevContext);
		m_pd3dDevContext->Draw(6,24);
	}
//Set Default State
	m_pd3dDevContext->OMSetBlendState(0,fBlendFactors,0xffffffff);


	// Shadow
//Set Buffer
    m_pd3dDevContext->IASetVertexBuffers(0,1,&VertexIndexControl::m_VBSkull,&iVertexStrides,&ioffset);
//Set State
	m_pd3dDevContext->OMSetBlendState(StateControl::m_TransparentBS,fBlendFactors,0xffffffff);
	m_pd3dDevContext->OMSetDepthStencilState(StateControl::m_DSSShadow,0);
//Set Matrix
	XMVECTOR shadowplane = XMVectorSet(0,1.0f,0,0);
	XMVECTOR LightPosition = XMVectorSet(m_vLightPos.x,m_vLightPos.y,m_vLightPos.z,1.0f);
	XMMATRIX matshadow = XMMatrixShadow(shadowplane,LightPosition);

	XMMATRIX shadowYtranslation = XMMatrixTranslation(0,0.01,0);

	XMMATRIX shadowWorld = m_matSkullWorld  * matshadow * shadowYtranslation;
	m_matWorldViewProj = shadowWorld * m_matView * m_matProj;

//Set FXVariable
	FXControl::m_fxMatWorldViewProj->SetMatrix((float*)&m_matWorldViewProj);
	FXControl::m_fxMatWorld->SetMatrix((float*)&shadowWorld);

//Set Different Tech
	D3DX11_TECHNIQUE_DESC shadowTechDesc;
	FXControl::m_fxShadowTech->GetDesc(&shadowTechDesc);

	for(int p = 0; p < shadowTechDesc.Passes; p++ )
	{
		FXControl::m_fxShadowTech->GetPassByIndex(p)->Apply(0,m_pd3dDevContext);
		m_pd3dDevContext->DrawIndexed(VertexIndexControl::m_iIndexCount,0,0);
	}*/
//Set Default State
	m_pd3dDevContext->OMSetBlendState(0,fBlendFactors,0xffffffff);
	m_pd3dDevContext->OMSetDepthStencilState(0,0);

	HRESULT hr = m_pSwapChain->Present(0,0);

	if(FAILED(hr))
	{
		MessageBoxA(0,"Fail to Swap Chain",0,0);
	}
}

void RenderObject::KeyDown(WPARAM wParam)
{
	if(wParam == 'W')
		m_v3WorldCoordinate.z += 0.05;
	if(wParam == 'D')
		m_v3WorldCoordinate.x += 0.05;
	if(wParam == 'S')
		m_v3WorldCoordinate.z -= 0.05;
	if(wParam == 'A')
		m_v3WorldCoordinate.x -= 0.05;

	if(wParam == 'K')
		m_SPHSystem.SPH_Init();

}

void RenderObject::MouseButtonDown(WPARAM wParam,float fMousePointx,float fMousePointy)
{
	m_MousePosLastx = fMousePointx;
	m_MousePosLasty = fMousePointy;
}

void RenderObject::MouseMove(WPARAM wParam,float fMousePointx,float fMousePointy)
{

	if( (wParam & MK_LBUTTON) != 0 )
	{

		m_MousePosNowx = fMousePointx;
		m_MousePosNowy = fMousePointy;
		// Make each pixel correspond to a quarter of a degree.
		float dx = XMConvertToRadians(0.25f*static_cast<float>(m_MousePosNowx - m_MousePosLastx));
		float dy = XMConvertToRadians(0.25f*static_cast<float>(m_MousePosNowy - m_MousePosLasty));

		// Update angles based on input to orbit camera around box.
		m_Thea += dx;
		m_Phi  += dy;

		if(m_Phi > PI - 0.1)
			m_Phi = PI - 0.1;
		if(m_Phi < 0.1f)
			m_Phi = 0.1f;


		// Restrict the angle mPhi.
	}
	if( (wParam & MK_RBUTTON) != 0 )
	{
		m_MousePosNowx = fMousePointx;
		m_MousePosNowy = fMousePointy;

		// Make each pixel correspond to 0.01 unit in the scene.
		float dx = 0.05f*static_cast<float>(m_MousePosNowx - m_MousePosLastx);
		float dy = 0.05f*static_cast<float>(m_MousePosNowy - m_MousePosLasty);

		// Update the camera radius based on input.
		m_fRadius += dx - dy;

	}

	m_MousePosLastx = m_MousePosNowx;
	m_MousePosLasty = m_MousePosNowy;

}



XMMATRIX RenderObject::MatrixInverseTranspose(XMMATRIX& matrix)
{
	XMMATRIX matTemp = matrix;
	//matTemp.r[3] = XMVectorSet(0,0,0,1.0f); //这里不知道为什么

	XMVECTOR det = XMMatrixDeterminant(matTemp);

	return XMMatrixTranspose(XMMatrixInverse(&det,matTemp));
}
