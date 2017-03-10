#include "StateControl.h"
ID3D11BlendState* StateControl::m_TransparentBS = 0;
//D3D11BlendState* StateControl::m_NoColorWriteBS = 0;

//ID3D11DepthStencilState* StateControl::m_DSSMirror = 0;
//ID3D11DepthStencilState* StateControl::m_DSSSkullReflection = 0;
//ID3D11DepthStencilState* StateControl::m_DSSShadow = 0;

ID3D11RasterizerState* StateControl::m_NoBackSolid = 0;
//ID3D11RasterizerState* StateControl::m_BackCullRsClockWise = 0;
ID3D11RasterizerState* StateControl::m_NoBackWireFrame = 0;

HRESULT StateControl::InitAllState(ID3D11Device * pd3dDevice)
{
//BlendState
//m_TransparentBS
	D3D11_BLEND_DESC transparent_blendDESC = {0};

	transparent_blendDESC.AlphaToCoverageEnable = false;
	transparent_blendDESC.IndependentBlendEnable = false;

	transparent_blendDESC.RenderTarget[0].BlendEnable = true;
	transparent_blendDESC.RenderTarget[0].BlendOp = D3D11_BLEND_OP_ADD;
	transparent_blendDESC.RenderTarget[0].DestBlend = D3D11_BLEND_INV_SRC_ALPHA;//
	transparent_blendDESC.RenderTarget[0].SrcBlend = D3D11_BLEND_SRC_ALPHA;
	transparent_blendDESC.RenderTarget[0].BlendOpAlpha = D3D11_BLEND_OP_ADD;
	transparent_blendDESC.RenderTarget[0].DestBlendAlpha = D3D11_BLEND_ZERO;
	transparent_blendDESC.RenderTarget[0].SrcBlendAlpha = D3D11_BLEND_ONE;
	transparent_blendDESC.RenderTarget[0].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;

	HRESULT hr = pd3dDevice->CreateBlendState(&transparent_blendDESC,&m_TransparentBS);

//m_NoColorWriteBS
	/*D3D11_BLEND_DESC noColorBlend_DESC = {0};

	noColorBlend_DESC.RenderTarget[0].BlendEnable = false;
	noColorBlend_DESC.RenderTarget[0].BlendOp = D3D11_BLEND_OP_ADD;
	noColorBlend_DESC.RenderTarget[0].DestBlend = D3D11_BLEND_ONE;
	noColorBlend_DESC.RenderTarget[0].SrcBlend = D3D11_BLEND_ZERO;
	noColorBlend_DESC.RenderTarget[0].BlendOpAlpha = D3D11_BLEND_OP_ADD;
	noColorBlend_DESC.RenderTarget[0].DestBlendAlpha = D3D11_BLEND_ONE;
	noColorBlend_DESC.RenderTarget[0].SrcBlendAlpha = D3D11_BLEND_ZERO;
	noColorBlend_DESC.RenderTarget[0].RenderTargetWriteMask = 0;

	pd3dDevice->CreateBlendState(&noColorBlend_DESC,&m_NoColorWriteBS);
// DepthStencilState
//MirrorDSS
	D3D11_DEPTH_STENCILOP_DESC mirrorstencilop;
	mirrorstencilop.StencilFailOp = D3D11_STENCIL_OP_KEEP;
	mirrorstencilop.StencilDepthFailOp = D3D11_STENCIL_OP_KEEP;
	mirrorstencilop.StencilPassOp = D3D11_STENCIL_OP_REPLACE;
	mirrorstencilop.StencilFunc = D3D11_COMPARISON_ALWAYS;

	D3D11_DEPTH_STENCIL_DESC mirrorstencil;
	mirrorstencil.DepthEnable = TRUE;
	mirrorstencil.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ZERO;
	mirrorstencil.DepthFunc = D3D11_COMPARISON_LESS;

	mirrorstencil.StencilEnable = TRUE;
	mirrorstencil.StencilReadMask = 0xff;
	mirrorstencil.StencilWriteMask = 0xff;
	mirrorstencil.FrontFace = mirrorstencilop;
	mirrorstencil.BackFace = mirrorstencilop;

	hr = pd3dDevice->CreateDepthStencilState(&mirrorstencil,&m_DSSMirror);
	if(FAILED(hr))
	{
		MessageBoxA(0,"Create DepthStencilState Failed",0,0);
		return E_FAIL;
	}
//ReflectionDSS
	D3D11_DEPTH_STENCILOP_DESC reflectionstencilop;
	reflectionstencilop.StencilFailOp = D3D11_STENCIL_OP_KEEP;
	reflectionstencilop.StencilDepthFailOp = D3D11_STENCIL_OP_KEEP;
	reflectionstencilop.StencilPassOp = D3D11_STENCIL_OP_KEEP;
	reflectionstencilop.StencilFunc = D3D11_COMPARISON_EQUAL;

	D3D11_DEPTH_STENCIL_DESC reflectinostencil;

	reflectinostencil.DepthEnable = TRUE;
	reflectinostencil.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL; //attention
	reflectinostencil.DepthFunc = D3D11_COMPARISON_LESS;

	reflectinostencil.StencilEnable = TRUE;
	reflectinostencil.StencilReadMask = 0xff;
	reflectinostencil.StencilWriteMask = 0xff;
	reflectinostencil.FrontFace = reflectionstencilop;
	reflectinostencil.BackFace = reflectionstencilop;

	hr = pd3dDevice->CreateDepthStencilState(&reflectinostencil,&m_DSSSkullReflection);
	if(FAILED(hr))
	{
		MessageBoxA(0,"Create DepthStencilState Failed",0,0);
		return E_FAIL;
	}
//ShadowDSS
	D3D11_DEPTH_STENCILOP_DESC shadowopdesc;
	shadowopdesc.StencilFailOp = D3D11_STENCIL_OP_KEEP;
	shadowopdesc.StencilDepthFailOp = D3D11_STENCIL_OP_KEEP;
	shadowopdesc.StencilPassOp = D3D11_STENCIL_OP_INCR;
	shadowopdesc.StencilFunc = D3D11_COMPARISON_EQUAL;

	D3D11_DEPTH_STENCIL_DESC shadowdesc;

	shadowdesc.DepthEnable = true;
	shadowdesc.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL;
	shadowdesc.DepthFunc = D3D11_COMPARISON_LESS;

	shadowdesc.StencilEnable = TRUE;
	shadowdesc.StencilReadMask = 0xff;
	shadowdesc.StencilWriteMask = 0xff;
	shadowdesc.FrontFace = shadowopdesc;
	shadowdesc.BackFace = shadowopdesc;

	hr = pd3dDevice->CreateDepthStencilState(&shadowdesc,&m_DSSShadow);
	if(FAILED(hr))
	{
		MessageBoxA(0,"Create DepthStencilState Failed",0,0);
		return E_FAIL;
	}
*/
//Rasterizer State
	//m_BackCullRs
	D3D11_RASTERIZER_DESC mRasterizerDesc;
	ZeroMemory(&mRasterizerDesc,sizeof(D3D11_RASTERIZER_DESC));
	mRasterizerDesc.CullMode = D3D11_CULL_NONE;
	mRasterizerDesc.FillMode = D3D11_FILL_SOLID;
	mRasterizerDesc.DepthClipEnable = TRUE;


	hr = pd3dDevice->CreateRasterizerState(&mRasterizerDesc,&m_NoBackSolid);
	if(FAILED(hr))
	{
		MessageBoxA(0,"Fail to Create Rasterizer State",0,0);
		return E_FAIL;
	}


	//m_BackCullRsClockWise
// 	D3D11_RASTERIZER_DESC mRasterizerWiseDesc;
// 	ZeroMemory(&mRasterizerWiseDesc,sizeof(D3D11_RASTERIZER_DESC));
// 	mRasterizerWiseDesc.CullMode = D3D11_CULL_BACK;
// 	mRasterizerWiseDesc.FillMode = D3D11_FILL_SOLID;
// 	mRasterizerWiseDesc.DepthClipEnable = TRUE;
// 	mRasterizerWiseDesc.FrontCounterClockwise = TRUE;
// 
// 
// 	hr = pd3dDevice->CreateRasterizerState(&mRasterizerWiseDesc,&m_BackCullRsClockWise);
// 	if(FAILED(hr))
// 	{
// 		MessageBoxA(0,"Fail to Create Rasterizer State",0,0);
// 		return E_FAIL;
// 	}



	D3D11_RASTERIZER_DESC mWireRasterizerDesc;
	ZeroMemory(&mWireRasterizerDesc,sizeof(D3D11_RASTERIZER_DESC));
	mWireRasterizerDesc.CullMode = D3D11_CULL_NONE;
	mWireRasterizerDesc.FillMode = D3D11_FILL_WIREFRAME;
	mWireRasterizerDesc.DepthClipEnable = TRUE;


	hr = pd3dDevice->CreateRasterizerState(&mWireRasterizerDesc,&m_NoBackWireFrame);
	if(FAILED(hr))
	{
		MessageBoxA(0,"Fail to Create Rasterizer State",0,0);
		return E_FAIL;
	}
	return S_OK;

}

void StateControl::Release()
{


}


