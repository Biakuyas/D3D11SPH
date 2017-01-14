#ifndef VERTEXCONTROL_H_
#define VERTEXCONTROL_H_
#include"MyD3D11app.h"
#include"FXControl.h"
#include"SPH_System.h"
#include<vector>
#include<fstream>
#include<string>

struct Vertex
{
	XMFLOAT3 Pos;
	XMFLOAT3 Normal;
	XMFLOAT2 Tex;

	Vertex(
		float px, float py, float pz, 
		float u, float v)
		: Pos(px,py,pz),
		Tex(u,v){}
	Vertex(){}
};

namespace Colors{
	XMGLOBALCONST XMVECTORF32 White = {1.0,1.0,1.0,1.0};
	XMGLOBALCONST XMVECTORF32 Blue = {0,0,1.0,1.0};
	XMGLOBALCONST XMVECTORF32 Red = {1.0,0,0,1.0};
	XMGLOBALCONST XMVECTORF32 LightRed = {0.5f,0,0,1.0f};
	XMGLOBALCONST XMVECTORF32 LightBlue = {0.25f,0.25f,1.0f,1.0f};
}

class VertexIndexControl{

private:
		void Release();
public:
	VertexIndexControl(){}
	~VertexIndexControl(){Release();}


	static void VertexIndexInit(ID3D11Device * pd3dDevice);
	static HRESULT SetLayout(ID3D11Device * pd3dDevice);

	static ID3D11Buffer * m_IB;
	static ID3D11Buffer * m_VB;
	static ID3D11Buffer * m_VBSkull;
	static UINT m_iIndexCount;
	static ID3D11InputLayout* m_InputLayout;
};
#endif