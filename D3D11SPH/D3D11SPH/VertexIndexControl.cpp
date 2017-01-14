#include"VertexIndexControl.h"


ID3D11Buffer * VertexIndexControl::m_VBSkull = 0;

ID3D11Buffer * VertexIndexControl::m_VB = 0;

ID3D11Buffer * VertexIndexControl::m_IB = 0;

UINT VertexIndexControl::m_iIndexCount = 0;

ID3D11InputLayout* VertexIndexControl::m_InputLayout = 0;


void VertexIndexControl::VertexIndexInit(ID3D11Device * pd3dDevice)
{
	Vertex v[Particle_number]; 
    for(int i = 0; i < Particle_number; i++)
	{
		v[i] = Vertex(0, 0, 0, 0.0f, 0.0f);
	}
	/*Vertex v[Particle_number];
	for(int i = 0; i < Particle_number / 100 ; i++)
		for(int j = 0; j < 10; j++)
			for(int p = 0; p < 10; p++)
			{
				int cc =  i * 100 + j * 10 + p;
				v[cc] = Vertex(p, i, j, 0.0f, 0.0f);
			}*/
	/*// Fill in the front face vertex data.
	// Floor: Observe we tile texture coordinates.
	v[0] = Vertex(-3.5f, 0.0f, -10.0f, 0.0f, 1.0f, 0.0f, 0.0f, 4.0f);
	v[1] = Vertex(-3.5f, 0.0f,   0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f);
	v[2] = Vertex( 7.5f, 0.0f,   0.0f, 0.0f, 1.0f, 0.0f, 4.0f, 0.0f);

	v[3] = Vertex(-3.5f, 0.0f, -10.0f, 0.0f, 1.0f, 0.0f, 0.0f, 4.0f);
	v[4] = Vertex( 7.5f, 0.0f,   0.0f, 0.0f, 1.0f, 0.0f, 4.0f, 0.0f);
	v[5] = Vertex( 7.5f, 0.0f, -10.0f, 0.0f, 1.0f, 0.0f, 4.0f, 4.0f);

	// Wall: Observe we tile texture coordinates, and that we
	// leave a gap in the middle for the mirror.
	v[6]  = Vertex(-3.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 2.0f);
	v[7]  = Vertex(-3.5f, 4.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f);
	v[8]  = Vertex(-2.5f, 4.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.5f, 0.0f);

	v[9]  = Vertex(-3.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 2.0f);
	v[10] = Vertex(-2.5f, 4.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.5f, 0.0f);
	v[11] = Vertex(-2.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.5f, 2.0f);

	v[12] = Vertex(2.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 2.0f);
	v[13] = Vertex(2.5f, 4.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f);
	v[14] = Vertex(7.5f, 4.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, 0.0f);

	v[15] = Vertex(2.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 2.0f);
	v[16] = Vertex(7.5f, 4.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, 0.0f);
	v[17] = Vertex(7.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 2.0f, 2.0f);

	v[18] = Vertex(-3.5f, 4.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 1.0f);
	v[19] = Vertex(-3.5f, 6.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f);
	v[20] = Vertex( 7.5f, 6.0f, 0.0f, 0.0f, 0.0f, -1.0f, 6.0f, 0.0f);

	v[21] = Vertex(-3.5f, 4.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 1.0f);
	v[22] = Vertex( 7.5f, 6.0f, 0.0f, 0.0f, 0.0f, -1.0f, 6.0f, 0.0f);
	v[23] = Vertex( 7.5f, 4.0f, 0.0f, 0.0f, 0.0f, -1.0f, 6.0f, 1.0f);

	// Mirror
	v[24] = Vertex(-2.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 1.0f);
	v[25] = Vertex(-2.5f, 4.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f);
	v[26] = Vertex( 2.5f, 4.0f, 0.0f, 0.0f, 0.0f, -1.0f, 1.0f, 0.0f);

	v[27] = Vertex(-2.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 1.0f);
	v[28] = Vertex( 2.5f, 4.0f, 0.0f, 0.0f, 0.0f, -1.0f, 1.0f, 0.0f);
	v[29] = Vertex( 2.5f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 1.0f, 1.0f);*/


	D3D11_BUFFER_DESC mVertexBufferDesc;
	mVertexBufferDesc.Usage = D3D11_USAGE_IMMUTABLE;
	mVertexBufferDesc.MiscFlags = 0;
	mVertexBufferDesc.CPUAccessFlags = 0;
	mVertexBufferDesc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	mVertexBufferDesc.ByteWidth = sizeof(Vertex) * Particle_number;
	mVertexBufferDesc.StructureByteStride = 0;


	D3D11_SUBRESOURCE_DATA mVertexSubReData;
	mVertexSubReData.pSysMem = v;

	pd3dDevice->CreateBuffer(&mVertexBufferDesc,&mVertexSubReData,&m_VB);



	/*std::ifstream fin("Models/skull.txt");

	if(!fin)
	{
		MessageBox(0, L"Models/skull.txt not found.", 0, 0);
		return;
	}

	UINT vcount = 0;
	UINT tcount = 0;
	std::string ignore;

	fin >> ignore >> vcount;
	fin >> ignore >> tcount;
	fin >> ignore >> ignore >> ignore >> ignore;


	std::vector<Vertex> vertices(vcount);
	for(UINT i = 0; i < vcount; ++i)
	{
		fin >> vertices[i].Pos.x >> vertices[i].Pos.y >> vertices[i].Pos.z;

		//vertices[i].Color = black;

		// Normal not used in this demo.
		fin >> vertices[i].Normal.x >> vertices[i].Normal.y >> vertices[i].Normal.z;
	}

	fin >> ignore;
	fin >> ignore;
	fin >> ignore;

	m_iIndexCount = 3 * tcount;
	std::vector<UINT> indices(m_iIndexCount);
	for(UINT i = 0; i < tcount; ++i)
	{
		fin >> indices[i*3+0] >> indices[i*3+1] >> indices[i*3+2];
	}

	fin.close();

	D3D11_BUFFER_DESC mVertexBufferDesc_Skull;
	mVertexBufferDesc_Skull.Usage = D3D11_USAGE_IMMUTABLE;
	mVertexBufferDesc_Skull.MiscFlags = 0;
	mVertexBufferDesc_Skull.CPUAccessFlags = 0;
	mVertexBufferDesc_Skull.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	mVertexBufferDesc_Skull.ByteWidth = sizeof(Vertex) * vcount;
	mVertexBufferDesc_Skull.StructureByteStride = 0;


	D3D11_SUBRESOURCE_DATA mVertexSubReData_Skull;
	mVertexSubReData_Skull.pSysMem = &vertices[0];

	pd3dDevice->CreateBuffer(&mVertexBufferDesc_Skull,&mVertexSubReData_Skull,&m_VBSkull);


	D3D11_BUFFER_DESC mIndexBufferDesc;
	mIndexBufferDesc.Usage =  D3D11_USAGE_IMMUTABLE;
	mIndexBufferDesc.BindFlags = D3D11_BIND_INDEX_BUFFER;
	mIndexBufferDesc.ByteWidth = sizeof(UINT) * m_iIndexCount;
	mIndexBufferDesc.CPUAccessFlags = 0;
	mIndexBufferDesc.MiscFlags = 0;
	mIndexBufferDesc.StructureByteStride = 0;

	D3D11_SUBRESOURCE_DATA mIndexSubReData;
	mIndexSubReData.pSysMem = &indices[0];


	pd3dDevice->CreateBuffer(&mIndexBufferDesc,&mIndexSubReData,&m_IB);*/
}

void VertexIndexControl::Release()
{
	RELEASE_COM(m_VB);
	RELEASE_COM(m_IB);
	RELEASE_COM(m_VBSkull);
	RELEASE_COM(m_InputLayout);
}

HRESULT VertexIndexControl::SetLayout(ID3D11Device * pd3dDevice)
{
	//ÉèÖÃLayout
	D3D11_INPUT_ELEMENT_DESC d3dInputDesc[] =
	{
		{"POSITION",0,DXGI_FORMAT_R32G32B32_FLOAT,0,0,D3D11_INPUT_PER_VERTEX_DATA,0},
		{"TEXCOORD",0,DXGI_FORMAT_R32G32_FLOAT,0,24,D3D11_INPUT_PER_VERTEX_DATA,0},
		{"NORMAL",0,DXGI_FORMAT_R32G32B32_FLOAT,0,12,D3D11_INPUT_PER_VERTEX_DATA,0},


	};

	D3DX11_PASS_DESC passDesc;
	FXControl::m_fxTech->GetPassByIndex(0)->GetDesc(&passDesc);


	HRESULT hr = pd3dDevice->CreateInputLayout(d3dInputDesc,3,passDesc.pIAInputSignature,passDesc.IAInputSignatureSize,&m_InputLayout);

	if(FAILED(hr))
	{
		MessageBoxA(0,"Fail to Create Input Layout",0,0);
		return E_FAIL;
	}

	return S_OK;
}



