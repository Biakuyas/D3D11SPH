cbuffer cbObject{
	float4x4 matWorldViewProj;
	float4x4 matWorld;
	float4x4 matWorldInverseTranspose;
	float3   EyePos;
	float3   LightDir;
};

Texture2D gTexture;

SamplerState samAnisotropic{
	Filter = ANISOTROPIC;
	MaxAnisotropy = 4;

	AddressU = WRAP;
	AddressV = WRAP;
};

struct VertexIn{
	float3 Pos : POSITION;
	float3 Normal : NORMAL;
	float2 Tex : TEXCOORD;
};

struct VertexOut{
	float4 Pos : SV_POSITION;
	float3 PosW : POSITION;
	float3 Normal: NORMAL;
	float2 Tex : TEXCOORD;
};

VertexOut VS(VertexIn vin){
	VertexOut vout;

	vout.Pos = mul(float4(vin.Pos,1.0),matWorldViewProj);
	vout.PosW = mul(float4(vin.Pos,1.0),matWorld).xyz;
	vout.Normal = normalize(mul(vin.Normal,(float3x3)matWorldInverseTranspose));
    vout.Tex =  vin.Tex;
	//vout.Color = vin.Color;

	return vout;
}

float4 PS(VertexOut pin, uniform bool istransparent,uniform bool isshadow) : SV_Target{

	//默认各种光强都为1.0f

/*	float4 DiffuseMaterial = float4(0.5f,0.5f,0.5f,1.0f);
	float4 SpecularMaterial = float4(0.5f,0.5f,0.5f,1.0f);
	float4 AmbientMaterial = float4(0.5f,0.5f,0.5f,1.0f);


	float3 PostoEye = normalize(EyePos - pin.PosW);

	float AmbientPower = 0.3f;

	float PrepareDiffuse = dot(normalize(-LightDir),pin.Normal);
	float DiffusePower = clamp(PrepareDiffuse,0,1);

	float SpecularPower = 0;

	if(PrepareDiffuse > 0)
	SpecularPower = pow(max(dot(reflect(normalize(LightDir),pin.Normal),PostoEye),0),10);

	float4 TexColor;
		TexColor = gTexture.Sample(samAnisotropic,pin.Tex);


	float4 Color = TexColor * (AmbientPower * AmbientMaterial +  SpecularPower * SpecularMaterial  + DiffusePower * DiffuseMaterial);

	if(istransparent)
		Color.w = 0.5f;

    if(isshadow)
			Color = float4(0,0,0,0.5f);*/

	float4 Color = float4(10.0f,10.0f,10.0f,5.0f);

	return Color;

}
 
technique11  BasicTech{

	pass P0{
		SetVertexShader( CompileShader( vs_5_0, VS() ) );
		SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_5_0, PS(false,false) ) );
	}

}

technique11 TransparentTech{

	pass P0{
		SetVertexShader( CompileShader( vs_5_0, VS() ) );
		SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_5_0, PS(true,false) ) );
	}

}

technique11 ShadowTech{

	pass P0{
		SetVertexShader( CompileShader( vs_5_0, VS() ) );
		SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_5_0, PS(false,true) ) );
	}

}


/*cbuffer cbPerObject
{
	float4x4 gWorldViewProj; 
};

struct VertexIn
{
	float3 PosL  : POSITION;
    float4 Color : COLOR;
};

struct VertexOut
{
	float4 PosH  : SV_POSITION;
    float4 Color : COLOR;
};

VertexOut VS(VertexIn vin)
{
	VertexOut vout;
	
	// Transform to homogeneous clip space.
	vout.PosH = mul(float4(vin.PosL, 1.0f), gWorldViewProj);
	
	// Just pass vertex color into the pixel shader.
    vout.Color = vin.Color;
    
    return vout;
}

float4 PS(VertexOut pin) : SV_Target
{
    return pin.Color;
}

technique11 ColorTech
{
    pass P0
    {
        SetVertexShader( CompileShader( vs_5_0, VS() ) );
		SetGeometryShader( NULL );
        SetPixelShader( CompileShader( ps_5_0, PS() ) );
    }
}*/