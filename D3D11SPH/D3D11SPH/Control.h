#ifndef CONTROL_H_
#define CONTROL_H_
#include "SPH_Math.h"



/*
1. 加入UI //1，2 2/3
2. 加入一些碰撞体 //3
3. 多做几个例子进行演示 // 3
4. 把粒子模式下的点换成公告牌或球体//4 Clear
5. 考虑要不要把渲染部分也加进来
6. 程序优化//4
7. GPGPU // 待定
8. 
//可能有内存泄露

*/


struct SPHBox
{
	float min_x,min_y,min_z;
	float max_x,max_y,max_z;

	SPHBox(float minx,float miny,float minz,float maxx,float maxy, float maxz)
	{
		min_x = minx;
		min_y = miny;
		min_z = minz;
		max_x = maxx;
		max_y = maxy;
		max_z = maxz;
	}
	SPHBox()
	{
		min_x = 0;
		min_y = 0;
		min_z = 0;
		max_x = 0;
		max_y = 0;
		max_z = 0;
	}

};

struct SPHVertex
{
	XMFLOAT3 pos;
	float val;

	SPHVertex()
	{
		pos = XMFLOAT3(0,0,0);
		val = 0;
	}
};

struct SPHCube
{
	XMFLOAT3 pos[8];
	float val[8];
	XMFLOAT3 normal[8];


	SPHCube()
	{
		for(int i = 0; i < 8; i++)
		{
		  pos[i] = XMFLOAT3(0,0,0);
		  val[i] = 0;
		  normal[i] = XMFLOAT3(0,0,0);
		}
	}
};

struct SPHTriangle
{
	XMFLOAT3 v1,v2,v3;
	XMFLOAT3 n1,n2,n3;
	SPHTriangle(XMFLOAT3 v_1,XMFLOAT3 v_2,XMFLOAT3 v_3,
		XMFLOAT3 n_1,XMFLOAT3 n_2,XMFLOAT3 n_3)
	{
		v1 = v_1;
		v2 = v_2;
		v3 = v_3;
		n1 = n_1;
		n2 = n_2;
		n3 = n_3;
	}
};
#endif