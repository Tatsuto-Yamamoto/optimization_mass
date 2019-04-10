# main
import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

import particle_mass as pmas #自作の、粒子の質量が入ったライブラリ
import gd #自作の、勾配降下法を行うライブラリ

# 以下の値が、"いい感じ" な値らしい
# x = np.array([5,-14.48,8.6,0.28,-280000,34.9,17])

# 質量の実験値	
mex_top = 173000
mex_up = 2.2

def f(x): #評価関数を返す
	f = (1 - pmas.ups(x))**2
	return f

def df(x): #fの勾配
	const = -2 * (1 - pmas.ups(x))
	df = const * pmas.dups(x)
	# df = np.array(df)
	# print(df)
	return df

########
initial = np.array([5.0,-14.48,8.6,0.28,-280000.0,34.9,17.0]) #初期値	
########

algo = gd.GradientDescent(f,df) # 勾配降下法のインスタンスを生成
algo.solve(initial) # 勾配降下法を計算

print(algo.x__) # 計算結果で得られた xx を出力
print(algo.opt__) # 最小化された目的関数 f を出力

# plt.scatter(initial[0],initial[1])
# plt.plot(algo.path__[:,0], algo.path__[:,1], linewidth = 1.5)

# xs = np.linspace(-2,2,300)
# ys = np.linspace(-2,2,300)
# xmesh, ymesh = np.meshgrid(xs,ys)
# xx = np.r_[xmesh.reshape(1,-1), ymesh.reshape(1,-1)]
# levels = [-3.0,-2.5,-2.0,0,1,2]
# # plt.contour(xs,ys,f(xx).reshape(xmesh,ymesh),levels = levels, colors = "k", linestyles = "dotted")

# plt.show()

